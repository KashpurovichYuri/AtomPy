"""Analysis: term classification of determinant-basis eigenstates.

PURPOSE: given the operator matrices of
:func:`operators.build_operators`, classify every eigenstate of H by
its exact content -- L, S, and the term symbol ``^(2S+1)L`` from
``<L^2> = L(L+1)``, ``<S^2> = S(S+1)`` -- and quantify how exactly
the CI space preserves rotational symmetry.

Two meters: the Frobenius norms of ``[H, Lz]`` and ``[H, L^2]`` --
exact rotational invariance means both vanish (up to numerics /
basis-set closure of the active space). The Lz eigenvalue spectrum
reports the M_L content of the space.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
import json

import numpy as np

from atomterm.detbasis import commutator, determinant_basis, frobenius_norm
from atomterm.windows import lm_sectors, shell_weights

__all__ = [
    "quantum_number",
    "spin_quantum_number",
    "term_symbol",
    "StateSummary",
    "AnalysisResult",
    "analyze",
    "configuration_summary",
    "format_report",
    "save_report",
]

_L_LETTERS = "SPDFGHIKLMNOQRTUVWXYZ"  # J skipped, and the second P and
# second S skipped, by spectroscopic convention (S P D F G H I K L M N
# O Q R T ...); length 21 covers L <= 20: f^3 d^1 reaches L = 8, f^4
# terms reach L = 12 (Q) -- both manifolds live in the uranium window


# --------------------------------------------------------- the term labels


def quantum_number(observed: float) -> int:
    """Quantum number N with N(N+1) closest to ``observed`` (integer
    grid).

    Examples
    --------
    >>> quantum_number(6.0)
    2
    """
    n_real = (np.sqrt(1.0 + 4.0 * observed) - 1.0) / 2.0
    return int(round(n_real))


def spin_quantum_number(s2: float) -> float:
    """Spin S with S(S+1) closest to ``s2``, on the HALF-integer grid.

    Orbital L is always an integer, but S of an odd-electron state is
    a half-integer (S^2 = 0.75 -> S = 1/2 -> multiplicity 2): rounding
    S to the integer grid silently relabels doublets/quartets as
    singlets/triplets (first exposed by Sc, the project's first
    doublet).
    """
    s_real = (np.sqrt(1.0 + 4.0 * s2) - 1.0) / 2.0
    return round(2.0 * s_real) / 2.0


def term_symbol(s2: float, l2: float) -> str:
    """Term symbol ``^(2S+1)L`` from <S^2> and <L^2> expectation
    values.

    Raises
    ------
    ValueError
        If L exceeds the letter table (L <= 20).

    Examples
    --------
    >>> term_symbol(2.0, 2.0)
    '3P'
    """
    s = spin_quantum_number(s2)
    l = quantum_number(l2)
    if l >= len(_L_LETTERS):
        raise ValueError(f"L = {l} beyond the letter table "
                         f"({_L_LETTERS})")
    return f"{int(round(2 * s + 1))}{_L_LETTERS[l]}"


# -------------------------------------------------------------- the records


@dataclass
class StateSummary:
    """One eigenstate of H with its labels.

    Attributes
    ----------
    index : position in the ascending-energy list.
    energy : E in Eh.
    l2 : ``<L^2>``.
    l : the snapped L.
    s2 : ``<S^2>``.
    s : the snapped S -- half-integer allowed (odd-electron states).
    term : the term symbol string.
    lz : ``<L_z>`` of the root (the M_L label).
    lz2 : ``<L_z^2>`` of the root: M_L^2 on broken folds, gauge on
        isotropic ones (README, the <Lz^2> note).
    """

    index: int
    energy: float
    l2: float
    l: int
    s2: float
    s: float
    term: str
    lz: float = 0.0
    lz2: float = 0.0


@dataclass
class AnalysisResult:
    """The full result of one sector analysis.

    Attributes
    ----------
    sector : the label, e.g. ``"Ms=+1 (na=2, nb=0)"``.
    n_det : determinants in the sector.
    states : one StateSummary per eigenstate, ascending energy.
    lz_spectrum : sorted Lz eigenvalues (the M_L content).
    comm_norm_h_lz, comm_norm_h_l2 : the two symmetry meters
        (Frobenius norms of the commutators).
    h_norm : ``||H||_F`` -- the scale the meters are read against.
    sz_spectrum : sorted Sz eigenvalues, exactly ``M_S`` times the
        identity -- reported for the record (H is spin-free,
        ``[H, Sz] = 0`` identically), or None.
    """

    sector: str
    n_det: int
    states: list[StateSummary]
    lz_spectrum: list[float]
    comm_norm_h_lz: float
    comm_norm_h_l2: float
    h_norm: float
    sz_spectrum: list[float] | None = None


# -------------------------------------------------------------- the analysis


def _expectation(evecs: np.ndarray, op: np.ndarray) -> np.ndarray:
    """Diagonal of evecs^dag op evecs (i.e. <state|op|state>)."""
    mat = evecs.conj().T @ op @ evecs
    return np.real(np.diag(mat))


def analyze(h: np.ndarray, lz: np.ndarray, l2: np.ndarray, s2: np.ndarray, sector: str,
            sz: np.ndarray | None = None, lz2: np.ndarray | None = None) -> AnalysisResult:
    """Classify the eigenstates of H and measure the [H, L] meters.

    Parameters
    ----------
    h, lz, l2, s2 : the full sector matrices from
        ``operators.build_operators``.
    sector : the label stored on the result.
    sz : optional -- the delivered Sz matrix; its spectrum (exactly
        the sector's M_S) is reported for the record, certifying
        nothing.
    lz2 : optional -- the delivered Lz^2 matrix; its per-root
        expectation is recorded on every state.

    Returns
    -------
    AnalysisResult

    Raises
    ------
    ValueError
        If H carries non-noise imaginary parts (a construction error
        upstream, not physics to diagonalize).
    """
    if frobenius_norm(h.imag) > 1e-9 * max(frobenius_norm(h.real), 1.0):
        raise ValueError("Hamiltonian matrix should be real")
    h = np.real(h)
    evals, evecs = np.linalg.eigh(h)

    l2_exp = _expectation(evecs, l2)
    s2_exp = _expectation(evecs, s2)
    lz_exp = _expectation(evecs, lz)
    lz2_exp = _expectation(evecs, lz2) if lz2 is not None else None
    states = [
        StateSummary(
            index=i,
            energy=float(evals[i]),
            l2=float(l2_exp[i]),
            l=quantum_number(l2_exp[i]),
            s2=float(s2_exp[i]),
            s=spin_quantum_number(s2_exp[i]),
            term=term_symbol(s2_exp[i], l2_exp[i]),
            lz=float(lz_exp[i]),
            lz2=(float(lz2_exp[i]) if lz2_exp is not None else 0.0),
        )
        for i in range(len(evals))
    ]
    lz_eigs = np.sort(np.real(np.linalg.eigvals(lz)))
    sz_eigs = (np.sort(np.real(np.linalg.eigvals(sz)))
               if sz is not None else None)

    h_norm = max(frobenius_norm(h), 1e-300)
    return AnalysisResult(
        sector=sector,
        n_det=len(evals),
        states=states,
        lz_spectrum=[float(v) for v in lz_eigs],
        sz_spectrum=([float(v) for v in sz_eigs] if sz_eigs is not None
                     else None),
        comm_norm_h_lz=frobenius_norm(commutator(h, lz)),
        comm_norm_h_l2=frobenius_norm(commutator(h, l2)),
        h_norm=h_norm,
    )


# --------------------------------------------------------------- the reports


def format_report(res: AnalysisResult, title: str) -> str:
    """Human-readable report: the state table, the M_L / M_S content,
    and the two commutator norms (absolute and relative)."""
    lines = [title, "-" * len(title)]
    lines.append(f"sector: {res.sector}   ({res.n_det} determinants)")
    lines.append("")
    lines.append(f"{'#':>2} {'E (Eh)':>16} {'<L^2>':>10} {'L':>2} {'<S^2>':>10} "
                  f"{'S':>4} {'term':>5}")
    for st in res.states:
        lines.append(
            f"{st.index:2d} {st.energy:16.8f} {st.l2:10.6f} {st.l:2d} "
            f"{st.s2:10.6f} {st.s:4g} {st.term:>5}"
        )
    lines.append("")
    lines.append(f"Lz eigenvalues (M_L content): "
                 f"{np.round(res.lz_spectrum, 6).tolist()}")
    if res.sz_spectrum is not None:
        uniq, counts = np.unique(np.round(res.sz_spectrum, 6),
                                 return_counts=True)
        parts = []
        for v, c in zip(uniq, counts):
            txt = "0" if v == 0 else f"{v:+g}"
            parts.append(f"[{txt}] x{c}" if c > 1 else f"[{txt}]")
        lines.append("Sz eigenvalues (M_S content): " + ", ".join(parts)
                     + "   (exactly M_S·1 in a fixed sector)")
    lines.append(f"||[H, Lz]||_F  = {res.comm_norm_h_lz:.3e}   "
                 f"(relative {res.comm_norm_h_lz / res.h_norm:.3e})")
    lines.append(f"||[H, L^2]||_F = {res.comm_norm_h_l2:.3e}   "
                 f"(relative {res.comm_norm_h_l2 / res.h_norm:.3e})")
    return "\n".join(lines)


def save_report(res: AnalysisResult, path: str) -> None:
    """The full AnalysisResult as JSON (asdict dump)."""
    payload = asdict(res)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)


# ------------------------------------------------- the configuration readers


def configuration_summary(mol, s, ci_list, c_win: np.ndarray, na: int, nb: int,
                          top: int = 2) -> list[str]:
    """Approximate valence-configuration label of each CI root.

    Works for ANY window and CI vectors. Each window column's
    l-character comes from the S-metric projection
    (:func:`atomterm.windows.shell_weights`); a determinant's
    configuration is the l-mass of its occupied columns; each root is
    labeled by its ``top`` dominant determinants.

    Parameters
    ----------
    ci_list : CI vectors in PYSCF order -- the ordering bridge runs
        inside (fact 24: PySCF's FCI string order differs from
        detbasis' from ``norb >= 4``).
    top : determinants reported per root.

    Returns
    -------
    list[str], one per root, e.g. ``"p^2.00 (|c|=0.98)"``.

    Notes
    -----
    Exact for raw shell-AO windows (determinants ARE AO occupations);
    approximate for optimized windows, whose columns are delocalized
    over radials/l -- the label then reports the dominant character
    (enough to tell a 3p^2 state from an escaped 3s^2 one).
    """
    from atomterm.cells import ci_in_detbasis_order  # lazy on purpose:
    # cells -> stages -> analysis is a module-level cycle

    names, wl = shell_weights(mol, c_win, s)  # (ncas, n_l)
    norb = c_win.shape[1]
    dets = determinant_basis(norb, na, nb)

    labels = []
    for v in ci_list:
        # pyscf FCI ordering -> detbasis ordering (they differ for
        # norb >= 4 with >= 2 electrons in a block; the bridge certifies)
        v = ci_in_detbasis_order(v, norb, na, nb)
        order = np.argsort(-np.abs(v))[:top]
        parts = []
        for i in order:
            occ_cols = sorted(d % norb for d in dets[i])
            mass = {n: 0.0 for n in names}
            for col in occ_cols:
                for k, n in enumerate(names):
                    mass[n] += wl[col, k]
            label = "".join(
                f"{n}^{mass[n]:.2f}" for n in names if mass[n] > 0.01
            )
            parts.append(f"{label} (|c|={abs(v[i]):.2f})")
        labels.append("; ".join(parts))
    return labels


def shell_occupations(mol, s, c_win: np.ndarray, ci_ground, na: int,
                      nb: int) -> dict:
    """Electron occupations per angular (l) AO subspace, ground root
    -- ``{'s': .., 'p': .., ...}``.

    The ELECTRON scale of the window analysis (the span counterpart
    is :func:`~atomterm.windows.shell_weights`): gamma = the ground
    root's spin-summed 1-RDM in the window-MO basis of ``c_win``
    (S-orthonormal columns); then ``n_l = Tr[gamma C^dag S_l C]``
    with ``S_l`` the l-block of the AO overlap -- the expectation of
    the l-subspace number operator.

    Parameters
    ----------
    ci_ground : the ground root's CI vector in PYSCF order (the
        ``cas.ci`` / cell-store convention); mapped through the
        ordering bridge (fact 24).

    Notes
    -----
    Different-l AO blocks are exactly S-orthogonal (one center), so
    ``sum_l S_l = S`` and the ``n_l`` partition ``na + nb`` EXACTLY.
    The bare AO diagonal of ``C gamma C^dag`` is NOT it: same-l
    shells are non-orthogonal and overshoot by the cross terms.
    Term-stable inside a degenerate cluster (same configuration
    content).
    """
    from atomterm.natural import ci_rdm1
    from atomterm.cells import ci_in_detbasis_order  # lazy: the
    # cells -> stages -> analysis module cycle

    norb = c_win.shape[1]
    dets = determinant_basis(norb, na, nb)
    v = ci_in_detbasis_order(np.asarray(ci_ground).ravel(), norb, na, nb)
    gamma = ci_rdm1(dets, v, norb)
    g_ao = c_win @ gamma @ c_win.conj().T          # the AO-basis 1-RDM
    ds = g_ao @ s                                   # Mulliken rows
    # l per AO row from the shell walker (lm_sectors) -- NO label
    # parsing (fact 1: labels carry member suffixes; the walker is the
    # fused-shell-honest source)
    l_of_row = {row: l
                for (l, m), rows in lm_sectors(mol).items() for row in rows}
    letters = "spdfgh"
    out: dict[str, float] = {}
    for i in range(mol.nao_nr()):
        key = letters[l_of_row[i]]
        out[key] = out.get(key, 0.0) + float(np.real(ds[i, i]))
    return out
