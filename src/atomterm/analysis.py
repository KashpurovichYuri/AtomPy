"""L-multiplicity (term-symbol) analysis of many-electron states.

Input: full matrices in the determinant basis (atomterm.operators) for the
Hamiltonian H, Lz, L^2, S^2. The states are the eigenstates of H; each is
classified by

    <L^2> = L(L+1),  <S^2> = S(S+1)  =>  term symbol  ^(2S+1)L

and the Lz eigenvalue spectrum reveals the M_L content of the space.
Symmetry preservation is quantified by the Frobenius norms of the
commutators [H, Lz] and [H, L^2]: exact rotational invariance means both
vanish (up to numerics / basis-set closure of the active space).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
import json

import numpy as np

from atomterm.detbasis import commutator, frobenius_norm

__all__ = [
    "quantum_number",
    "spin_quantum_number",
    "term_symbol",
    "AnalysisResult",
    "analyze",
    "configuration_summary",
    "format_report",
    "save_report",
]

_L_LETTERS = "SPDFGHIKLMNOP"  # J skipped by spectroscopic convention;
# length 13 covers L <= 12: f^3 d^1 reaches L = 11, f^4 reaches L = 12
# (both manifolds live in the uranium window)


def quantum_number(observed: float) -> int:
    """Quantum number N with N(N+1) closest to ``observed`` (integer grid)."""
    n_real = (np.sqrt(1.0 + 4.0 * observed) - 1.0) / 2.0
    return int(round(n_real))


def spin_quantum_number(s2: float) -> float:
    """Spin S with S(S+1) closest to ``s2``, on the HALF-integer grid.

    Orbital L is always an integer, but S of an odd-electron state is a
    half-integer (S^2 = 0.75 -> S = 1/2 -> multiplicity 2): rounding S to
    the integer grid silently relabels doublets/quartets as singlets/
    triplets (first exposed by Sc, the project's first doublet).
    """
    s_real = (np.sqrt(1.0 + 4.0 * s2) - 1.0) / 2.0
    return round(2.0 * s_real) / 2.0


def term_symbol(s2: float, l2: float) -> str:
    """Term symbol ^(2S+1)L from expectation values of S^2 and L^2."""
    s = spin_quantum_number(s2)
    l = quantum_number(l2)
    return f"{int(round(2 * s + 1))}{_L_LETTERS[l]}"


@dataclass
class StateSummary:
    index: int
    energy: float
    l2: float
    l: int
    s2: float
    s: float  # half-integer allowed (odd-electron states)
    term: str


@dataclass
class AnalysisResult:
    sector: str  # e.g. "Ms=1 (na=2, nb=0)"
    n_det: int
    states: list[StateSummary]
    lz_spectrum: list[float]
    comm_norm_h_lz: float
    comm_norm_h_l2: float
    h_norm: float


def _expectation(evecs: np.ndarray, op: np.ndarray) -> np.ndarray:
    """Diagonal of evecs^dag op evecs (i.e. <state|op|state>)."""
    mat = evecs.conj().T @ op @ evecs
    return np.real(np.diag(mat))


def analyze(h: np.ndarray, lz: np.ndarray, l2: np.ndarray, s2: np.ndarray, sector: str) -> AnalysisResult:
    """Classify the eigenstates of H and test the [H, L] commutators."""
    if frobenius_norm(h.imag) > 1e-9 * max(frobenius_norm(h.real), 1.0):
        raise ValueError("Hamiltonian matrix should be real")
    h = np.real(h)
    evals, evecs = np.linalg.eigh(h)

    l2_exp = _expectation(evecs, l2)
    s2_exp = _expectation(evecs, s2)
    states = [
        StateSummary(
            index=i,
            energy=float(evals[i]),
            l2=float(l2_exp[i]),
            l=quantum_number(l2_exp[i]),
            s2=float(s2_exp[i]),
            s=spin_quantum_number(s2_exp[i]),
            term=term_symbol(s2_exp[i], l2_exp[i]),
        )
        for i in range(len(evals))
    ]
    lz_eigs = np.sort(np.real(np.linalg.eigvals(lz)))

    h_norm = max(frobenius_norm(h), 1e-300)
    return AnalysisResult(
        sector=sector,
        n_det=len(evals),
        states=states,
        lz_spectrum=[float(v) for v in lz_eigs],
        comm_norm_h_lz=frobenius_norm(commutator(h, lz)),
        comm_norm_h_l2=frobenius_norm(commutator(h, l2)),
        h_norm=h_norm,
    )


def format_report(res: AnalysisResult, title: str) -> str:
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
    lines.append(f"||[H, Lz]||_F  = {res.comm_norm_h_lz:.3e}   "
                 f"(relative {res.comm_norm_h_lz / res.h_norm:.3e})")
    lines.append(f"||[H, L^2]||_F = {res.comm_norm_h_l2:.3e}   "
                 f"(relative {res.comm_norm_h_l2 / res.h_norm:.3e})")
    return "\n".join(lines)


def save_report(res: AnalysisResult, path: str) -> None:
    payload = asdict(res)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)


def configuration_summary(mol, s, ci_list, c_win: np.ndarray, na: int, nb: int,
                          top: int = 2) -> list[str]:
    """Approximate valence-configuration label of each CI root in a window.

    Method-agnostic: works for ANY window and CI vectors. Each window
    column's l-character comes from the S-metric projection
    (:func:`atomterm.drivers.shell_weights`); a determinant's configuration
    is the l-mass of its occupied columns; each root is labeled by its
    ``top`` dominant determinants.

    Approximation level: exact when the window columns are raw shell AOs
    (determinants ARE AO occupations, e.g. "p^2.00"); approximate for
    optimized windows, whose columns are delocalized over radials/l -- the
    label then reports the dominant character (enough to tell a 3p^2 state
    from an escaped 3s^2 one). Returns one string per root, e.g.
    ``"p^2.00 (|c|=0.98)"``.
    """
    from atomterm.drivers import shell_weights
    from atomterm.detbasis import determinant_basis

    names, wl = shell_weights(mol, c_win, s)  # (ncas, n_l)
    norb = c_win.shape[1]
    dets = determinant_basis(norb, na, nb)

    labels = []
    for v in ci_list:
        v = np.asarray(v).ravel()
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
