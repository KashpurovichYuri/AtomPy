"""Windows: construction, selection, and the window-kind registry.

A *window* is the active space. Three construction kinds, one guarantee
each: contracted shells and primitive bands (AO rows -- m-complete by
construction), raw-SCF MO spans (never rotation-closed), and
natural-orbital occupation groups (exact from a term-averaged density;
`atomterm.natural` builds those from a parent CASCI). A
:class:`Window` object carries the construction result and materializes
to a full MO substrate through :meth:`Window.columns`; the
bounded-slice contract (``ncas`` from ``mc.ncas``; active columns at
``[ncore : ncore + ncas]``) is owned here.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from atomterm.natural import (ci_rdm1, natural_orbitals,
                              select_no_window)

__all__ = [
    "Window",
    "occupied_core_idx",
    "select_shell_ao",
    "select_valence_ao",
    "primitive_chain",
    "select_primitives",
    "shell_weights",
    "avas_window",
    "span_overlap",
    "complete_orthonormal",
    "build_full_mo",
    "lm_sectors",
    "radial_basis",
    "radial_mask",
    "window_shell",
    "window_primitive",
    "window_scf_mo",
    "window_no",
    "SELECTORS",
    "lowdin",
    "lowdin",
]


def ao_span_projector(s_mat: np.ndarray, rows: list[int]) -> np.ndarray:
    """S-metric projector onto an AO row span: ``P = S[:, o] S_oo^{-1} S[o, :]``.

    The weight of an S-orthonormal MO column ``c`` on the span is
    ``c^dag P c`` -- bounded by 1 (a true projector weight), unlike the
    bare quadratic form ``c_o^dag S_oo c_o`` which is unbounded for
    diffuse bases (fact 15c).
    """
    s_oo = s_mat[np.ix_(rows, rows)]
    return s_mat[:, rows] @ np.linalg.inv(s_oo) @ s_mat[rows, :]


def lowdin(c: np.ndarray, s: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """S-metric Loewdin orthonormalization of the columns of ``c``.

    Raises ValueError when the column overlap is near-singular (the
    needle class of fact 16e: a nearly-linearly-dependent admission)
    instead of silently NaN-ing everything downstream -- pad or drop
    the needle directions, or raise ``tol``.
    """
    ov = c.T @ s @ c
    ev, u = np.linalg.eigh(ov)
    if ev.min() <= tol:
        raise ValueError(
            f"Lowdin on a near-linearly-dependent block: min overlap "
            f"eigenvalue {ev.min():.2e} <= tol {tol:.0e}; pad/drop the "
            f"needle directions (fact 16e) or raise tol")
    return c @ (u / np.sqrt(ev)) @ u.T


def occupied_core_idx(mf, ncore: int) -> list[int]:
    """Indices of the ``ncore`` doubly-occupied MOs -- the honest core
    block for every build_full_mo/radial_basis call.

    Canonical MO ORDER is meaningless on near-degenerate uncontracted
    bases (U: the first canonical MOs are virtual junk at eps ~ 0.02,
    while the closed block sits at 22/61/88/120/143); ``range(ncore)``
    then silently folds five garbage columns (measured: max|dE| = 14.5
    Eh). On contracted bases the two coincide -- this makes the
    coincidence a checked invariant instead of an assumption.
    """
    closed = np.flatnonzero(np.asarray(mf.mo_occ) == 2.0)
    if len(closed) < ncore:
        raise ValueError(
            f"asked for {ncore} core orbitals; ROHF closed block has "
            f"{len(closed)}")
    if not np.array_equal(closed[:ncore], np.arange(ncore)):
        import warnings

        warnings.warn(
            f"core block is NOT the first {ncore} canonical MOs "
            f"(closed at {closed[:ncore].tolist()}); folding the "
            f"occupied block instead", stacklevel=2)
    return closed[:ncore].tolist()


def span_rows(mol, pattern: str) -> list[int]:
    """AO rows of one occupation-span pattern: ``"l0".."l9"`` selects a
    WHOLE angular-momentum subspace (the honest vocabulary for uncontracted
    bases, where labels are shell counters -- fact 1); anything else is a
    label prefix (:func:`select_shell_ao`)."""
    if len(pattern) == 2 and pattern[0] == "l" and pattern[1].isdigit():
        l = int(pattern[1])
        return sorted(r for (ll, m), rws in lm_sectors(mol).items()
                      if ll == l for r in rws)
    rows = select_shell_ao(mol, pattern)
    if not rows:
        raise ValueError(
            f"span pattern {pattern!r} matched no AO labels; available: "
            f"{[line.split()[-1] for line in mol.ao_labels()]}")
    return rows


def select_shell_ao(mol, pattern: str) -> list[int]:
    """AO indices of a named valence shell, e.g. pattern="2p" -> [3, 4, 5].

    Plain prefix matching on PySCF's AO labels (``search_ao_label`` patterns
    are regex-flavored and over-match: "2p*" also hits "2s").
    """
    labels = [line.split()[-1] for line in mol.ao_labels()]
    return [i for i, lbl in enumerate(labels) if lbl.startswith(pattern)]


def select_valence_ao(mol, patterns: list[str]) -> list[int]:
    """AO indices of several valence shells, concatenated in pattern order,
    e.g. ``["3d", "4s"]`` -> the five 3d AOs followed by the 4s AO.

    The concatenation defines the active space; it is rotation-closed when
    each pattern contributes a complete m-set of one (n, l) shell.
    """
    idx: list[int] = []
    for pattern in patterns:
        matches = select_shell_ao(mol, pattern)
        if not matches:
            raise ValueError(
                f"shell pattern {pattern!r} matched no AO labels; "
                f"available: {[line.split()[-1] for line in mol.ao_labels()]}"
            )
        idx += matches
    return idx


def primitive_chain(mol, l: int) -> list[tuple[float, list[int]]]:
    """The radial primitives of angular momentum ``l`` as ``(exponent,
    ao_rows)`` pairs, sorted by exponent ASCENDING (diffuse end first) --
    the honest radial coordinate of an uncontracted basis, where PySCF's
    AO labels number shells from the TIGHT end starting at the valence
    name ("5f" = the exp-36.8 primitive on U; fact-1 amendment).

    Walks the basis exactly as :func:`lm_sectors` does (``bas_angular``
    checked, ``nctr`` unpacked): a primitive's m-set occupies the rows
    ``loc[ib] + c*(2l+1) + m`` of its shell. A contraction slot is a
    RAW primitive iff its ``bas_ctr_coeff`` column carries exactly one
    exponent; true contractions (several exponents combined into one
    function) have no single exponent and RAISE -- for them use
    :func:`select_valence_ao` (a contracted shell is one pre-converged
    radial) or MO/NO-based window selection.
    """
    loc = mol.ao_loc_nr()
    width = 2 * l + 1
    chain: list[tuple[float, list[int]]] = []
    for ib in range(mol.nbas):
        if abs(mol.bas_angular(ib)) != l:
            continue
        exps = np.asarray(mol.bas_exp(ib))
        coeff = mol.bas_ctr_coeff(ib)  # (nprim, nctr)
        for c in range(coeff.shape[1]):
            nz = np.flatnonzero(np.abs(coeff[:, c]) > 1e-12)
            if len(nz) != 1:
                raise ValueError(
                    f"shell {ib} (l={l}) contraction slot {c} combines "
                    f"{len(nz)} exponents -- exponent-ranked selection is "
                    "undefined for contractions; use select_valence_ao "
                    "(a contracted shell is one pre-converged radial) or "
                    "MO/NO-based window selection"
                )
            rows = [loc[ib] + c * width + m for m in range(width)]
            chain.append((float(exps[nz[0]]), rows))
    chain.sort(key=lambda t: t[0])
    return chain


def select_primitives(mol, l: int, count: int = None,
                      exponent_band: tuple[float, float] = None) -> list[int]:
    """AO rows of radial primitives of angular momentum ``l``, selected
    by EXPONENT on the :func:`primitive_chain` -- m-complete by
    construction (each primitive contributes its whole 2l+1 m-set), so
    the window is exactly rotation-closed.

    Exactly one of ``count`` / ``exponent_band``:
    * ``exponent_band=(lo, hi)``: every primitive with ``lo <= exp <=
      hi``. The scan protocol: bands are CONTIGUOUS sets built from the
      measured open-MO primitive weights, never hand-picked intervals,
      and conclusions require neighboring-band stability (admit one
      more / one fewer primitive on each side -> same ordering).
    * ``count``: the ``count`` most diffuse primitives (convenience
      shorthand for the diffuse end of the chain).
    """
    chain = primitive_chain(mol, l)
    if (count is None) == (exponent_band is None):
        raise ValueError("give exactly one of count / exponent_band")
    if exponent_band is not None:
        lo, hi = exponent_band
        picked = [rows for exp, rows in chain if lo <= exp <= hi]
        if not picked:
            raise ValueError(
                f"no l={l} primitive has an exponent in "
                f"{exponent_band!r}; the chain (ascending): "
                f"{[round(e, 6) for e, _ in chain]}")
    else:
        if count > len(chain):
            raise ValueError(
                f"asked for {count} primitives of l={l}; the chain has "
                f"{len(chain)}")
        picked = [rows for _, rows in chain[:count]]
    return sorted(r for rows in picked for r in rows)


def avas_window(mf, labels: list[str], ncore: int, threshold: float = None):
    """Window chosen by AVAS (pyscf.mcscf.avas, arXiv:1701.07862): an
    independent ecosystem window selector for cross-validating our raw
    shell AOs (verdicts: fact 19, README §7.3).

    AVAS projects the target AO span (``labels``, e.g. ``["Sc 3d",
    "Sc 4s"]``) out of the converged MOs and keeps eigenvectors of the
    projected overlap above ``threshold`` (its default 0.2). Because
    the construction PROJECTS onto the target span, a full
    shell-dimension AVAS window IS the shell span exactly -- closed by
    construction; a threshold that drops members of an m-set gives a
    rotation-DEFECTIVE span (the negative-control case for the closure
    meters).

    Returns ``(mo_full, ncas, nelecas)`` with the AVAS active span at
    ``[ncore : ncore+ncas]`` (ready for run_casci; electron bookkeeping
    stays ours -- pass the sector you want).
    """
    from pyscf.mcscf import avas

    kw = {} if threshold is None else {"threshold": threshold}
    ncas, nelecas, mo = avas.avas(mf, labels, ncore=ncore, **kw)
    nelecas = np.atleast_1d(np.asarray(nelecas)).ravel()
    return mo, int(ncas), nelecas


def span_overlap(c_a: np.ndarray, c_b: np.ndarray, s: np.ndarray) -> float:
    """Smallest principal cosine between two spans in the S metric
    (1.0 = identical spans, 0.0 = orthogonal direction present).

    The identity-change gauge of the ladder decomposition: sv_min of
    ``C_a^dag S C_b`` is 1 exactly when the converged window kept the
    shell span, and drops first when the window identity changes.
    """
    sv = np.linalg.svd(c_a.conj().T @ s @ c_b, compute_uv=False)
    return float(sv.min())


def shell_weights(mol, c: np.ndarray, s: np.ndarray):
    """S-metric projection norms of columns onto angular (l) AO subspaces.

    ``w[i, k] = c[:, i]^dag S_l c[:, i]`` over ALL AOs of angular momentum
    l_k. S-orthonormal columns make every row of ``w`` sum to 1 -- valid
    because different-l subspaces are exactly orthogonal; grouping by
    (n, l) instead would NOT partition (2p and 3p AO subspaces overlap,
    weights double-count). Same AO-metric quantity as the
    ``contamination`` of :func:`purified_s_core`. l-purity is NECESSARY
    but not SUFFICIENT for rotational closure -- only the commutators
    measure closure itself (facts 8, 19).
    Returns ``(names, w)`` with ``names = ['s', 'p', 'd', ...]``.
    """
    import re

    labels = [line.split()[-1] for line in mol.ao_labels()]
    groups: dict[str, list[int]] = {}
    for i, lbl in enumerate(labels):
        m = re.match(r"\d+([spdfgh])", lbl)
        groups.setdefault(m.group(1) if m else lbl, []).append(i)

    names = list(groups)
    w = np.zeros((c.shape[1], len(names)))
    for k, name in enumerate(names):
        idx = groups[name]
        ss = s[np.ix_(idx, idx)]
        w[:, k] = np.einsum("ij,ij->j", c[idx, :], ss @ c[idx, :]).real
    return names, w


def complete_orthonormal(v: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Append an orthonormal-in-metric-``s`` complement to the columns of ``v``."""
    ev, u = np.linalg.eigh(s)
    s_half = (u * np.sqrt(ev)) @ u.T
    s_ihalf = (u / np.sqrt(ev)) @ u.T
    vx = s_half @ v
    resid = np.eye(v.shape[0]) - vx @ vx.conj().T
    p, sv, _ = np.linalg.svd(resid)
    c = s_ihalf @ p[:, sv > 1e-8]
    return np.hstack([v, c])


def build_full_mo(mf, core_idx: list[int], act_idx: list[int], c_core: np.ndarray = None) -> np.ndarray:
    """Full orthonormal MO matrix with chosen core and active columns placed
    so that CASCI/CASSCF's window ``[ncore : ncore+ncas]`` selects the active
    space.

    Core columns are canonical SCF MOs by default; ``c_core`` supplies
    alternative columns (e.g. fractional-occupation core). Active columns
    are the RAW AO functions of the target shell -- exactly rotation-closed,
    essential for the L-quantum-number analysis.

    Construction is block-safe in the AO metric S (fact 14): (1) the core
    block is Loewdin-orthonormalized within its own span (span- and
    density-preserving); (2) the active columns are projected against the
    core span and renormalized (raw valence AOs overlap the canonical
    core's same-l members: <2p_MO|3p_AO> can reach ~0.3) -- a same-l,
    m-diagonal correction, so window l-purity and m-set structure
    survive; (3) the complement is appended through an exact projector --
    the naive `I - P` filter is exact only for S-orthonormal v and
    silently returns extra columns otherwise (observed: (nao, nao+6)
    on Si).
    """
    v = _core_window_block(mf, core_idx, act_idx, c_core)
    if v.shape[1] == mf.mo_coeff.shape[0]:
        return v
    return complete_orthonormal(v, mf.get_ovlp())


def _core_window_block(mf, core_idx: list[int], act_idx: list[int],
                       c_core: np.ndarray = None) -> np.ndarray:
    """Steps (1)+(2) of :func:`build_full_mo`: the S-orthonormal
    core+window block that both it and :func:`radial_basis` start from."""
    n = mf.mo_coeff.shape[0]
    s = mf.get_ovlp()
    core = mf.mo_coeff[:, core_idx] if c_core is None else c_core
    act = np.eye(n)[:, act_idx]

    core = lowdin(core, s)

    p_core = core @ core.conj().T @ s
    act = act - p_core @ act
    act = lowdin(act, s)

    return np.hstack([core, act])


def lm_sectors(mol) -> dict[tuple[int, int], list[int]]:
    """AO rows grouped by (l, member): the EXACTLY S-orthogonal sectors of
    a one-center spherical basis -- the AO space is a direct sum of radial
    chains, one per (l, member); same-sector functions of different n
    OVERLAP, different-sector functions never do (angular orthogonality,
    fact 16e). This decomposition is what makes a per-sector treatment
    (:func:`radial_basis`) globally S-orthonormal by construction.
    """
    sectors: dict[tuple[int, int], list[int]] = {}
    loc = mol.ao_loc_nr()
    for ib in range(mol.nbas):
        l = abs(mol.bas_angular(ib))
        nctr = mol.bas_nctr(ib)
        for c in range(nctr):
            for member in range(2 * l + 1):
                # fused shells (nctr > 1) are contraction-outer, m-inner
                sectors.setdefault((l, member), []).append(
                    loc[ib] + c * (2 * l + 1) + member)
    return sectors


def radial_basis(mf, core_idx: list[int], act_idx: list[int],
                 c_core: np.ndarray = None):
    """Full MO matrix plus per-slot (l, member) labels -- the substrate for
    :func:`radial_mask`. ``c_core`` overrides the core block; pass
    :func:`purified_core_orbitals` output for the exactly-closed
    substrate (the fact-17 fix).

    Everything non-core is built in a MEMBER-PRESERVING gauge (sequential
    S-metric Gram-Schmidt, not a symmetric Loewdin): the Loewdin gauge
    within a degenerate member set is an arbitrary rotation --
    span-preserving, but it reassigns which column is which (l, m), and
    here column identity IS the product (fact 16e, README §3.3). Window
    slots carry their AO's sector label; STRONG sector residuals
    (>= 8e-2 of S-weight; core tails <= 3e-6) join with their label; the
    near-linearly-dependent needles (Sc 6p: 1.7e-4) pad through
    :func:`complete_orthonormal` with label ``None`` -- never rotation
    partners.

    Returns ``(mo_full, labels)`` with ``labels[i] = (l, member) | None``.
    """
    mol = mf.mol
    n = mf.mo_coeff.shape[0]
    s = mf.get_ovlp()
    sec = lm_sectors(mol)
    row_sector = {row: key for key, rows in sec.items() for row in rows}

    core = mf.mo_coeff[:, core_idx] if c_core is None else c_core
    core = lowdin(core, s)

    def sgs(cols, base):
        """S-metric Gram-Schmidt (two refinement passes) of the columns in
        ``cols`` on top of the S-orthonormal block(s) in ``base``;
        preserves column identity (unlike a symmetric Loewdin)."""
        out = []
        cur = np.hstack(base) if base else None
        for c in cols:
            c = np.asarray(c).reshape(-1, 1)  # uniform (n, 1) columns
            for _ in range(2):
                if cur is not None:
                    c = c - cur @ (cur.conj().T @ s @ c)
                c = c / float(np.real(np.sqrt((c.conj().T @ s @ c)[0, 0])))
            out.append(c)
            cur = c if cur is None else np.hstack([cur, c])
        return out

    # window in member order: raw AO minus the core span, sequentially
    # orthogonalized -- column j IS the (l, m) of act_idx[j] by construction
    p_core = core @ core.conj().T @ s
    window = sgs([np.eye(n)[:, ao] - p_core @ np.eye(n)[:, ao]
                  for ao in act_idx], [core])
    v = np.hstack([core] + window)
    labels: list = [None] * len(core_idx) + [row_sector[ao] for ao in act_idx]

    # strong sector residuals, same member-preserving treatment. Strength
    # cutoff 8e-2 / needle class 1.7e-4: docstring + fact 16e (the 3s of
    # cc-pVDZ C survives at 0.085; the Sc 6p needle does not)
    p_v = v @ v.conj().T @ s
    strong, strong_labels = [], []
    for key, rows in sorted(sec.items()):
        block = np.eye(n)[:, rows]
        r = block - p_v @ block
        gr = r.conj().T @ s @ r
        ev, u = np.linalg.eigh(gr)
        keep = ev > 1e-4
        if keep.any():
            strong += [r @ (u[:, j] / np.sqrt(ev[j])) for j in np.where(keep)[0]]
            strong_labels += [key] * int(keep.sum())
    if strong:
        v = np.hstack([v] + sgs(strong, [v]))
        labels += strong_labels

    mo = v if v.shape[1] == n else complete_orthonormal(v, s)
    labels += [None] * (mo.shape[1] - len(labels))

    # gauge check: every labeled column must be dominated by its own
    # sector. A column that fails is a NEEDLE -- its sector content was
    # consumed by the orthogonalization against the substrate -- and is
    # DEMOTED to unlabeled padding (fact 16e: strong directions labeled,
    # needles padded; never a rotation partner). The needle helper of
    # ticket 20/26: deterministic, no flake.
    demoted = []
    for i, lab in enumerate(labels):
        if lab is None:
            continue
        rows = sec[lab]
        c = mo[rows, i]
        w = float(np.real(c.conj().T @ s[np.ix_(rows, rows)] @ c))
        if w <= 0.9:
            labels[i] = None
            demoted.append((i, lab, round(w, 3)))
    if demoted:
        print(f"radial_basis: needle directions demoted to unlabeled "
              f"(fact 16e): {demoted}")
    assert mo.shape[1] == n
    return mo, labels


def radial_mask(labels: list, ncore: int, ncas: int):
    """``uniq_var_indices`` mask allowing ONLY same-(l, member) rotations
    into a window slot: pure RADIAL mixing (3d_m may borrow 4d_m/5d_m; the
    4s any s), l-content and m-membership of every slot untouchable, and
    no core<->window exchange at all (core slots carry no label).

    Physics caveat (measured, fact 16): exact rotational closure of the
    window additionally requires the radial rotation applied IDENTICALLY
    to all m partners; independent per-member radial mixing leaves
    [H, Lz] nonzero at the member-difference level. SA over the term's
    M_L set is the symmetry-restoring companion rung.
    """
    lab = list(labels)

    def mask(self, nmo, ncore_, ncas_, frozen):
        m = np.zeros((nmo, nmo), dtype=bool)
        for w in range(ncore, ncore + ncas):
            same = np.array([lab[u] is not None and lab[u] == lab[w]
                             for u in range(nmo)])
            same[w] = False
            m[same, w] = True
        return m

    return mask


class Window(NamedTuple):
    """A candidate active space in one of two representations.

    ``ao_rows``: AO indices -- m-complete by construction (shell /
    primitive selectors); materialize through :meth:`columns`.
    ``c_win``: orthonormal orbital columns (SCF-MO / NO selectors);
    used directly as the active columns. Exactly one is not None.
    ``m_complete``: the construction's closure GUARANTEE (shell /
    primitive: exact; scf_mo: never; no: only if the source density was
    m-scalar -- verify with window_closure, this flag records intent).
    """

    label: str
    ao_rows: list[int] | None
    c_win: np.ndarray | None
    m_complete: bool

    def columns(self, mf, ncore: int, closed: list[int] | None = None,
                c_core: np.ndarray = None) -> np.ndarray:
        """Full MO substrate with this window's active columns placed at
        ``[ncore : ncore + ncas]`` -- the object CASCI/CASSCF consume.

        closed : the OCCUPIED core indices (canonical order is
            meaningless on near-degenerate bases — pass
            ``occupied_core_idx(mf, ncore)`` for AO windows).
        c_core : alternative core columns (e.g. the pur substrate);
            defaults to the first ``ncore`` columns of ``mf``.
        """
        if self.ao_rows is not None:
            closed = closed if closed is not None else list(range(ncore))
            return build_full_mo(mf, closed, self.ao_rows, c_core=c_core)
        core = c_core if c_core is not None else mf.mo_coeff[:, :ncore]
        return complete_orthonormal(np.hstack([core, self.c_win]),
                                    mf.get_ovlp())


def window_shell(mol, mf, ncore: int, *, patterns) -> Window:
    """Contracted-shell window (AO rows; drivers.select_valence_ao)."""
    return Window(f"shell {patterns}", select_valence_ao(mol, patterns),
                  None, True)


def window_primitive(mol, mf, ncore: int, *, bands) -> Window:
    """Primitive-band window (AO rows), ``bands = {l: count | (lo, hi)}``.

    Bands must be contiguous sets of the measured open-MO primitive
    weights (the scan protocol); select_primitives enforces m-completeness
    per primitive."""
    rows: list[int] = []
    for l in sorted(bands):
        sel = bands[l]
        if isinstance(sel, tuple):
            rows += select_primitives(mol, l, exponent_band=sel)
        else:
            rows += select_primitives(mol, l, count=sel)
    return Window(f"primitive {bands}", rows, None, True)


def window_scf_mo(mol, mf, ncore: int, *, open_cols, rank, k: int) -> Window:
    """Raw-SCF-MO window (ticket 24 methodology B): the open reference
    MOs plus the ``k`` top-``rank`` virtuals (the caller owns the
    ranking metric -- atom-specific valence-band rows live in the
    scripts). NEVER rotation-closed: the span certificate measures 7-9
    at every practical size (README 7.5)."""
    cols = list(open_cols) + list(rank[:k])
    return Window(f"scf_mo +{k}", None, mf.mo_coeff[:, cols], False)


def window_no(mol, mf, ncore: int, *, c_act, dets, vec, n_max: int,
              weights=None, deg_tol: float = 1e-7) -> Window:
    """NO window: term-average roots -> 1-RDM -> top-occupation NOs.

    ``vec``/``weights`` select and average the roots (pass the whole
    term's multiplet for an m-scalar density; a single root yields an
    m-split window -- legal but NOT rotation-closed, see the module
    docstring). Verify every output with operators.window_closure.
    """
    norb = c_act.shape[1]
    gamma = ci_rdm1(dets, vec, norb, weights)
    c_no, occ = natural_orbitals(c_act, gamma)
    c_win, kept = select_no_window(c_no, occ, n_max, deg_tol)
    n_cols = sum(len(g) for g in kept)
    return Window(f"no {n_cols} cols", None, c_win, False)


SELECTORS = {
    "shell": window_shell,
    "primitive": window_primitive,
    "scf_mo": window_scf_mo,
    "no": window_no,
}
