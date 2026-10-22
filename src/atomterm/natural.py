"""Natural-orbital windows and the window-selector registry (ticket 22).

The closure theorem (fact 20) licenses ANY per-l radial subspace
replicated over ALL members as window units. NOs are legal units
exactly when they arrive as complete m-sets -- i.e. when the 1-RDM they
diagonalize is m-scalar per l, which holds iff the CI roots averaged
into it form a rotation-invariant set (a term average over the whole
L multiplet, all in one M_S sector). A single |L, M_L> root carries an
m-DIAGONAL but m-UNEQUAL density (p^2 3P M_L=0: occ(m=0)=2,
occ(+-1)=0 -- pinned in tests/test_natural.py), so occupation selection
on a single root is NOT symmetry-safe. The one-particle certificate
(:func:`atomterm.operators.window_closure`) is the arbiter: run it on
every NO window, never assume.

Spin-orbital conventions (repo-wide): alpha = 0..norb-1,
beta = norb..2norb-1; determinants are sorted spin-orbital tuples.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from atomterm.secondquant import annihilate, create

__all__ = [
    "Window",
    "ci_rdm1",
    "natural_orbitals",
    "occupation_groups",
    "select_no_window",
    "window_shell",
    "window_primitive",
    "window_scf_mo",
    "window_no",
    "SELECTORS",
]


class Window(NamedTuple):
    """A candidate active space in one of two representations.

    ``ao_rows``: AO indices -- m-complete by construction (shell /
    primitive selectors); feed through drivers.build_full_mo.
    ``c_win``: orthonormal orbital columns (SCF-MO / NO selectors);
    use directly as the active columns. Exactly one is not None.
    ``m_complete``: the construction's closure GUARANTEE (shell /
    primitive: exact; scf_mo: never; no: only if the source density was
    m-scalar -- verify with window_closure, this flag records intent).
    """

    label: str
    ao_rows: list[int] | None
    c_win: np.ndarray | None
    m_complete: bool


def ci_rdm1(dets, vec, norb: int, weights=None) -> np.ndarray:
    """Spin-summed spatial 1-RDM of CI root(s), active-orbital basis.

    vec: (ndet,) or (ndet, nroots) in the ``dets`` ordering;
    weights: per-root averaging weights (default: the first root only).
    gamma[p, q] = <a+_p a_q> summed over alpha and beta. Signs come from
    the tested secondquant rules; for a single determinant gamma is the
    projector onto its occupied orbitals (tr = na + nb).
    """
    vec = np.asarray(vec)
    if vec.ndim == 1:
        vec = vec[:, None]
    if weights is None:
        weights = np.eye(1, vec.shape[1])[0]
    weights = np.asarray(weights, dtype=float)

    def split(det):
        return (tuple(p for p in det if p < norb),
                tuple(p - norb for p in det if p >= norb))

    index = {split(det): i for i, det in enumerate(dets)}

    # transition table (root-independent, built once): j -> annih(q),
    # create(p) -> i with the combined fermionic sign
    js, ps, qs, sgs, tgt = [], [], [], [], []
    for j, det in enumerate(dets):
        occ_a, occ_b = split(det)
        for occ, new_a in ((occ_a, False), (occ_b, True)):
            for q in occ:
                sg_q, o1 = annihilate(occ, q)
                for p in range(norb):
                    made = create(o1, p)
                    if made is None:  # p already occupied
                        continue
                    sg_p, o2 = made
                    key = ((occ_a, o2) if new_a else (o2, occ_b))
                    i = index.get(key)
                    if i is None:
                        continue
                    js.append(j)
                    ps.append(p)
                    qs.append(q)
                    sgs.append(sg_q * sg_p)
                    tgt.append(i)
    js = np.asarray(js)
    ps = np.asarray(ps)
    qs = np.asarray(qs)
    sgs = np.asarray(sgs, dtype=float)
    tgt = np.asarray(tgt)

    gamma = np.zeros((norb, norb))
    for r in range(vec.shape[1]):
        w = float(weights[r])
        if w == 0.0:
            continue
        c = vec[:, r]
        contrib = sgs * c[tgt] * c[js]
        flat = np.bincount(ps * norb + qs, weights=contrib,
                           minlength=norb * norb)
        gamma += w * flat.reshape(norb, norb)
    return gamma


def natural_orbitals(c_act: np.ndarray,
                     gamma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """NOs of a spin-summed spatial density in an orthonormal basis.

    Returns (c_no, occ): NO columns in the AO basis and occupations,
    sorted by DESCENDING occupation. c_act must be orthonormal (the
    active space is, fact 4): the density eigenvectors ARE the NOs, no
    metric dance needed inside the span.
    """
    w, v = np.linalg.eigh(gamma)
    order = np.argsort(w)[::-1]
    return c_act @ v[:, order], w[order]


def occupation_groups(occ, deg_tol: float = 1e-7) -> list[list[int]]:
    """Cluster descending occupations into degenerate groups.

    For an m-scalar density each group IS an m-set (2l+1 equal
    occupations per radial); window selection must never split one
    (the fact-19 AVAS lesson: per-member thresholds build
    rotation-defective spans)."""
    occ = np.asarray(occ)
    groups: list[list[int]] = [[0]]
    for k in range(1, len(occ)):
        if occ[k] >= occ[groups[-1][0]] - deg_tol:
            groups[-1].append(k)
        else:
            groups.append([k])
    return groups


def select_no_window(c_no: np.ndarray, occ, n_max: int,
                     deg_tol: float = 1e-7) -> tuple[np.ndarray, list[list[int]]]:
    """Top-occupation NO window, whole occupation groups only.

    Takes groups by descending occupation while they fit ``n_max`` (the
    first group is always taken -- an m-set is the smallest legal
    window unit). Returns (c_win, groups): orthonormal NO columns and
    the kept group index lists (for reporting the m-set structure).
    """
    groups = occupation_groups(occ, deg_tol)
    kept: list[list[int]] = []
    size = 0
    for g in groups:
        if kept and size + len(g) > n_max:
            break
        kept.append(g)
        size += len(g)
    cols = [k for g in kept for k in g]
    return c_no[:, cols], kept


# -------------------------------------------------------------- registry
# One signature per methodology: (mol, mf, ncore, **spec) -> Window.
# Candidates are combinable with any fold / treatment (the axes product
# ticket 20 consumes). AO-row windows carry m_complete=True by
# construction; orbital-column windows carry their measured/declared
# status in the flag (verify with the certificate either way).

def window_shell(mol, mf, ncore: int, *, patterns) -> Window:
    """Contracted-shell window (AO rows; drivers.select_valence_ao)."""
    from atomterm.drivers import select_valence_ao
    return Window(f"shell {patterns}", select_valence_ao(mol, patterns),
                  None, True)


def window_primitive(mol, mf, ncore: int, *, bands) -> Window:
    """Primitive-band window (AO rows), ``bands = {l: count | (lo, hi)}``.

    Bands must be contiguous sets of the measured open-MO primitive
    weights (the scan protocol); select_primitives enforces m-completeness
    per primitive."""
    from atomterm.drivers import select_primitives
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
