"""Natural-orbital mathematics: the CI 1-RDM, the NOs, occupation
groups, and the NO window selection (ticket 22). The Window object
and the window-kind registry (shell / primitive / scf_mo / no) live
in atomterm.windows.

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

import numpy as np

from atomterm.secondquant import annihilate, create

__all__ = [
    "ci_rdm1",
    "natural_orbitals",
    "occupation_groups",
    "select_no_window",
]


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

