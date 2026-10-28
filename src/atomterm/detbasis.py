"""First-quantized many-electron operator matrices over a
Slater-determinant basis: operators enter as one-body matrices and
two-electron integrals in the spin-orbital basis, and the matrix
elements follow from the explicit Slater-Condon rules.

Conventions
-----------
* Spatial orbitals ``0..norb-1``. Spin-orbital ``p`` = spatial ``p``
  with alpha spin; ``norb + p`` = the same spatial orbital with beta
  spin.
* A determinant is a sorted tuple of occupied spin-orbitals -- the
  sparse form of an occupation bit string over the ``2*norb`` slots.
* One-body operator: matrix ``h[p, q] = <p|h|q>``, shape
  ``(2*norb, 2*norb)``.
* Two-body integrals, chemist notation:
  ``G[p, q, r, s] = (pq|rs) = Integral p(1) q(2) r(1) s(2)``, shape
  ``(2*norb,)*4``; the operator is ``Sum_{i<j} g(i, j)``.

The determinant ordering here (lexicographic over the (alpha, beta)
combinations) is NOT PySCF's FCI string order -- they differ from
``norb >= 4``; the mapping lives in ``cells.pyscf_ci_permutation``.
This module is deliberately NumPy-only (no pyscf import), testable in
complete isolation.
"""

from __future__ import annotations

from itertools import combinations
import bisect

import numpy as np

__all__ = [
    "determinant_basis",
    "one_body_matrix",
    "two_body_matrix",
    "slater_condon_matrix",
    "commutator",
    "frobenius_norm",
]


def determinant_basis(norb: int, na: int, nb: int) -> list[tuple[int, ...]]:
    """All determinants with ``na`` alpha and ``nb`` beta electrons:
    sorted tuples of occupied spin-orbitals.

    Returns
    -------
    Lexicographic over the (alpha combination, beta combination);
    alpha orbitals are ``0..norb-1``, beta are ``norb..2*norb-1``.

    Examples
    --------
    One alpha and one beta electron on two spatial orbitals:

    >>> determinant_basis(2, 1, 1)
    [(0, 2), (0, 3), (1, 2), (1, 3)]
    # (a0 b0)  (a0 b1)  (a1 b0)  (a1 b1)
    """
    dets = []
    for occ_a in combinations(range(norb), na):
        for occ_b in combinations(range(norb), nb):
            dets.append(tuple(sorted(occ_a)) + tuple(o + norb for o in occ_b))
    return dets


def _pop_sign(occ: list[int], orb: int) -> int:
    """Sign of annihilating ``orb`` from the ordered occupation
    ``occ``: (-1)^(its position)."""
    i = occ.index(orb)
    return -1 if i % 2 else 1


def _insert_sign(occ: list[int], orb: int) -> int:
    """Sign of creating ``orb`` into the ordered occupation ``occ``:
    (-1)^(# orbitals before it)."""
    n_before = sum(1 for o in occ if o < orb)
    return -1 if n_before % 2 else 1


def _transition_sign(occ_j: list[int], ann: list[int], cre: list[int]) -> int:
    """Sign of ``<I| a+_c1 a+_c2 ... a_i2 a_i1 |J>``.

    Annihilations run in ascending orbital order (rightmost operator
    acts first), creations in descending order, so that the paired
    Slater-Condon values follow the standard ``<ab||ij>`` ordering.
    """
    tmp = list(occ_j)
    sign = 1
    for o in sorted(ann):
        sign *= _pop_sign(tmp, o)
        tmp.remove(o)
    for o in sorted(cre, reverse=True):
        sign *= _insert_sign(tmp, o)
        bisect.insort(tmp, o)
    return sign


def one_body_matrix(
    dets: list[tuple[int, ...]], one_body: np.ndarray
) -> np.ndarray:
    """The ``O1 = Sum_pq h[p,q] a+_p a_q`` matrix over ``dets``
    (conventions and dtype: :func:`slater_condon_matrix`)."""
    return slater_condon_matrix(dets, one_body=one_body)


def two_body_matrix(
    dets: list[tuple[int, ...]], two_body: np.ndarray
) -> np.ndarray:
    """The ``O2 = 1/2 Sum_pqrs G[p,q,r,s] a+_p a+_q a_s a_r`` matrix
    over ``dets``, ``G`` in the module's chemist notation
    (conventions and dtype: :func:`slater_condon_matrix`)."""
    return slater_condon_matrix(dets, two_body=two_body)


def slater_condon_matrix(
    dets: list[tuple[int, ...]],
    one_body: np.ndarray | None = None,
    two_body: np.ndarray | None = None,
) -> np.ndarray:
    """Full many-electron matrix of ``Sum_i h(i) + Sum_{i<j} g(i, j)``
    in ONE fused pass -- the production path (handles the 10^5-detector
    sectors; :func:`_slater_condon_slow`, the O(ndet^2) Python loop it
    replaced, stays as the elementwise reference pinned by tests).

    Parameters
    ----------
    dets : determinant basis (see :func:`determinant_basis`).
    one_body : ``(2*norb, 2*norb)`` matrix ``h[p, q]`` (optional).
    two_body : ``(2*norb,)*4`` chemist-notation integrals (optional).

    Returns
    -------
    ``(len(dets), len(dets))`` matrix, ``M[I, J] = <I|O|J>``; dtype
    follows the inputs (real operators stay real -- a complex container
    would double memory traffic and slow the eigendecomposition). For
    the parts separately: :func:`one_body_matrix` /
    :func:`two_body_matrix`.
    """
    if one_body is None and two_body is None:
        raise ValueError("at least one of one_body / two_body is required")
    n = len(dets)
    nelec = len(dets[0])
    # the operator fixes the spin-orbital dimension: sectors without beta
    # electrons never index the beta block (max det index < norb)
    n2 = (one_body if one_body is not None else two_body).shape[0]
    occ = np.zeros((n, n2), dtype=bool)
    for i, d in enumerate(dets):
        occ[i, list(d)] = True
    # prefix counts: cum[j, t] = # occupied spin-orbitals < t in det j
    cum = np.zeros((n, n2), dtype=np.int64)
    np.cumsum(occ[:, :-1], axis=1, out=cum[:, 1:])
    dtype = complex if (np.iscomplexobj(one_body if one_body is not None
                                        else two_body)) else float
    mat = np.zeros((n, n), dtype=dtype)

    # ---- diagonal (k = 0) ------------------------------------------------
    diag = np.zeros(n, dtype=dtype)
    if one_body is not None:
        diag += occ @ np.diag(one_body).astype(dtype)
    if two_body is not None:
        if nelec >= 2:
            occ_idx = np.array([list(d) for d in dets])      # (n, nelec)
            comb = np.array(
                list(combinations(range(nelec), 2)))          # (nc, 2)
            p = occ_idx[:, comb[:, 0]]                        # (n, nc)
            q = occ_idx[:, comb[:, 1]]
            diag += 0.5 * (
                two_body[p, q, p, q] + two_body[q, p, q, p]
                - two_body[q, p, p, q] - two_body[p, q, q, p]
            ).sum(axis=1)
    mat[np.diag_indices(n)] = diag

    # ---- off-diagonal (k = 1, 2), chunked over rows ----------------------
    chunk = max(1, min(n, int(2e7) // (n * n2)))
    orb = np.arange(n2)
    for lo in range(0, n, chunk):
        hi = min(n, lo + chunk)
        occ_i = occ[lo:hi]                                    # (m, n2)
        k = nelec - (occ_i[:, None, :] & occ[None, :, :]).sum(-1)
        rows, cols = np.nonzero((k == 1) | (k == 2))
        if rows.size == 0:
            continue
        ii = rows + lo
        jj = cols
        occ_j_c = occ[jj]                                     # (P, n2)
        occ_i_c = occ[ii]
        ann_mask = occ_j_c & ~occ_i_c                         # k entries
        cre_mask = occ_i_c & ~occ_j_c
        k_c = k[rows, cols]

        # per-pair excitation indices, ascending within each set
        ann_r, ann_c = np.nonzero(ann_mask)
        cre_r, cre_c = np.nonzero(cre_mask)
        ann = np.zeros((ii.size, 2), dtype=np.int64)
        cre = np.zeros((ii.size, 2), dtype=np.int64)
        ann[ann_r, _pair_slot(ann_r, ann_c)] = ann_c
        cre[cre_r, _pair_slot(cre_r, cre_c)] = cre_c
        # _pair_slot: 0 for the first (smaller) member of the excitation
        # set, 1 for the second; k=1 pairs use only slot 0 (slot 1 stays
        # 0 = dummy)

        sign = np.ones(ii.size)
        pref_j = cum[jj]                                      # (P, n2)
        m1 = k_c == 1
        if m1.any():
            sign[m1] = _sign_single(pref_j[m1], ann[m1, 0], cre[m1, 0])
        m2 = ~m1
        if m2.any():
            sign[m2] = _sign_double(pref_j[m2], ann[m2, 0], ann[m2, 1],
                                    cre[m2, 0], cre[m2, 1])

        val = np.zeros(ii.size, dtype=dtype)
        if one_body is not None:
            val[m1] += one_body[cre[m1, 0], ann[m1, 0]]
        if two_body is not None:
            if m1.any():
                i1, a1 = ann[m1, 0], cre[m1, 0]
                common = occ_j_c[m1] & occ_i_c[m1]            # (P1, n2)
                g = two_body                                 # shorthand
                tb = 0.5 * (
                    g[a1[:, None], orb[None, :], i1[:, None], orb[None, :]]
                    + g[orb[None, :], a1[:, None], orb[None, :], i1[:, None]]
                    - g[orb[None, :], a1[:, None], i1[:, None], orb[None, :]]
                    - g[a1[:, None], orb[None, :], orb[None, :], i1[:, None]]
                )
                val[m1] += (tb * common).sum(axis=1)
            if m2.any():
                i2, j2 = ann[m2, 0], ann[m2, 1]
                a2, b2 = cre[m2, 0], cre[m2, 1]
                val[m2] += 0.5 * (
                    two_body[a2, b2, i2, j2] + two_body[b2, a2, j2, i2]
                    - two_body[b2, a2, i2, j2] - two_body[a2, b2, j2, i2]
                )
        nz = val != 0
        mat[ii[nz], jj[nz]] = sign[nz] * val[nz]
    return mat


def _sign_single(pj: np.ndarray, i: np.ndarray, a: np.ndarray) -> np.ndarray:
    """Closed-form k=1 transition sign (annihilate ``i``, create
    ``a``) from determinant J's prefix counts ``pj`` (``cum`` rows):

    ``(-1)^{cum[j,i] + cum[j,a] - (i < a)}``.
    """
    r = np.arange(len(i))
    return (-1.0) ** (pj[r, i] + pj[r, a] - (i < a).astype(np.int64))


def _sign_double(pj: np.ndarray, i: np.ndarray, j: np.ndarray,
                 a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Closed-form k=2 transition sign (annihilate ``i < j``, create
    ``a < b``, ``b`` inserted first) from determinant J's prefix
    counts:

    ``(-1)^{cum[j,i] + cum[j,j_o] - 1 + cnt(b) + cnt(a)}`` with
    ``cnt(y) = cum[j,y] - (i < y) - (j < y)``.
    """
    r = np.arange(len(i))
    cnt = (lambda y: pj[r, y]                        # noqa: E731
           - (i < y).astype(np.int64) - (j < y).astype(np.int64))
    return (-1.0) ** (pj[r, i] + pj[r, j] - 1 + cnt(b) + cnt(a))


def _pair_slot(rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    """Map each (pair, orbital) entry of ``np.nonzero`` output to its
    slot in the ``(..., 2)`` excitation-index arrays: 0 = smaller
    member of the pair, 1 = larger.

    The excitation sets have at most TWO members because the pairs
    were pre-filtered to k <= 2 (a two-body operator has zero matrix
    elements beyond double excitations); ``np.nonzero`` returns
    ascending column order within each row, so "first occurrence in
    the row" = the smaller orbital.

    Examples
    --------
    >>> _pair_slot(np.array([0, 0, 1, 2, 2]), np.array([3, 7, 5, 2, 9]))
    array([0, 1, 0, 0, 1])
    """
    out = np.zeros(len(rows), dtype=np.int64)
    order = np.lexsort((cols, rows))
    rs = rows[order]
    rk = np.zeros(len(rs), dtype=np.int64)
    rk[1:] = (rs[1:] == rs[:-1]).astype(np.int64)  # 2nd within row -> 1
    out[order] = rk
    return out


def _slater_condon_slow(
    dets: list[tuple[int, ...]],
    one_body: np.ndarray | None = None,
    two_body: np.ndarray | None = None,
) -> np.ndarray:
    """The original O(ndet^2) Python loop -- kept as the elementwise
    reference for the vectorized :func:`slater_condon_matrix`
    (pinned by tests; sign engine: :func:`_transition_sign`)."""
    ndet = len(dets)
    mat = np.zeros((ndet, ndet), dtype=complex)

    for j, occ_j in enumerate(dets):
        occ_j_list = list(occ_j)
        for i, occ_i in enumerate(dets):
            common = sorted(set(occ_i) & set(occ_j))
            ann = sorted(set(occ_j) - set(occ_i))
            cre = sorted(set(occ_i) - set(occ_j))
            k = len(ann)
            if k > 2:
                continue

            if k == 0:
                val = 0.0 + 0.0j
                if one_body is not None:
                    val += sum(one_body[p, p] for p in occ_i)
                if two_body is not None:
                    # fully antisymmetrized diagonal element (valid for any
                    # two-body kernel, not only particle-symmetric ones)
                    for p, q in combinations(occ_i, 2):
                        val += 0.5 * (
                            two_body[p, q, p, q]
                            + two_body[q, p, q, p]
                            - two_body[q, p, p, q]
                            - two_body[p, q, q, p]
                        )
            elif k == 1:
                i_o, a_o = ann[0], cre[0]
                val = 0.0 + 0.0j
                if one_body is not None:
                    val += one_body[a_o, i_o]
                if two_body is not None:
                    for o in common:
                        val += 0.5 * (
                            two_body[a_o, o, i_o, o]
                            + two_body[o, a_o, o, i_o]
                            - two_body[o, a_o, i_o, o]
                            - two_body[a_o, o, o, i_o]
                        )
            else:  # k == 2
                if two_body is None:
                    continue
                i_o, j_o = ann
                a_o, b_o = cre
                val = 0.5 * (
                    two_body[a_o, b_o, i_o, j_o]
                    + two_body[b_o, a_o, j_o, i_o]
                    - two_body[b_o, a_o, i_o, j_o]
                    - two_body[a_o, b_o, j_o, i_o]
                )

            if val != 0.0:
                mat[i, j] = _transition_sign(occ_j_list, ann, cre) * val

    return mat


def commutator(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """``[A, B] = AB - BA`` for plain matrices."""
    return a @ b - b @ a


def frobenius_norm(a: np.ndarray) -> float:
    """Frobenius norm as a float (magnitude of a matrix)."""
    return float(np.linalg.norm(a))
