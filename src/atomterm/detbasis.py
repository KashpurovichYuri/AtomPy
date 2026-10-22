"""First-quantized many-electron operator matrices over a Slater-determinant basis.

No second quantization: operators are specified by one-body matrices and
two-electron integrals in the spin-orbital basis, and many-electron matrix
elements follow from the explicit Slater-Condon rules applied to listed
determinants.

Conventions
-----------
* Spatial orbitals ``0..norb-1``. Spin-orbital ``p`` = spatial ``p`` with alpha
  spin; ``norb + p`` = the same spatial orbital with beta spin.
* A determinant is a sorted tuple of occupied spin-orbitals.
* One-body operator: matrix ``h[p, q] = <p|h|q>``, shape ``(2*norb, 2*norb)``.
* Two-body "integrals": chemist notation ``G[p, q, r, s] = (pq|rs) =
 Integral p(1) q(2) r(1) s(2)``, shape ``(2*norb,)*4``. The corresponding
  two-electron operator is ``Sum_{i<j} g(i, j)``; ``G`` may be complex and
  non-symmetric (e.g. the composite L.L operator), the rules only need the
  matrix elements.

This module is deliberately NumPy-only (no pyscf import) so that the
Slater-Condon machinery is testable in complete isolation.
"""

from __future__ import annotations

from itertools import combinations
import bisect

import numpy as np

__all__ = [
    "determinant_basis",
    "slater_condon_matrix",
    "commutator",
    "frobenius_norm",
]


def determinant_basis(norb: int, na: int, nb: int) -> list[tuple[int, ...]]:
    """All determinants with ``na`` alpha and ``nb`` beta electrons.

    Returns sorted tuples of occupied spin-orbitals (alpha block first in
    index value, i.e. alpha orbitals are 0..norb-1, beta are norb..2norb-1).
    """
    dets = []
    for occ_a in combinations(range(norb), na):
        for occ_b in combinations(range(norb), nb):
            dets.append(tuple(sorted(occ_a)) + tuple(o + norb for o in occ_b))
    return dets


def _pop_sign(occ: list[int], orb: int) -> int:
    """Sign of annihilating ``orb`` from the ordered occupation ``occ``."""
    i = occ.index(orb)
    return -1 if i % 2 else 1


def _insert_sign(occ: list[int], orb: int) -> int:
    """Sign of creating ``orb`` into the ordered occupation ``occ``."""
    n_before = sum(1 for o in occ if o < orb)
    return -1 if n_before % 2 else 1


def _transition_sign(occ_j: list[int], ann: list[int], cre: list[int]) -> int:
    """Sign of <I| a+_c1 a+_c2 ... a_i2 a_i1 |J>.

    Annihilations run in ascending orbital order (rightmost operator acts
    first), creations in descending order, so that the paired Slater-Condon
    values below follow the standard ``<ab||ij>`` ordering.
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


def slater_condon_matrix(
    dets: list[tuple[int, ...]],
    one_body: np.ndarray | None = None,
    two_body: np.ndarray | None = None,
) -> np.ndarray:
    """Full many-electron matrix of ``Sum_i h(i) + Sum_{i<j} g(i, j)``.

    Parameters
    ----------
    dets : determinant basis (see :func:`determinant_basis`).
    one_body : ``(2*norb, 2*norb)`` matrix ``h[p, q]`` (optional).
    two_body : ``(2*norb,)*4`` chemist-notation integrals (optional).

    Returns
    -------
    complex ndarray of shape ``(len(dets), len(dets))`` with
    ``M[I, J] = <I| O |J>``.

    Vectorized: pair classification by chunked occupation broadcasting,
    closed-form Slater signs from per-determinant prefix counts, and
    gathered Slater-Condon values. Elementwise-identical to the
    reference loop (:func:`_slater_condon_reference`, pinned by tests);
    the loop version is O(ndet^2) Python per pair and dies in the
    ~10^4-detector range, this handles 10^5 (the U campaign's sector
    sizes).
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
    # real operators stay real: every operator of this project is
    # real-valued, and the complex dtype costs 2x memory traffic plus a
    # ~2.5x eigendecomposition at the 10^4-detector scale
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

        # per-pair excitation indices, ascending within each pair
        ann_r, ann_c = np.nonzero(ann_mask)
        cre_r, cre_c = np.nonzero(cre_mask)
        ann = np.zeros((ii.size, 2), dtype=np.int64)
        cre = np.zeros((ii.size, 2), dtype=np.int64)
        ann[ann_r, _within_rank(ann_r, ann_c)] = ann_c
        cre[cre_r, _within_rank(cre_r, cre_c)] = cre_c
        # _within_rank: 0 for the first (smaller) index of a pair, 1 for
        # the second; k=1 pairs use only slot 0 (slot 1 stays 0 = dummy)

        # closed-form transition signs from prefix counts (cum):
        #  k=1 (ann i, cre a): (-1)^{cum[j,i] + cum[j,a] - (i<a)}
        #  k=2 (ann i<j, cre a<b, b inserted first):
        #      (-1)^{cum[j,i] + cum[j,j_o] - 1 + cnt(i_or_j removed)(a|b)}
        sign = np.ones(ii.size)
        pref_j = cum[jj]                                      # (P, n2)
        m1 = k_c == 1
        if m1.any():
            i1, a1 = ann[m1, 0], cre[m1, 0]
            sign[m1] = (-1.0) ** (pref_j[m1, i1] + pref_j[m1, a1]
                                  - (i1 < a1).astype(np.int64))
        m2 = ~m1
        if m2.any():
            i2, j2 = ann[m2, 0], ann[m2, 1]                   # i2 < j2
            a2, b2 = cre[m2, 0], cre[m2, 1]                   # a2 < b2
            pj = pref_j[m2]
            ar = np.arange(len(i2))
            cnt = (lambda y: pj[ar, y]                        # noqa: E731
                   - (i2 < y).astype(np.int64) - (j2 < y).astype(np.int64))
            sign[m2] = (-1.0) ** (pj[ar, i2] + pj[ar, j2] - 1
                                  + cnt(b2) + cnt(a2))

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


def _within_rank(rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    """0 for an entry's first occurrence within its row, 1 for the
    second (np.nonzero returns ascending column order per row; rows
    have at most two entries because k <= 2)."""
    out = np.zeros(len(rows), dtype=np.int64)
    order = np.lexsort((cols, rows))
    rs = rows[order]
    rk = np.zeros(len(rs), dtype=np.int64)
    rk[1:] = (rs[1:] == rs[:-1]).astype(np.int64)  # 2nd within row -> 1
    out[order] = rk
    return out


def _slater_condon_reference(
    dets: list[tuple[int, ...]],
    one_body: np.ndarray | None = None,
    two_body: np.ndarray | None = None,
) -> np.ndarray:
    """Original loop implementation -- kept as the elementwise-identity
    reference for the vectorized :func:`slater_condon_matrix` (tests)."""
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
