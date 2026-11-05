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
* Two-body integrals, the physicist tensor:
  ``G[p, q, r, s] = <pq|rs> = Integral p(1) q(2) r(1) s(2)``
  (p, r ride electron 1, q, s electron 2; equal to chemist
  ``(pr|qs)``), shape ``(2*norb,)*4``; the operator is
  ``Sum_{i<j} g(i, j)``.

The determinant ordering here (lexicographic over the (alpha, beta)
combinations) is NOT PySCF's FCI string order -- they differ from
``norb >= 4``; the mapping lives in
``operators.pyscf_ci_permutation`` (this module stays NumPy-only).
This module is deliberately NumPy-only (no pyscf import), testable in
complete isolation.
"""

from __future__ import annotations

from itertools import combinations
from typing import NamedTuple
import bisect
import math

import numpy as np

__all__ = [
    "determinant_basis",
    "sector_dimension",
    "one_body_matrix",
    "two_body_matrix",
    "slater_condon_matrix",
    "operator_apply",
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


def sector_dimension(norb: int, na: int, nb: int) -> int:
    """The number of determinants of a sector: C(norb, na) * C(norb, nb).

    Equals ``len(determinant_basis(norb, na, nb))`` without building the
    basis -- the size budget every enumeration caller needs first.
    """
    return math.comb(norb, na) * math.comb(norb, nb)


def _pop_sign(occ: list[int], orb: int) -> int:
    """Sign of annihilating ``orb`` from the ordered occupation
    ``occ``: ``(-1)^(its position)``."""
    i = occ.index(orb)
    return -1 if i % 2 else 1


def _insert_sign(occ: list[int], orb: int) -> int:
    """Sign of creating ``orb`` into the ordered occupation ``occ``:
    ``(-1)^(# orbitals before it)``."""
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


# ------------------------------------------- the shared pair classification

# elements per pair chunk (the memory budget); tests shrink it to pin
# the multi-chunk offsets -- at test sizes the budget yields ONE chunk
_PAIR_BUDGET = 2e7


class _Pairs(NamedTuple):
    """One memory-bounded chunk of classified determinant pairs -- the
    machinery both rule-set evaluators consume."""

    ii: np.ndarray        # (P,) row determinants (absolute index)
    jj: np.ndarray        # (P,) column determinants
    m1: np.ndarray        # (P,) k == 1 mask
    ann: np.ndarray       # (P, 2) annihilated orbitals, ascending
    cre: np.ndarray       # (P, 2) created orbitals, ascending
    sign: np.ndarray      # (P,) closed-form transition signs
    occ_i: np.ndarray     # (P, n2) occupation rows of I
    occ_j: np.ndarray     # (P, n2) occupation rows of J


def _pair_chunks(occ: np.ndarray, cum: np.ndarray, nelec: int):
    """Yield the determinant pairs a two-body operator connects, in
    memory-bounded chunks: every pair with excitation rank ``k <= 2``
    (beyond double excitations the elements vanish), with per-pair
    excitation indices and closed-form transition signs.
    """
    n, n2 = occ.shape
    chunk = max(1, min(n, int(_PAIR_BUDGET) // (n * n2)))
    for lo in range(0, n, chunk):
        hi = min(n, lo + chunk)
        occ_i = occ[lo:hi]
        k = nelec - (occ_i[:, None, :] & occ[None, :, :]).sum(-1)
        rows, cols = np.nonzero((k == 1) | (k == 2))
        if rows.size == 0:
            continue
        ii = rows + lo
        jj = cols
        occ_j_c = occ[jj]
        occ_i_c = occ[ii]
        ann_mask = occ_j_c & ~occ_i_c
        cre_mask = occ_i_c & ~occ_j_c
        k_c = k[rows, cols]

        # per-pair excitation indices, ascending within each set;
        # _pair_slot: 0 = smaller member, 1 = larger (k=1 leaves the
        # second slot as a dummy)
        ann_r, ann_c = np.nonzero(ann_mask)
        cre_r, cre_c = np.nonzero(cre_mask)
        ann = np.zeros((ii.size, 2), dtype=np.int64)
        cre = np.zeros((ii.size, 2), dtype=np.int64)
        ann[ann_r, _pair_slot(ann_r, ann_c)] = ann_c
        cre[cre_r, _pair_slot(cre_r, cre_c)] = cre_c

        sign = np.ones(ii.size)
        m1 = k_c == 1
        if m1.any():
            sign[m1] = _sign_single(cum[jj][m1], ann[m1, 0], cre[m1, 0])
        m2 = ~m1
        if m2.any():
            sign[m2] = _sign_double(cum[jj][m2], ann[m2, 0], ann[m2, 1],
                                    cre[m2, 0], cre[m2, 1])

        yield _Pairs(ii, jj, m1, ann, cre, sign, occ_i_c, occ_j_c)


# ------------------------------------------------------- the two rule sets

def _one_body_pairs(pairs: _Pairs, h: np.ndarray):
    """The SIGNED one-body couplings of a pair chunk as
    ``(rows, cols, values)`` -- ``<I| h[a, i] a+_a a_i |J>`` on the
    k=1 pairs.  The one home for the rule: the matrix builder and the
    matrix-free application consume the same triples."""
    if not pairs.m1.any():
        z = np.zeros(0, dtype=int)
        return z, z, np.zeros(0)
    val = h[pairs.cre[pairs.m1, 0], pairs.ann[pairs.m1, 0]]
    nz = val != 0
    return (pairs.ii[pairs.m1][nz], pairs.jj[pairs.m1][nz],
            pairs.sign[pairs.m1][nz] * val[nz])


def _two_body_pairs(pairs: _Pairs, g: np.ndarray):
    """The SIGNED two-body couplings of a pair chunk as
    ``(rows, cols, values)``, both off-diagonal rules (chemists' ``G``
    throughout) -- the one home for the rules:

    * k=1: the contraction over the orbitals common to I and J (the
      four-term antisymmetrized gather);
    * k=2: the direct-exchange pair value.
    """
    n2 = pairs.occ_i.shape[1]
    orb = np.arange(n2)
    iis, jjs, vals = [], [], []
    if pairs.m1.any():
        i1, a1 = pairs.ann[pairs.m1, 0], pairs.cre[pairs.m1, 0]
        common = pairs.occ_j[pairs.m1] & pairs.occ_i[pairs.m1]
        tb = 0.5 * (
            g[a1[:, None], orb[None, :], i1[:, None], orb[None, :]]
            + g[orb[None, :], a1[:, None], orb[None, :], i1[:, None]]
            - g[orb[None, :], a1[:, None], i1[:, None], orb[None, :]]
            - g[a1[:, None], orb[None, :], orb[None, :], i1[:, None]]
        )
        val = (tb * common).sum(axis=1)
        nz = val != 0
        iis.append(pairs.ii[pairs.m1][nz])
        jjs.append(pairs.jj[pairs.m1][nz])
        vals.append(pairs.sign[pairs.m1][nz] * val[nz])
    m2 = ~pairs.m1
    if m2.any():
        i2, j2 = pairs.ann[m2, 0], pairs.ann[m2, 1]
        a2, b2 = pairs.cre[m2, 0], pairs.cre[m2, 1]
        val = 0.5 * (
            g[a2, b2, i2, j2] + g[b2, a2, j2, i2]
            - g[b2, a2, i2, j2] - g[a2, b2, j2, i2]
        )
        nz = val != 0
        iis.append(pairs.ii[m2][nz])
        jjs.append(pairs.jj[m2][nz])
        vals.append(pairs.sign[m2][nz] * val[nz])
    if not iis:
        z = np.zeros(0, dtype=int)
        return z, z, np.zeros(0)
    return np.concatenate(iis), np.concatenate(jjs), np.concatenate(vals)


def _diag_one_body(dets, occ, h, dtype) -> np.ndarray:
    """Diagonal one-body elements: ``Sum_p <p|h|p>`` per determinant."""
    return occ @ np.diag(h).astype(dtype)


def _diag_two_body(dets, g, dtype) -> np.ndarray:
    """Diagonal two-body elements: the antisymmetrized sum over all
    occupied pairs ``1/2 (G[pqpq] + G[qpqp] - G[qppq] - G[pqqp])``."""
    nelec = len(dets[0])
    if nelec < 2:
        return np.zeros(len(dets), dtype=dtype)
    occ_idx = np.array([list(d) for d in dets])          # (n, nelec)
    comb = np.array(list(combinations(range(nelec), 2)))
    p = occ_idx[:, comb[:, 0]]                            # (n, nc)
    q = occ_idx[:, comb[:, 1]]
    return 0.5 * (
        g[p, q, p, q] + g[q, p, q, p]
        - g[q, p, p, q] - g[p, q, q, p]
    ).sum(axis=1)


# ----------------------------------------------------------- the public API

def one_body_matrix(
    dets: list[tuple[int, ...]], one_body: np.ndarray
) -> np.ndarray:
    """The matrix of ``O1 = Sum_pq h[p,q] a+_p a_q`` over ``dets`` --
    the one-body rule set alone (the name-twin of
    :func:`atomterm.secondquant.one_body_matrix`).

    Parameters
    ----------
    dets : determinant basis (see :func:`determinant_basis`).
    one_body : ``(2*norb, 2*norb)`` matrix ``h[p, q]``.

    Returns
    -------
    ``(len(dets), len(dets))`` matrix ``M[I, J] = <I|O1|J>``; dtype
    follows the input (real ``h`` stays real).
    """
    return _build(dets, one_body, None)


def two_body_matrix(
    dets: list[tuple[int, ...]], two_body: np.ndarray
) -> np.ndarray:
    """The matrix of ``O2 = 1/2 Sum_pqrs G[p,q,r,s] a+_p a+_q a_s a_r``
    over ``dets``, ``G`` in the module's chemist notation -- the
    two-body rule set alone (the name-twin of
    :func:`atomterm.secondquant.two_body_matrix`).

    Parameters
    ----------
    dets : determinant basis (see :func:`determinant_basis`).
    two_body : ``(2*norb,)*4`` chemist-notation integrals.

    Returns
    -------
    ``(len(dets), len(dets))`` matrix; dtype follows the input.
    """
    return _build(dets, None, two_body)


def slater_condon_matrix(
    dets: list[tuple[int, ...]],
    one_body: np.ndarray | None = None,
    two_body: np.ndarray | None = None,
) -> np.ndarray:
    """The full matrix of ``Sum_i h(i) + Sum_{i<j} g(i, j)`` over
    ``dets``: both rule sets scattered over ONE pair classification --
    the production path (the 10^5-detector sectors;
    :func:`_slater_condon_slow`, the O(ndet^2) Python loop it
    replaced, stays as the elementwise reference pinned by tests).

    Parameters
    ----------
    dets : determinant basis (see :func:`determinant_basis`).
    one_body : ``(2*norb, 2*norb)`` matrix ``h[p, q]`` (optional).
    two_body : ``(2*norb,)*4`` chemist-notation integrals (optional).

    Returns
    -------
    ``(len(dets), len(dets))`` matrix, ``M[I, J] = <I|O|J>``; dtype
    follows the inputs (real operators stay real -- a complex
    container would double memory traffic and slow the
    eigendecomposition). For the parts separately:
    :func:`one_body_matrix` / :func:`two_body_matrix`.

    Raises
    ------
    ValueError
        If both operator parts are ``None``.
    """
    if one_body is None and two_body is None:
        raise ValueError("at least one of one_body / two_body is required")
    return _build(dets, one_body, two_body)


def _singles_couplings(dets, occ: np.ndarray, cum: np.ndarray):
    """All SINGLE-excitation couplings as ``(rows, cols, ann, cre)`` --
    the one-body operator's entire nonzero pattern, enumerated by
    index arithmetic (bitmask searchsorted) instead of the pairwise
    classification: O(n * nelec * nvirt), no n^2 work.  The one-body
    closure at big sectors rides this; the pair classifier serves the
    two-body class.
    """
    n, n2 = occ.shape
    nelec = len(dets[0])
    nvirt = n2 - nelec
    if nelec == 0 or nvirt == 0:
        z = np.zeros(0, dtype=np.int64)
        return z, z, z, z
    if n2 > 63:
        raise ValueError(
            f"the bitmask enumeration addresses {n2} spin-orbitals in "
            "int64 -- the safe limit is 63; wider sectors route to "
            "_singles_couplings_wide (the packed-key lookup), not a "
            "wider mask")
    masks = (occ * (np.int64(1) << np.arange(n2))).sum(axis=1)
    order = np.argsort(masks)
    sorted_masks = masks[order]
    occ_of = np.nonzero(occ)[1].reshape(n, nelec)
    virt_of = np.nonzero(~occ)[1].reshape(n, nvirt)
    ii = np.repeat(np.arange(n, dtype=np.int64), nelec * nvirt)
    ann = np.repeat(occ_of, nvirt, axis=1).ravel()
    cre = np.tile(virt_of, (1, nelec)).ravel()
    j_mask = masks[ii] ^ (np.int64(1) << ann) ^ (np.int64(1) << cre)
    pos = np.searchsorted(sorted_masks, j_mask)
    np.clip(pos, 0, n - 1, out=pos)
    valid = sorted_masks[pos] == j_mask
    ii, ann, cre = ii[valid], ann[valid], cre[valid]
    jj = order[pos[valid]]
    return ii, jj, ann, cre


def _singles_couplings_wide(dets, occ: np.ndarray):
    """The singles enumeration past int64 (``n2 > 63``): the SAME
    ``O(n * nelec * nvirt)`` triples as :func:`_singles_couplings`,
    with the target lookup on packed occupation rows (fixed-width
    byte keys, lexicographic searchsorted) instead of int64 bitmask
    arithmetic -- width-unlimited in the spin-orbital count, the
    triples table the next ceiling.  The signs ride the prefix counts
    (:func:`_sign_single`), the same rule as the bitmask path; the
    validity filter stays (a target must EXIST -- the sector's
    determinant table is complete, so a miss is a bug, never a
    silent mis-map).
    """
    n, n2 = occ.shape
    nelec = len(dets[0])
    nvirt = n2 - nelec
    if nelec == 0 or nvirt == 0:
        z = np.zeros(0, dtype=np.int64)
        return z, z, z, z
    occ_of = np.nonzero(occ)[1].reshape(n, nelec)
    virt_of = np.nonzero(~occ)[1].reshape(n, nvirt)
    ii = np.repeat(np.arange(n, dtype=np.int64), nelec * nvirt)
    ann = np.repeat(occ_of, nvirt, axis=1).ravel()
    cre = np.tile(virt_of, (1, nelec)).ravel()
    # the target rows: occ[ii] with the annihilated bit cleared and the
    # created bit set (disjoint by construction), packed straight to
    # the byte keys -- a uint8 gather, never the bool table
    src = np.packbits(occ, axis=1, bitorder="little")
    k = src.shape[1]
    rows = src[ii]
    shift_a = (ann & 7).astype(np.uint8)
    shift_c = (cre & 7).astype(np.uint8)
    rows[np.arange(ii.size), ann >> 3] &= ~np.left_shift(
        np.uint8(1), shift_a)
    rows[np.arange(ii.size), cre >> 3] |= np.left_shift(np.uint8(1),
                                                        shift_c)
    dkeys = src.view(f"S{k}").ravel()
    order = np.argsort(dkeys)
    sorted_keys = dkeys[order]
    query = rows.view(f"S{k}").ravel()
    pos = np.clip(np.searchsorted(sorted_keys, query), 0, n - 1)
    valid = sorted_keys[pos] == query
    ii, ann, cre = ii[valid], ann[valid], cre[valid]
    jj = order[pos[valid]]
    return ii, jj, ann, cre


def operator_apply(
    dets: list[tuple[int, ...]],
    one_body: np.ndarray | None = None,
    two_body: np.ndarray | None = None,
):
    """The MATRIX-FREE application of a (one-body, two-body) operator
    pair: returns ``apply(v) -> O @ v`` -- the operator exists only
    through its action, the ``len(dets)^2`` container never forms.

    The one-body class rides the singles enumeration: bitmask
    searchsorted at up to 63 spin-orbitals (int64), the packed-key
    lookup past that (:func:`_singles_couplings_wide`) -- the same
    signed triples either way, wider sectors never silently
    truncated.

    The nonzero pattern is exact, not heuristic: an operator of this
    form connects determinants differing by at most two spin-orbitals
    (beyond that the strings vanish by Pauli, or the result is
    orthogonal to the bra) -- :func:`_pair_chunks` enumerates exactly
    those pairs, and the SAME signed couplings
    (:func:`_one_body_pairs` / :func:`_two_body_pairs`) that build the
    matrix (:func:`slater_condon_matrix`) accumulate the output
    vector here, so ``apply(v) == slater_condon_matrix(...) @ v`` is
    an identity, certifiable elementwise.

    Cost: the FIRST application pays the pair classification once and
    caches the signed couplings (~16 bytes per nonzero -- budget the
    sector); every application after is O(#nonzeros), no quadratic
    work and no dense operator in memory.  For sectors whose cached
    couplings would not fit, the PySCF contraction boundary
    (:func:`atomterm.operators.pyscf_apply`) applies the same operator
    pair at C speed with no caching.

    Parameters
    ----------
    dets, one_body, two_body
        As :func:`slater_condon_matrix` (at least one part required).

    Returns
    -------
    ``apply`` -- a closure mapping a sector vector (real or complex)
    to ``O @ v``.

    Examples
    --------
    >>> dets = determinant_basis(3, 2, 1)
    >>> h = np.eye(6)
    >>> f = operator_apply(dets, one_body=h)
    >>> v = np.zeros(len(dets)); v[0] = 1.0
    >>> np.allclose(f(v), slater_condon_matrix(dets, one_body=h) @ v)
    True
    """
    if one_body is None and two_body is None:
        raise ValueError("at least one of one_body / two_body is required")
    n = len(dets)
    n2 = (one_body if one_body is not None else two_body).shape[0]
    occ = np.zeros((n, n2), dtype=bool)
    for i, d in enumerate(dets):
        occ[i, list(d)] = True
    # prefix counts: cum[j, t] = # occupied spin-orbitals < t in det j
    cum = np.zeros((n, n2), dtype=np.int64)
    np.cumsum(occ[:, :-1], axis=1, out=cum[:, 1:])
    dtype = complex if (np.iscomplexobj(one_body if one_body is not None
                                        else two_body)) else float
    diag = np.zeros(n, dtype=dtype)
    if one_body is not None:
        diag += _diag_one_body(dets, occ, one_body, dtype)
    if two_body is not None:
        diag += _diag_two_body(dets, two_body, dtype)

    couplings = None  # the signed (rows, cols, values) chunks, built once

    def _build_couplings():
        # one-body-only operators ride the singles enumeration (no
        # pairwise classification -- the big-sector one-body path);
        # anything with a two-body part uses the chunked classifier
        out = []
        if two_body is None:
            n2 = occ.shape[1]
            pick = (_singles_couplings_wide(dets, occ) if n2 > 63
                    else _singles_couplings(dets, occ, cum))
            ii, jj, ann, cre = pick
            # ii enumerates the SOURCE determinants (the kets the
            # operator acts on), jj their single-excitation partners;
            # the sign is the source's transition sign, and the
            # scatter accumulates into the PARTNER rows
            sign = (_sign_single(cum[ii], ann, cre)
                    if ii.size else np.zeros(0))
            val = sign * one_body[cre, ann]
            nz = val != 0
            out.append((jj[nz], ii[nz], val[nz]))
            return out
        for pairs in _pair_chunks(occ, cum, len(dets[0])):
            rows_i, cols_i, vals_i = [], [], []
            if one_body is not None:
                ii, jj, val = _one_body_pairs(pairs, one_body)
                rows_i.append(ii); cols_i.append(jj); vals_i.append(val)
            ii, jj, val = _two_body_pairs(pairs, two_body)
            rows_i.append(ii); cols_i.append(jj); vals_i.append(val)
            out.append((np.concatenate(rows_i),
                        np.concatenate(cols_i),
                        np.concatenate(vals_i)))
        return out

    def apply(v):
        nonlocal couplings
        v = np.asarray(v)
        out_dtype = np.result_type(v.dtype, dtype)
        if couplings is None:
            couplings = _build_couplings()
        out = diag * v.astype(out_dtype)
        for ii, jj, val in couplings:
            prod = val * v[jj]
            if np.iscomplexobj(prod):
                out += np.bincount(ii, weights=prod.real, minlength=n)
                out += 1j * np.bincount(ii, weights=prod.imag, minlength=n)
            else:
                out += np.bincount(ii, weights=prod, minlength=n)
        return out

    return apply


def _build(
    dets: list[tuple[int, ...]],
    one_body: np.ndarray | None,
    two_body: np.ndarray | None,
) -> np.ndarray:
    """Assemble the matrix: occupation table and prefix counts once,
    diagonals from the active rule sets, then ONE chunked pair
    classification consumed by both off-diagonal evaluators."""
    n = len(dets)
    # the operator fixes the spin-orbital dimension: sectors without
    # beta electrons never index the beta block (max det index < norb)
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

    diag = np.zeros(n, dtype=dtype)
    if one_body is not None:
        diag += _diag_one_body(dets, occ, one_body, dtype)
    if two_body is not None:
        diag += _diag_two_body(dets, two_body, dtype)
    mat[np.diag_indices(n)] = diag

    for pairs in _pair_chunks(occ, cum, len(dets[0])):
        if one_body is not None:
            ii, jj, val = _one_body_pairs(pairs, one_body)
            mat[ii, jj] += val
        if two_body is not None:
            ii, jj, val = _two_body_pairs(pairs, two_body)
            mat[ii, jj] += val
    return mat


# ------------------------------------------------------ the sign machinery

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
    reference for the vectorized path (pinned by tests; sign engine:
    :func:`_transition_sign`)."""
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
