"""Second-quantized operator construction via explicit a+/a
application -- the second, independent construction of the
many-electron operator matrices (the first being the Slater-Condon
rules in :mod:`atomterm.detbasis`). Matrix elements are obtained by
*applying the operator to each basis ket once* and scattering the
results: a different algorithmic path entirely, pinned elementwise
against the first construction in tests.

Conventions (identical to :mod:`atomterm.detbasis`)
---------------------------------------------------
* Spin-orbitals: alpha = 0..norb-1, beta = norb..2norb-1.
* Determinants: sorted tuples of occupied spin-orbitals.
* Fock state: |occ> = prod_{p in occ, ascending application}
  a+_p |0>, i.e. the lowest orbital sits rightmost.
* One-body operator:  O1 = sum_pq h[p,q] a+_p a_q.
* Two-body operator:  O2 = 1/2 sum_pqrs G[p,q,r,s] a+_p a+_q a_s a_r,
  with G in the detbasis convention G[p,q,r,s] = <pq|rs> =
  Integral p(1) q(2) r(1) s(2)  (physicist arrangement = chemist
  (pr|qs)). Between determinants I = {p,q}, J = {r,s} the double
  excitation element is
  ``1/2 (G[pqrs] - G[pqsr] - G[qprs] + G[qpsr])``; the string
  application evaluates this exactly, including the contraction terms
  inside the strings -- random non-particle-symmetric G is part of
  the elementwise pin.

This module is deliberately NumPy-only (no pyscf import).
"""

from __future__ import annotations

import bisect

import numpy as np

__all__ = [
    "annihilate",
    "create",
    "apply_one_body",
    "apply_two_body",
    "one_body_matrix",
    "two_body_matrix",
    "secondquant_matrix",
]


# the elementary operators -------------------------------------------

def annihilate(occ: tuple[int, ...], p: int):
    """``a_p |occ>`` -> ``(sign, occ')`` or None if p is not occupied."""
    occ = list(occ)
    try:
        i = occ.index(p)
    except ValueError:
        return None
    del occ[i]
    return (-1 if i % 2 else 1), tuple(occ)


def create(occ: tuple[int, ...], p: int):
    """``a+_p |occ>`` -> ``(sign, occ')`` or None if p is already
    occupied."""
    if p in occ:
        return None
    i = bisect.bisect_left(occ, p)
    occ = list(occ)
    occ.insert(i, p)
    return (-1 if i % 2 else 1), tuple(occ)


# the ket applications -----------------------------------------------

def apply_one_body(occ: tuple[int, ...], h: np.ndarray) -> dict[tuple, complex]:
    """``(sum_pq h[p,q] a+_p a_q) |occ>`` -> ``{occ': amplitude}``."""
    out: dict[tuple, complex] = {}
    for q in occ:
        killed = annihilate(occ, q)
        if killed is None:  # pragma: no cover - q is occupied by construction
            continue
        sg_q, o1 = killed
        col = h[:, q]
        for p in np.nonzero(col)[0]:
            made = create(o1, int(p))
            if made is None:  # p was occupied: a+_p a_q gives zero
                continue
            sg_p, o2 = made
            amp = sg_q * sg_p * h[p, q]
            out[o2] = out.get(o2, 0j) + amp
    return out


def apply_two_body(
    occ: tuple[int, ...], g: np.ndarray, nz: np.ndarray | None = None
) -> dict[tuple, complex]:
    """``(1/2 sum_pqrs G[p,q,r,s] a+_p a+_q a_s a_r) |occ>`` ->
    ``{occ': amplitude}``.

    Parameters
    ----------
    occ : the ket's occupation.
    g : the two-body integrals (module convention).
    nz : optional precomputed nonzero indices of ``g`` (from
        ``np.argwhere``) to skip structurally zero terms.
    """
    out: dict[tuple, complex] = {}
    terms = np.argwhere(g != 0) if nz is None else nz
    for p, q, r, s in terms:
        p, q, r, s = int(p), int(q), int(r), int(s)
        kill_r = annihilate(occ, r)  # rightmost operator acts first
        if kill_r is None:
            continue
        sg_r, o1 = kill_r
        kill_s = annihilate(o1, s)
        if kill_s is None:
            continue
        sg_s, o2 = kill_s
        make_q = create(o2, q)
        if make_q is None:
            continue
        sg_q, o3 = make_q
        make_p = create(o3, p)
        if make_p is None:
            continue
        sg_p, o4 = make_p
        amp = 0.5 * sg_r * sg_s * sg_q * sg_p * g[p, q, r, s]
        out[o4] = out.get(o4, 0j) + amp
    return out


# the matrices -------------------------------------------------------

def one_body_matrix(
    dets: list[tuple[int, ...]], one_body: np.ndarray
) -> np.ndarray:
    """Matrix of ``O1 = sum_pq h[p,q] a+_p a_q`` over ``dets``.

    Parameters
    ----------
    dets : a COMPLETE determinant basis (states outside it are
        dropped from the scattering).
    one_body : ``(2*norb, 2*norb)`` matrix ``h[p, q]``.

    Returns
    -------
    ``(len(dets), len(dets))`` complex matrix; column J holds
    ``O1|det_J>`` scattered onto the basis. The name-twin of
    :func:`atomterm.detbasis.one_body_matrix` (elementwise equal).
    """
    index = {det: i for i, det in enumerate(dets)}
    mat = np.zeros((len(dets), len(dets)), dtype=complex)
    for j, occ in enumerate(dets):
        for target, amp in apply_one_body(occ, one_body).items():
            i = index.get(target)
            if i is not None:
                mat[i, j] += amp
    return mat


def two_body_matrix(
    dets: list[tuple[int, ...]], two_body: np.ndarray
) -> np.ndarray:
    """Matrix of ``O2 = 1/2 sum_pqrs G[p,q,r,s] a+_p a+_q a_s a_r``
    over ``dets``.

    Parameters
    ----------
    dets : a COMPLETE determinant basis (same scattering contract as
        :func:`one_body_matrix`).
    two_body : ``(2*norb,)*4`` integrals in the module convention.

    Returns
    -------
    ``(len(dets), len(dets))`` complex matrix. The name-twin of
    :func:`atomterm.detbasis.two_body_matrix` (elementwise equal).
    """
    index = {det: i for i, det in enumerate(dets)}
    mat = np.zeros((len(dets), len(dets)), dtype=complex)
    nz = np.argwhere(two_body != 0)
    for j, occ in enumerate(dets):
        for target, amp in apply_two_body(occ, two_body, nz).items():
            i = index.get(target)
            if i is not None:
                mat[i, j] += amp
    return mat


def secondquant_matrix(
    dets: list[tuple[int, ...]],
    one_body: np.ndarray | None = None,
    two_body: np.ndarray | None = None,
) -> np.ndarray:
    """Full many-electron matrix of ``O1 + O2`` in the determinant
    basis: :func:`one_body_matrix` + :func:`two_body_matrix`.
    Elementwise identical (to fp noise) to
    :func:`atomterm.detbasis.slater_condon_matrix` on the same inputs.
    """
    if one_body is None and two_body is None:
        raise ValueError("at least one of one_body / two_body is required")
    if one_body is None:
        return two_body_matrix(dets, two_body)
    if two_body is None:
        return one_body_matrix(dets, one_body)
    return one_body_matrix(dets, one_body) + two_body_matrix(dets, two_body)
