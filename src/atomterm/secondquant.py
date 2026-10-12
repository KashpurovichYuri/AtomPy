"""Second-quantized operator construction via explicit a+/a application.

This is the second independent construction of the many-electron operator
matrices (the first being the Slater-Condon rules in atomterm.detbasis).
Operators are fermionic strings acting on occupation-number vectors; matrix
elements are obtained by *applying the operator to each basis ket once* and
scattering the results -- a different algorithmic path entirely.

Conventions (identical to atomterm.detbasis)
--------------------------------------------
* Spin-orbitals: alpha = 0..norb-1, beta = norb..2norb-1.
* Determinants: sorted tuples of occupied spin-orbitals.
* Fock state: |occ> = prod_{p in occ, ascending application} a+_p |0>,
  i.e. the lowest orbital sits rightmost.
* One-body operator:  O1 = sum_pq h[p,q] a+_p a_q.
* Two-body operator:  O2 = 1/2 sum_pqrs G[p,q,r,s] a+_p a+_q a_s a_r,
  with G in the detbasis convention G[p,q,r,s] = <pq|rs> =
  Integral p(1) q(2) r(1) s(2)  (physicist arrangement = chemist (pr|qs)).

With this convention the double-excitation element between determinants
I = {p,q}, J = {r,s} evaluates to

    <I|O2|J> = 1/2 (G[pqrs] - G[pqsr] - G[qprs] + G[qpsr])

which reduces to the textbook `G[pqrs] - G[pqsr]` for particle-symmetric
kernels (all composites used in atomterm are particle-symmetric). The string
application below evaluates the full expression exactly -- including the
one-body-like contraction terms a+_p n_q a_r inside the strings -- with no
additional approximations. Correctness is pinned elementwise against the
first-quantized Slater-Condon matrices in tests/test_secondquant.py.

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
    "secondquant_matrix",
]


def annihilate(occ: tuple[int, ...], p: int):
    """a_p |occ>  ->  (sign, occ') or None if p is not occupied."""
    occ = list(occ)
    try:
        i = occ.index(p)
    except ValueError:
        return None
    del occ[i]
    return (-1 if i % 2 else 1), tuple(occ)


def create(occ: tuple[int, ...], p: int):
    """a+_p |occ>  ->  (sign, occ') or None if p is already occupied."""
    if p in occ:
        return None
    i = bisect.bisect_left(occ, p)
    occ = list(occ)
    occ.insert(i, p)
    return (-1 if i % 2 else 1), tuple(occ)


def apply_one_body(occ: tuple[int, ...], h: np.ndarray) -> dict[tuple, complex]:
    """(sum_pq h[p,q] a+_p a_q) |occ>  ->  {occ': amplitude}."""
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
    """(1/2 sum_pqrs G[p,q,r,s] a+_p a+_q a_s a_r) |occ>  ->  {occ': amp}.

    `nz` optionally carries the nonzero indices of G (from np.argwhere) to
    skip structurally zero terms.
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


def secondquant_matrix(
    dets: list[tuple[int, ...]],
    one_body: np.ndarray | None = None,
    two_body: np.ndarray | None = None,
) -> np.ndarray:
    """Full many-electron matrix of O1 + O2 in the determinant basis.

    Column J holds O|det_J>, accumulated by explicit operator application.
    States outside ``dets`` are dropped, so pass a complete basis for the full
    matrix. Output is elementwise identical (up to fp noise) to
    ``atomterm.detbasis.slater_condon_matrix`` on the same inputs.
    """
    if one_body is None and two_body is None:
        raise ValueError("at least one of one_body / two_body is required")
    index = {det: i for i, det in enumerate(dets)}
    mat = np.zeros((len(dets), len(dets)), dtype=complex)
    nz = np.argwhere(two_body != 0) if two_body is not None else None
    for j, occ in enumerate(dets):
        if one_body is not None:
            for target, amp in apply_one_body(occ, one_body).items():
                i = index.get(target)
                if i is not None:  # off-basis states are dropped; pass a
                    mat[i, j] += amp  # complete basis for the full matrix
        if two_body is not None:
            for target, amp in apply_two_body(occ, two_body, nz).items():
                i = index.get(target)
                if i is not None:
                    mat[i, j] += amp
    return mat
