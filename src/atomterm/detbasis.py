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
    """
    if one_body is None and two_body is None:
        raise ValueError("at least one of one_body / two_body is required")
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
