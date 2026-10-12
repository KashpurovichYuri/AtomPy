"""Tests of the pure-numpy determinant-basis machinery (no pyscf)."""

import numpy as np
import pytest

from atomterm.detbasis import (
    commutator,
    determinant_basis,
    frobenius_norm,
    slater_condon_matrix,
)
from atomterm.operators import (
    _composite_two_body,
    number_operator,
    s2_operator_parts,
    spin_block_one_body,
    spin_block_two_body,
)


def _eigs(mat):
    return np.sort(np.real(np.linalg.eigvals(mat)))


def test_number_operator():
    for norb, na, nb in [(3, 2, 0), (3, 1, 1), (4, 2, 1)]:
        dets = determinant_basis(norb, na, nb)
        n_mat = slater_condon_matrix(dets, one_body=number_operator(norb))
        assert np.allclose(np.diag(n_mat), na + nb)


@pytest.mark.parametrize(
    "na, nb, counts",
    [(2, 0, {2.0: 3}), (1, 1, {0.0: 6, 2.0: 3})],
)
def test_s2_spectrum_two_p_electrons(na, nb, counts):
    """2p^2 textbook: Ms=1 -> 3P only (S=1); Ms=0 -> 1S+1D (S=0) + 3P (S=1)."""
    dets = determinant_basis(3, na, nb)
    one, two = s2_operator_parts(3)
    s2 = slater_condon_matrix(dets, one_body=one, two_body=two)
    assert frobenius_norm(s2 - s2.conj().T) < 1e-12
    vals = _eigs(s2)
    for target, multiplicity in counts.items():
        assert np.sum(np.abs(vals - target) < 1e-9) == multiplicity, (vals, target)


def test_sz_diagonal_matches_ms():
    norb = 3
    for na, nb in [(2, 0), (1, 1)]:
        dets = determinant_basis(norb, na, nb)
        sz_so = np.zeros((2 * norb, 2 * norb), dtype=complex)
        sz_so[:norb, :norb] = 0.5 * np.eye(norb)
        sz_so[norb:, norb:] = -0.5 * np.eye(norb)
        sz = slater_condon_matrix(dets, one_body=sz_so)
        assert np.allclose(np.diag(sz), (na - nb) / 2)


def test_hermiticity_random_integrals():
    rng = np.random.default_rng(7)
    norb, na, nb = 4, 2, 1
    a = rng.normal(size=(2 * norb, 2 * norb))
    h1 = 0.5 * (a + a.T)
    g = rng.normal(size=(2 * norb,) * 4)
    # Hermitian two-body operator: <pq|g|rs>* = <rs|g|pq>
    g = 0.5 * (g + g.transpose(2, 3, 0, 1).conj())
    dets = determinant_basis(norb, na, nb)
    h = slater_condon_matrix(dets, one_body=h1, two_body=g)
    assert frobenius_norm(h - h.conj().T) < 1e-12


@pytest.mark.parametrize("na, nb", [(2, 0), (1, 1)])
def test_ladder_algebra_in_det_space(synthetic_p_operators, na, nb):
    """[Lx, Ly] = i Lz must hold for the many-electron matrices as well."""
    l_ops = synthetic_p_operators
    dets = determinant_basis(3, na, nb)
    mats = {
        a: slater_condon_matrix(dets, one_body=spin_block_one_body(l))
        for a, l in l_ops.items()
    }
    comm = commutator(mats["x"], mats["y"]) - 1j * mats["z"]
    assert frobenius_norm(comm) < 1e-12


def test_composite_two_body_hermitian():
    """Composite integrals 2 (Lz x Lz) from an imaginary Hermitian Lz."""
    i = 1j
    lz = np.array([[0, -i, 0], [i, 0, 0], [0, 0, 0]])
    two = 2.0 * _composite_two_body(lz, lz)
    assert frobenius_norm(two - two.conj().T) < 1e-12
    one = lz @ lz
    dets = determinant_basis(3, 2, 0)
    m = slater_condon_matrix(
        dets,
        one_body=spin_block_one_body(one),
        two_body=spin_block_two_body(two),
    )
    assert frobenius_norm(m - m.conj().T) < 1e-12
