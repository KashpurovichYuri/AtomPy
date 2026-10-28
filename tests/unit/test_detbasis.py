"""Tests of the pure-numpy determinant-basis machinery (no pyscf)."""

import numpy as np
import pytest

from atomterm.detbasis import (
    _pair_slot,
    _insert_sign,
    _pop_sign,
    _slater_condon_slow,
    _transition_sign,
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
from conftest import random_integrals

pytestmark = pytest.mark.fast


def _eigs(mat):
    return np.sort(np.real(np.linalg.eigvals(mat)))


def test_number_operator_counts_electrons():
    for norb, na, nb in [(3, 2, 0), (3, 1, 1), (4, 2, 1)]:
        dets = determinant_basis(norb, na, nb)
        n_mat = slater_condon_matrix(dets, one_body=number_operator(norb))
        assert np.allclose(np.diag(n_mat), na + nb)


@pytest.mark.parametrize(
    "na, nb, counts",
    [(2, 0, {2.0: 3}), (1, 1, {0.0: 6, 2.0: 3})],
)
def test_s2_spectrum_two_p_electrons(na, nb, counts):
    """2p^2 textbook: Ms=1 -> 3P only; Ms=0 -> 1S+1D + 3P."""
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


def test_matrix_hermitian_on_hermitian_random_integrals():
    rng = np.random.default_rng(7)
    norb, na, nb = 4, 2, 1
    h1, g = random_integrals(rng, norb, "hermitian")
    dets = determinant_basis(norb, na, nb)
    h = slater_condon_matrix(dets, one_body=h1, two_body=g)
    assert frobenius_norm(h - h.conj().T) < 1e-12


@pytest.mark.parametrize("na, nb", [(2, 0), (1, 1)])
def test_ladder_algebra_in_det_space(synthetic_p_operators, na, nb):
    """[Lx, Ly] = i Lz for the many-electron matrices."""
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


def test_vectorized_matches_slow_loop():
    rng = np.random.default_rng(11)
    norb = 5
    h, g = random_integrals(rng, norb, "complex")
    for na, nb in [(2, 0), (2, 1), (2, 2), (3, 1), (4, 0), (1, 0)]:
        dets = determinant_basis(norb, na, nb)
        for one, two in [(h, g), (h, None), (None, g)]:
            a = slater_condon_matrix(dets, one_body=one, two_body=two)
            b = _slater_condon_slow(dets, one_body=one, two_body=two)
            # <1e-13 = summation-reassociation noise; the slow loop
            # covers complex NON-symmetric two-body kernels (the L.L
            # composite class) and the beta-free sector whose n2 comes
            # from the operator, not the dets
            assert np.allclose(a, b, atol=1e-13), (na, nb, one is None)


def test_vectorized_matches_slow_loop_multichunk(monkeypatch):
    """The chunked pair loop's multi-chunk offsets (ii = rows + lo):
    the production memory budget yields ONE chunk at test sizes, so
    this test shrinks the budget to force eight chunks."""
    import atomterm.detbasis as db

    norb, na, nb = 4, 2, 1
    monkeypatch.setattr(db, "_PAIR_BUDGET", 3 * (2 * norb) * 24)  # chunk=3
    rng = np.random.default_rng(31)
    h, g = random_integrals(rng, norb, "complex")
    dets = determinant_basis(norb, na, nb)
    assert len(dets) == 24
    a = slater_condon_matrix(dets, one_body=h, two_body=g)
    b = _slater_condon_slow(dets, one_body=h, two_body=g)
    assert np.allclose(a, b, atol=1e-13)


# ------------------------------------------------------- the sign algebra

def test_elementary_sign_rules():
    # parity of the annihilated position
    assert _pop_sign([0, 1], 0) == 1
    assert _pop_sign([0, 1], 1) == -1
    assert _pop_sign([0, 1, 2], 2) == 1
    # parity of the creation position (# orbitals before it)
    assert _insert_sign([1, 2], 0) == 1
    assert _insert_sign([0, 2], 1) == -1
    # <I| a+_1 a_2 |J>, J = (0, 2): ann 2 at position 1 (-1), then
    # cre 1 with one orbital before it (-1) -> +1
    assert _transition_sign([0, 2], [2], [1]) == 1
    # double: J = (0, 1, 3), ann {1, 3}, cre {2, 4}: ann signs
    # (-1)(-1), then cre 4 (orbital 0 precedes it: -1), then cre 2
    # (orbital 0 precedes it: -1) -> +1
    assert _transition_sign([0, 1, 3], [1, 3], [2, 4]) == 1


def test_pair_slot_maps_difference_members():
    rows = np.array([0, 0, 1, 2, 2])
    cols = np.array([3, 7, 5, 2, 9])
    assert _pair_slot(rows, cols).tolist() == [0, 1, 0, 0, 1]
