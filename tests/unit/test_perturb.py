"""perturb.py rungs 0-1: the FOIS enumeration, the diagonal twin, the
EN construction twins, the degenerate anchors. Synthetic substrates
exercise every term class (random non-symmetric integrals); the
frozen carbon fixture pins the first real EN number -- NO SCF here.
"""

import numpy as np
import pytest

from helpers import random_integrals

from atomterm.detbasis import determinant_basis, slater_condon_matrix
from atomterm.perturb import (
    PT2Result,
    diagonal_energies,
    epstein_nesbet,
    epstein_nesbet_dense,
    fois_indices,
)

pytestmark = pytest.mark.fast


def test_en_vanishes_when_cas_is_whole_space():
    rng = np.random.default_rng(3)
    n = 3
    h1, g = random_integrals(rng, n, "hermitian")
    dets = determinant_basis(n, 1, 1)
    ci = rng.standard_normal(len(dets))
    ci /= np.linalg.norm(ci)
    res = epstein_nesbet(h1, g, n_act=n, na=1, nb=1, ci=ci)
    assert res.e2 == 0.0
    assert res.n_fois == 0
    assert res.floor == float("inf")
    assert res.w_ref == 1.0


def test_en_matches_dense_twin_on_random_integrals():
    # every term class: complex non-particle-symmetric integrals over
    # two sector/window shapes
    for seed, n_val, n_act, na, nb in ((11, 4, 2, 1, 1),
                                       (12, 4, 3, 2, 1),
                                       (13, 4, 3, 2, 0)):
        rng = np.random.default_rng(seed)
        h1, g = random_integrals(rng, n_val, "complex")
        dets_act = determinant_basis(n_act, na, nb)
        ci = rng.standard_normal(len(dets_act)) \
            + 1j * rng.standard_normal(len(dets_act))
        a = epstein_nesbet(h1, g, n_act, na, nb, ci)
        b = epstein_nesbet_dense(h1, g, n_act, na, nb, ci)
        assert isinstance(a, PT2Result) and isinstance(b, PT2Result)
        assert a.n_fois == b.n_fois > 0
        for f in ("e0", "e2", "floor", "w_ref"):
            assert np.allclose(getattr(a, f), getattr(b, f)), (f, seed)


def test_diagonal_energies_match_slater_condon_diagonal():
    # the occupation-sum diagonal vs the pair-classification diagonal,
    # on complex non-symmetric integrals and a mixed det list
    rng = np.random.default_rng(21)
    n = 4
    h1, g = random_integrals(rng, n, "complex")
    dets = determinant_basis(n, 2, 1)[::3]  # a strided det sublist
    ref = np.diag(slater_condon_matrix(dets, one_body=h1, two_body=g))
    got = diagonal_energies(dets, h1, g)
    assert np.allclose(got, ref)


def test_fois_cardinality_matches_brute_enumeration():
    # the production coupled pattern == the brute substitution count
    # on generic couplings (random hermitian integrals)
    rng = np.random.default_rng(31)
    n_val, n_act, na, nb = 4, 2, 1, 1
    h1, g = random_integrals(rng, n_val, "hermitian")
    dets_act = determinant_basis(n_act, na, nb)
    ci = rng.standard_normal(len(dets_act))
    ci /= np.linalg.norm(ci)
    res = epstein_nesbet(h1, g, n_act, na, nb, ci)
    brute = fois_indices(n_val, na, nb, n_act)
    assert res.n_fois == len(brute) > 0


def test_en_one_electron_closed_answer():
    # one electron, two valence orbitals, one active: every ingredient
    # hand-computable -- H_10 = h[1, 0] (no two-body channel), so
    # e2 = h10^2 / (h00 - h11) exactly
    h1 = np.array([[0.0, 0.1], [0.1, 1.0]])
    h1 = np.block([[h1, np.zeros((2, 2))], [np.zeros((2, 2)),
                                            np.zeros((2, 2))]])
    g = np.zeros((4, 4, 4, 4))
    res = epstein_nesbet(h1, g, n_act=1, na=1, nb=0,
                         ci=np.array([1.0]))
    assert res.e0 == pytest.approx(0.0)
    assert res.e2 == pytest.approx(0.1 ** 2 / (0.0 - 1.0))
    assert res.floor == pytest.approx(1.0)
    c1sq = (0.1 / (0.0 - 1.0)) ** 2
    assert res.w_ref == pytest.approx(1.0 / (1.0 + c1sq))


def test_en_raises_on_zero_denominator():
    # a degenerate partition must fail fast, not divide through: the
    # off-diagonal h1 gives det (1,) a NONZERO coupling while its
    # diagonal equals E0 = 0 -- a zero denominator with a live channel
    h1 = np.zeros((4, 4))
    h1[0, 1] = h1[1, 0] = 1.0
    g = np.zeros((4, 4, 4, 4))
    with pytest.raises(ValueError, match="zero EN denominator"):
        epstein_nesbet(h1, g, n_act=1, na=1, nb=0, ci=np.array([1.0]))


def test_en_carbon_frozen_fixture_pins_real_number(carbon):
    # the first real EN number: C/cc-pVDZ ROHF (frozen), valence =
    # the canonical MOs above the folded core (He . 2s), active = the
    # two open-shell MOs, sector (1,1); the reference CI is the CAS
    # ground. The pin is the OBSERVED fresh value (2026-11-18); the
    # twin agreement is the certification. The spatial tensors from
    # active_space_integrals are spin-blocked into the detbasis
    # convention (operators.spin_block_*).
    _, mol, mf, _ = carbon
    from atomterm.operators import (
        active_space_integrals,
        spin_block_one_body,
        spin_block_two_body,
    )

    ncore = (mol.nelectron - 2) // 2
    h1v_s, gv_s = active_space_integrals(mol, mf, mf.mo_coeff[:, 2:],
                                         ncore)
    h1v = spin_block_one_body(h1v_s)
    gv = spin_block_two_body(gv_s)

    n_act = 2
    na, nb = 1, 1
    dets_act = determinant_basis(n_act, na, nb)
    h1a_s, ga_s = active_space_integrals(mol, mf, mf.mo_coeff[:, 2:4],
                                         ncore)
    h_cas = slater_condon_matrix(dets_act, one_body=spin_block_one_body(h1a_s),
                                 two_body=spin_block_two_body(ga_s))
    evals, evecs = np.linalg.eigh(h_cas)
    ci = evecs[:, 0]

    a = epstein_nesbet(h1v, gv, n_act, na, nb, ci)
    b = epstein_nesbet_dense(h1v, gv, n_act, na, nb, ci)
    assert np.allclose(a.e2, b.e2)
    assert a.n_fois == b.n_fois
    # observed fresh (2026-11-18): the valence-space ruler (the core
    # energy constant is NOT in e0 -- the folded-field convention,
    # mc.e_cas-shaped); e2 = -7.689 mEh, floor 0.054, w_ref 0.996
    assert a.e0 == pytest.approx(-1.351391845908347)
    assert a.e2 == pytest.approx(-0.00768920099075394)
    assert a.floor == pytest.approx(0.054071241799460434)
    assert a.w_ref == pytest.approx(0.9959577555697562)
