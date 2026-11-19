"""perturb.py rungs 0-1: the FOIS enumeration, the diagonal twin, the
EN construction twins, the degenerate anchors. Synthetic substrates
exercise every term class (random non-symmetric integrals); the
frozen carbon fixture pins the first real EN number -- NO SCF here.
"""

import numpy as np
import pytest
from itertools import combinations

from helpers import random_integrals

from atomterm.detbasis import determinant_basis, slater_condon_matrix
from atomterm.perturb import (
    PT2Result,
    diagonal_energies,
    dyall_e2,
    dyall_e2_dense,
    dyall_h0,
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
    with pytest.raises(ValueError, match="zero PT2 denominator"):
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


# ------------------------------------------------------- rung 2: Dyall

def _physical_integrals(rng, norb: int):
    """hermitian + particle-symmetric random integrals: the Dyall
    Fock's input class -- the chemist 8-fold symmetry projected by
    averaging over its generated group (the detbasis physicist
    generators: 0<->2, 1<->3, and the pair swap (0<->1)(2<->3))."""
    h, g = random_integrals(rng, norb, "hermitian")
    gens = [(2, 1, 0, 3), (0, 3, 2, 1), (1, 0, 3, 2)]
    seen = {(0, 1, 2, 3)}
    frontier = [(0, 1, 2, 3)]
    while frontier:
        cur = frontier.pop()
        for gen in gens:
            nxt = tuple(cur[j] for j in gen)
            if nxt not in seen:
                seen.add(nxt)
                frontier.append(nxt)
    assert len(seen) == 8  # the chemist 8-fold, closed
    return h, sum(g.transpose(p) for p in seen) / len(seen)


def _canonical_substrate(rng, n_val: int, n_act: int, na: int, nb: int):
    """A synthetic CANONICAL mean field: the Dyall Fock comes out
    diagonal by construction (``h1 = chosen diag - the gamma field``),
    so the canonical energies are the chosen diagonal and every
    Dyall ingredient has a closed-form surface.  Returns the tensors,
    the single-determinant reference CI, its occupied spin-orbitals
    (VALENCE labels), and the chosen diagonal."""
    h, g = _physical_integrals(rng, n_val)
    dets_act = determinant_basis(n_act, na, nb)
    ci = np.zeros(len(dets_act))
    ci[0] = 1.0  # the first CAS determinant as the reference
    # the CAS det in VALENCE labels: the active-first relabel maps
    # the beta block through n_val (alpha q -> q, beta q -> n_val+q)
    occ = tuple(o if o < n_act else n_val + (o - n_act)
                for o in dets_act[0])
    gamma = np.zeros(2 * n_val)
    gamma[list(occ)] = 1.0
    field = (np.einsum("k,pkqk->pq", gamma, g)
             - np.einsum("k,pkkq->pq", gamma, g))
    chosen = np.array([1.1 + 0.37 * i for i in range(2 * n_val)])
    h1 = np.diag(chosen) - field
    return h1, g, ci, occ, chosen


def test_dyall_closed_form_denominators():
    # the closed-form surface of the ticket's MP2 identity, in the
    # reachable in-convention form: a canonical mean field + a
    # single-determinant reference.  The Dyall-specific ingredient --
    # the denominators -- is certified in CLOSED FORM (the chosen
    # diagonal minus the canonical energies, per excited column);
    # the numerators ride the dense coupling column (detbasis's own
    # certified machinery).  The singles are exactly uncoupled HERE
    # (the constructed mean field is globally diagonal -- Brillouin's
    # cancellation); on non-canonical mean fields the [V] singles
    # couple and the meters read them (the campaign cells' class).
    rng = np.random.default_rng(51)
    n_val, n_act, na, nb = 4, 1, 1, 1
    h1, g, ci, occ, chosen = _canonical_substrate(rng, n_val, n_act,
                                                  na, nb)
    h0 = dyall_h0(h1, g, n_act, na, nb, ci)
    res = dyall_e2(h0, n_act, na, nb, ci)

    dets = determinant_basis(n_val, na, nb)
    # the CAS dets in VALENCE labels (the beta block maps through
    # n_val); the map is order-preserving, so ci rides in order
    cas_list = [tuple(o if o < n_act else n_val + (o - n_act)
                      for o in d)
                for d in determinant_basis(n_act, na, nb)]
    cas_set = set(cas_list)
    v = np.zeros(len(dets))
    for cj, dv in enumerate(cas_list):
        v[dets.index(dv)] = ci[cj]
    H = slater_condon_matrix(dets, one_body=h0.h1, two_body=h0.g)
    couplings = H @ v

    act_set = set(h0.act.tolist())
    virt_list = h0.virt.tolist()

    def val(p):
        # the diagonal value of spin-orbital p: the chosen diagonal
        # on the active columns (unrotated), the canonical energy on
        # the virtual ones
        return (chosen[p] if p in act_set
                else h0.eps[virt_list.index(p)])

    tr_mg = sum(chosen[p] for p in occ)  # gamma: point occupations
    total = 0.0
    n_fois = 0
    for mu in range(len(dets)):
        if dets[mu] in cas_set:
            continue
        num = couplings[mu]
        if num == 0.0:
            continue
        n_fois += 1
        gap = tr_mg - sum(val(p) for p in dets[mu])
        total += float(num * np.conj(num)) / gap
    assert res.n_fois == n_fois > 0
    assert res.e2 == pytest.approx(total, abs=1e-10)

    # the doubles' denominators are the full orbital-gap sums (the
    # MP2 form): both occupied spin-orbitals replaced
    for pair in combinations(virt_list, 2):
        a, b = pair
        gap = tr_mg - val(a) - val(b)
        assert gap == pytest.approx(chosen[occ[0]] + chosen[occ[1]]
                                    - val(a) - val(b))


def test_dyall_e0_invariant_under_canonicalization():
    # the virtual rotation is span-preserving: the CAS block (hence
    # E0) survives bit-exactly -- the rotation's correctness pin
    rng = np.random.default_rng(52)
    n_val, n_act, na, nb = 4, 2, 2, 0
    h1, g = _physical_integrals(rng, n_val)
    dets_act = determinant_basis(n_act, na, nb)
    ci = rng.standard_normal(len(dets_act))
    ci /= np.linalg.norm(ci)
    h0 = dyall_h0(h1, g, n_act, na, nb, ci)
    assert epstein_nesbet(h1, g, n_act, na, nb, ci).e0 \
        == pytest.approx(dyall_e2(h0, n_act, na, nb, ci).e0, abs=1e-11)


def test_dyall_matches_dense_twin():
    # the two-construction certificate on hermitian random integrals
    # (the Fock build rides eigh -- hermitian by construction; the
    # complex non-symmetric class stays with the EN twins)
    for seed, n_val, n_act, na, nb in ((61, 4, 2, 2, 0),
                                       (62, 4, 2, 1, 1),
                                       (63, 5, 3, 2, 1)):
        rng = np.random.default_rng(seed)
        h1, g = _physical_integrals(rng, n_val)
        dets_act = determinant_basis(n_act, na, nb)
        ci = rng.standard_normal(len(dets_act))
        ci /= np.linalg.norm(ci)
        h0 = dyall_h0(h1, g, n_act, na, nb, ci)
        a = dyall_e2(h0, n_act, na, nb, ci)
        b = dyall_e2_dense(h0, n_act, na, nb, ci)
        assert a.n_fois == b.n_fois > 0, (seed,)
        for f in ("e0", "e2", "floor", "w_ref"):
            assert np.allclose(getattr(a, f), getattr(b, f)), (f, seed)


def test_dyall_whole_space_vanishes():
    # CAS = the valence space: no external range, no FOIS, E(2) = 0
    rng = np.random.default_rng(71)
    n = 3
    h1, g = _physical_integrals(rng, n)
    dets = determinant_basis(n, 1, 1)
    ci = rng.standard_normal(len(dets))
    ci /= np.linalg.norm(ci)
    h0 = dyall_h0(h1, g, n, 1, 1, ci)
    res = dyall_e2(h0, n, 1, 1, ci)
    assert res.e2 == 0.0
    assert res.n_fois == 0
    assert res.floor == float("inf")
    assert res.w_ref == 1.0


def test_dyall_h0_block_structure():
    # the H0 matrix is block-diagonal by construction (zero cross
    # block); the virtual block of the mean-field Fock rebuilt from
    # the ROTATED tensors comes out diagonal with the canonical
    # energies -- the partition's own consistency twin
    rng = np.random.default_rng(72)
    n_val, n_act, na, nb = 4, 2, 1, 1
    h1, g = _physical_integrals(rng, n_val)
    dets_act = determinant_basis(n_act, na, nb)
    ci = rng.standard_normal(len(dets_act))
    ci /= np.linalg.norm(ci)
    h0 = dyall_h0(h1, g, n_act, na, nb, ci)
    act, virt = h0.act, h0.virt
    assert np.allclose(h0.m[np.ix_(act, virt)], 0.0)
    assert np.allclose(h0.m[np.ix_(virt, virt)],
                       np.diag(h0.eps), atol=1e-12)
    gdiag = np.real(np.diag(h0.gamma))
    f = h0.h1 + np.einsum("k,pkqk->pq", gdiag, h0.g) \
        - np.einsum("k,pkkq->pq", gdiag, h0.g)
    assert np.allclose(f[np.ix_(virt, virt)],
                       np.diag(h0.eps), atol=1e-10)


def test_dyall_refuses_non_physical_integrals():
    # the Fock build's contract: PHYSICAL integrals only -- the
    # "hermitian" kind carries the pair swap but not the within-pair
    # symmetries, and the built Fock would not be hermitian (eigh
    # would silently read one triangle).  Fail fast at the seam.
    rng = np.random.default_rng(73)
    n_val, n_act, na, nb = 4, 2, 1, 1
    h1, g = random_integrals(rng, n_val, "hermitian")
    assert not np.allclose(g, g.transpose(1, 0, 3, 2))
    dets_act = determinant_basis(n_act, na, nb)
    ci = rng.standard_normal(len(dets_act))
    ci /= np.linalg.norm(ci)
    with pytest.raises(ValueError, match="chemist 8-fold"):
        dyall_h0(h1, g, n_act, na, nb, ci)
