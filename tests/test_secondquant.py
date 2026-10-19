"""Second-quantized operator construction, validated against Carbon.

Three validation tiers:
1. elementwise identity with the first-quantized Slater-Condon matrices;
2. core-explicit full-space Hamiltonian (CAS(6,5), no folding) vs PySCF
   CASCI with ncore=0 -- identical integrals on both sides;
3. the fixed-core sector of the full-space H equals the folded active-space
   result plus the core energy (validates the folding conventions).
"""

import numpy as np
import pytest

from atomterm.detbasis import (
    commutator,
    determinant_basis,
    frobenius_norm,
    slater_condon_matrix,
)
pytestmark = pytest.mark.fast
from atomterm.secondquant import annihilate, create, secondquant_matrix
from atomterm.operators import (
    build_operators,
    full_space_integrals,
    lz_operator,
    l2_operator_parts,
    number_operator,
    s2_operator_parts,
    spin_block_one_body,
    spin_block_two_body,
)
from atomterm.drivers import build_full_mo, complete_orthonormal


# ---------------------------------------------------------------- sign rules

def test_creation_annihilation_basics():
    assert create((), 0) == (1, (0,))
    assert create((), 1) == (1, (1,))
    assert create((0,), 0) is None            # Pauli
    assert annihilate((), 0) is None          # nothing to kill
    assert annihilate((0,), 0) == (1, ())
    # parity: a_1 removes the higher orbital -> sign -1
    assert annihilate((0, 1), 1) == (-1, (0,))
    assert annihilate((0, 1), 0) == (1, (1,))


def test_canonical_anticommutation():
    """{a_p, a+_p} = 1 and {a_p, a+_q} = 0 on low Fock states."""
    for occ in [(), (0,), (1,), (0, 1), (0, 2)]:
        for p in range(3):
            total: dict[tuple, complex] = {}
            for first_op in ["create", "annihilate"]:
                if first_op == "create":
                    first = create(occ, p)
                    if first is None:
                        continue
                    sg, mid = first
                    second = annihilate(mid, p)
                else:
                    first = annihilate(occ, p)
                    if first is None:
                        continue
                    sg, mid = first
                    second = create(mid, p)
                if second is None:
                    continue
                sg2, o2 = second
                total[o2] = total.get(o2, 0j) + sg * sg2
            assert total == {occ: 1}


def test_annihilators_anticommute():
    """a_1 a_0 |01> = +|0> but a_0 a_1 |01> = -|0>."""
    occ = (0, 1)
    s1, o1 = annihilate(occ, 0)
    s2, o2 = annihilate(o1, 1)
    assert (s1 * s2, o2) == (1, ())
    s1, o1 = annihilate(occ, 1)
    s2, o2 = annihilate(o1, 0)
    assert (s1 * s2, o2) == (-1, ())


# ------------------------------------------------- tier 1: Carbon, folded core

@pytest.mark.parametrize("na, nb", [(2, 0), (1, 1)])
def test_elementwise_identical_to_slater_condon(carbon, na, nb):
    """Both constructions consume identical parts and produce identical
    matrices -- elementwise, for every operator."""
    _, mol, mf, act_idx = carbon
    c_act = np.eye(mol.nao_nr())[:, act_idx]
    m1 = build_operators(mol, mf, c_act, 2, na, nb, method="slater_condon")
    m2 = build_operators(mol, mf, c_act, 2, na, nb, method="secondquant")
    for key in ["H", "Lz", "L2", "S2", "N"]:
        assert np.allclose(m1[key], m2[key], atol=1e-12), key


@pytest.mark.parametrize("na, nb", [(2, 0), (1, 1)])
def test_textbook_spectra_through_secondquant(carbon, na, nb):
    """The 2p^2 spectra must also come out of the a+/a path."""
    _, mol, mf, act_idx = carbon
    c_act = np.eye(mol.nao_nr())[:, act_idx]
    ops = build_operators(mol, mf, c_act, 2, na, nb, method="secondquant")
    l2 = np.sort(np.real(np.linalg.eigvals(ops["L2"])))
    s2 = np.sort(np.real(np.linalg.eigvals(ops["S2"])))
    if (na, nb) == (2, 0):
        assert np.allclose(l2, 2.0) and np.allclose(s2, 2.0)
    else:
        assert np.allclose(np.sort(l2), [0.0] * 1 + [2.0] * 3 + [6.0] * 5)
        assert np.allclose(np.sort(s2), [0.0] * 6 + [2.0] * 3)


# ------------------------------- tier 2: core-explicit full space, CAS(6,5)

@pytest.fixture(scope="module")
def full_space(carbon):
    """(h1, g, w) for the all-electron CAS(6,5) space.

    Core columns are the canonical core MOs *purified* to exact l = 0
    (p/d tail dropped, pair Loewdin-orthonormalized): close to the SCF
    solution yet exactly orthonormal against the 2p shell AOs and exactly
    rotation-closed. (Raw s AOs work equally well; see wayfinder ticket 05.)
    """
    from atomterm.drivers import purified_s_core

    _, mol, mf, act_idx = carbon
    core, contamination = purified_s_core(mol, mf, 2)
    assert np.all(contamination < 1e-3)  # SCF noise, not physical mixing
    w = np.hstack([core, np.eye(mol.nao_nr())[:, act_idx]])
    assert np.allclose(w.conj().T @ mf.get_ovlp() @ w, np.eye(5), atol=1e-12)
    h1, g = full_space_integrals(mol, mf, w)
    return mol, mf, w, h1, g


@pytest.mark.parametrize("na, nb", [(4, 2), (3, 3)])
def test_full_space_vs_pyscf_casci_core_explicit(full_space, carbon, na, nb):
    """Core-explicit second-quantized H (raw integrals, no folding) vs PySCF
    CASCI with ncore=0: identical integrals -> solver-precision agreement."""
    from pyscf import mcscf

    mol, mf, w, h1, g = full_space
    _, _, _, act_idx = carbon
    enuc = mol.energy_nuc()

    dets = determinant_basis(5, na, nb)
    s2o, s2t = s2_operator_parts(5)
    l2o, l2t = l2_operator_parts(mol, w)
    H = secondquant_matrix(
        dets, one_body=spin_block_one_body(h1), two_body=spin_block_two_body(g)
    )
    Lz = secondquant_matrix(dets, one_body=lz_operator(mol, w))
    L2 = secondquant_matrix(dets, one_body=l2o, two_body=l2t)
    S2 = secondquant_matrix(dets, one_body=s2o, two_body=s2t)
    N = secondquant_matrix(dets, one_body=number_operator(5))
    assert np.allclose(np.diag(N), 6)

    evals, evecs = np.linalg.eigh(H)

    # commutators over the whole 50/100-determinant space
    assert frobenius_norm(commutator(H, Lz)) < 1e-9
    assert frobenius_norm(commutator(H, L2)) < 1e-9

    # ground state: closed s core (L = S = 0) x 2p^2 3P
    l2_0 = float(np.real(evecs[:, 0].conj() @ L2 @ evecs[:, 0]))
    s2_0 = float(np.real(evecs[:, 0].conj() @ S2 @ evecs[:, 0]))
    assert abs(l2_0 - 2.0) < 1e-7 and abs(s2_0 - 2.0) < 1e-7

    # PySCF reference: CAS(6,5) with ncore=0, same orbital columns
    s = mf.get_ovlp()
    mo_full = np.hstack([w, complete_orthonormal(w, s)])
    mc = mcscf.CASCI(mf, 5, (na, nb))
    mc.canonicalization = False
    mc.fcisolver.nroots = 5
    mc.kernel(mo_full)
    e_pyscf = np.sort(np.atleast_1d(np.asarray(mc.e_cas)))[: len(evals)]
    # both are electronic energies (nuclear repulsion excluded on both sides)
    assert np.allclose(evals[: len(e_pyscf)], e_pyscf, atol=1e-8), (
        evals[: len(e_pyscf)],
        e_pyscf,
    )


# ------------------------------------------------------- tier 3: folded core

def test_purified_core_is_pure(carbon):
    """The purified core orbitals are exactly l = 0 (L^2 kills them),
    orthonormal, and numerically almost identical to the canonical MOs.
    The dropped contamination quantifies the SCF symmetry breaking.
    """
    from atomterm.angular import l2_matrix
    from atomterm.drivers import purified_s_core

    _, mol, mf, _ = carbon
    core, contamination = purified_s_core(mol, mf, 2)

    # exact l = 0: L^2 annihilates pure s functions
    assert np.abs(l2_matrix(mol) @ core).max() < 1e-12
    # orthonormal in the AO metric
    s = mf.get_ovlp()
    assert np.allclose(core.conj().T @ s @ core, np.eye(2), atol=1e-12)
    # numerically close to the canonical core (overlap ~ 1)
    overlap = np.diag(mf.mo_coeff[:, :2].conj().T @ s @ core)
    assert np.allclose(np.abs(overlap), 1.0, atol=1e-4)
    # dropped tail: nonzero (SCF noise is real) but tiny
    assert np.all(contamination > 0) and np.all(contamination < 1e-3)
    print(f"\ncore MO contamination (1 - s weight): {contamination}")

def test_fixed_core_sector_matches_folded(full_space, carbon):
    """The doubly-occupied-core block of the full-space H equals the folded
    active-space Hamiltonian shifted by the core energy. The core orbitals
    are pure s functions, so both the explicit and the folded (actual
    density) treatments are exactly spherical -> machine precision.
    """
    from pyscf import mcscf

    _, mol, mf, act_idx = carbon
    _, _, w, h1, g = full_space
    dets = determinant_basis(5, 4, 2)
    H = secondquant_matrix(
        dets, one_body=spin_block_one_body(h1), two_body=spin_block_two_body(g)
    )

    core_sos = (0, 1, 5, 6)  # 1s alpha/beta, 2s alpha/beta (spin-orbitals)
    idx = [i for i, d in enumerate(dets) if all(c in d for c in core_sos)]
    assert len(idx) == 3  # the 2p^2 Ms=1 determinants
    e_block = np.linalg.eigvalsh(H[np.ix_(idx, idx)])

    # explicit core energy: <core-only determinant | H | core-only determinant>
    core_det = secondquant_matrix(
        [core_sos], one_body=spin_block_one_body(h1), two_body=spin_block_two_body(g)
    )[0, 0].real

    # folded reference 1: first quantization, same core columns
    eri = mol.intor("int2e", aosym="s1")
    d_core = 2.0 * w[:, :2] @ w[:, :2].conj().T
    veff = (
        np.einsum("ls,mnls->mn", d_core, eri, optimize=True)
        - 0.5 * np.einsum("ls,mlns->mn", d_core, eri, optimize=True)
    )
    p = act_idx
    h1_fold = (mf.get_hcore() + veff)[np.ix_(p, p)]
    g_fold = eri[np.ix_(p, p, p, p)].transpose(0, 2, 1, 3)
    H_fold = secondquant_matrix(
        determinant_basis(3, 2, 0),
        one_body=spin_block_one_body(h1_fold),
        two_body=spin_block_two_body(g_fold),
    )
    # total energy = folded active-space energy + explicit core energy
    assert np.allclose(e_block, np.linalg.eigvalsh(H_fold) + core_det, atol=1e-10)

    # folded reference 2: pyscf CASCI on the same orbital set
    s = mf.get_ovlp()
    mo5 = np.hstack([w, complete_orthonormal(w, s)])
    mc = mcscf.CASCI(mf, 3, (2, 0))
    mc.canonicalization = False
    mc.kernel(mo5)
    assert np.allclose(np.sort(e_block) + mol.energy_nuc(), mc.e_tot, atol=1e-8)
