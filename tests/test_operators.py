"""Tests of many-electron operator matrices for the carbon CAS(2,3) space."""

import numpy as np
import pytest

from atomterm.detbasis import commutator, frobenius_norm
from atomterm.drivers import build_full_mo
from atomterm.operators import build_operators
from atomterm.analysis import analyze, quantum_number, term_symbol


@pytest.fixture(scope="module")
def shell_ops(carbon):
    """Operators in the pure 2p-shell determinant spaces (Ms=1 and Ms=0)."""
    _, mol, mf, act_idx = carbon
    iao = np.eye(mol.nao_nr())
    c_act = iao[:, act_idx]
    ncore = (mol.nelectron - 2) // 2  # 2 for both (2,0) and (1,1)
    ops1 = build_operators(mol, mf, c_act, ncore, na=2, nb=0)
    ops0 = build_operators(mol, mf, c_act, ncore, na=1, nb=1)
    return mol, mf, ops1, ops0


def test_hamiltonian_matches_pyscf_fci(carbon):
    """Slater-Condon rules and phases, validated three independent ways.

    (a) my H spectrum == pyscf's FCI solver on the same integrals (pyscf's
        4-fold input is chemist (pq|rs); ours is physicist <pq|rs>, hence the
        transpose(0, 2, 1, 3));
    (b) my integrals fold the core like mcscf: eig(H) == mcscf's own CI energy
        (to ~1e-4: mcscf uses the actual core density, ours the spherical
        average -- they differ only by SCF numerical noise);
    The ground truth for both is the hand-check E(3P det) = 2 h_2p + J - K.
    """
    from pyscf import fci, mcscf

    _, mol, mf, act_idx = carbon
    iao = np.eye(mol.nao_nr())
    c_act = iao[:, act_idx]
    ncore = (mol.nelectron - 2) // 2
    mo = build_full_mo(mf, list(range(ncore)), act_idx)

    for na, nb in [(2, 0), (1, 1)]:
        ops = build_operators(mol, mf, c_act, ncore, na=na, nb=nb)
        # (a) same operator, independent solver: pyscf's 4-fold input is
        # chemist (pq|rs), ours physicist <pq|rs> = chemist (pr|qs) -> transpose
        e_fci = np.atleast_1d(
            fci.direct_spin1.kernel(
                ops["h1"], ops["g"].transpose(0, 2, 1, 3), 3, (na, nb), conv_tol=1e-13
            )[0]
        )
        e_mine = np.linalg.eigvalsh(ops["H"])[: len(e_fci)]
        assert np.allclose(e_mine, np.sort(e_fci), atol=1e-9), (e_mine, e_fci)

        # (b) folding convention vs mcscf (mo passed explicitly: CASCI methods
        # re-process mc.mo_coeff otherwise); mcscf keeps its chemist convention
        mc = mcscf.CASCI(mf, 3, (na, nb))
        mc.canonicalization = False
        h1eff, _ = mc.get_h1eff(mo)
        g_chem = ao2mo_restore(np.asarray(mc.get_h2eff(mo)), 3)
        e_mcscf = np.atleast_1d(
            fci.direct_spin1.kernel(
                np.asarray(h1eff), g_chem, 3, (na, nb), conv_tol=1e-13
            )[0]
        )
        # 1e-4: mcscf folds the actual core density, ours the spherical
        # average -- they differ only by SCF numerical noise (~1e-5 in h1)
        assert np.allclose(e_mine, np.sort(e_mcscf), atol=1e-4), (e_mine, e_mcscf)


def ao2mo_restore(eri_comp, norb):
    from pyscf import ao2mo

    return ao2mo.restore(1, eri_comp, norb)


def test_3p_threefold_degeneracy_ms1(shell_ops):
    """Ms=1 sector of pure 2p^2: exactly one 3P triplet, threefold M_L."""
    _, _, ops1, _ = shell_ops
    evals = np.linalg.eigvalsh(ops1["H"])
    assert np.allclose(evals - evals[0], 0.0, atol=1e-8), evals
    l2_vals = np.sort(np.real(np.linalg.eigvals(ops1["L2"])))
    assert np.allclose(l2_vals, 2.0, atol=1e-8)
    s2_vals = np.sort(np.real(np.linalg.eigvals(ops1["S2"])))
    assert np.allclose(s2_vals, 2.0, atol=1e-8)


def test_full_multiplet_ms0(shell_ops):
    """Ms=0 sector: 1S (x1, L=0), 1D (x5, L=2), 3P Ms=0 components (x3, L=1)."""
    _, _, _, ops0 = shell_ops
    evals, evecs = np.linalg.eigh(np.real(ops0["H"]))
    l2_exp = np.real(np.diag(evecs.conj().T @ ops0["L2"] @ evecs))
    s2_exp = np.real(np.diag(evecs.conj().T @ ops0["S2"] @ evecs))

    terms = [term_symbol(s2_exp[i], l2_exp[i]) for i in range(9)]
    from collections import Counter

    counts = Counter(terms)
    assert counts == {"3P": 3, "1D": 5, "1S": 1}, counts
    assert np.allclose(np.sort(l2_exp), [0] + [2] * 3 + [6] * 5, atol=1e-7)


def test_commutators_vanish_in_shell_space(shell_ops):
    """[H, Lz] and [H, L^2] vanish identically: the shell space is closed."""
    _, _, ops1, ops0 = shell_ops
    for ops in (ops1, ops0):
        assert frobenius_norm(commutator(ops["H"], ops["Lz"])) < 1e-9
        assert frobenius_norm(commutator(ops["H"], ops["L2"])) < 1e-9


def test_ground_state_is_3p(shell_ops):
    _, _, ops1, _ = shell_ops
    res = analyze(ops1["H"], ops1["Lz"], ops1["L2"], ops1["S2"], "Ms=1")
    assert res.states[0].term == "3P"
    assert res.states[0].l == 1 and res.states[0].s == 1


def test_lz_spectrum_ms0(shell_ops):
    """Lz eigenvalues in the Ms=0 space: 1S(0), 1D(5 x {-2..2}), 3P(3 x {-1..1})."""
    _, _, _, ops0 = shell_ops
    vals = np.sort(np.real(np.linalg.eigvals(ops0["Lz"])))
    expected = np.sort([0.0, -2, -1, 0, 1, 2, -1, 0, 1])
    assert np.allclose(vals, expected, atol=1e-8)


def test_quantum_number_and_term_helpers():
    assert quantum_number(2.0) == 1
    assert quantum_number(6.0) == 2
    assert quantum_number(0.0) == 0
    assert quantum_number(2.0 + 1e-9) == 1
    assert term_symbol(2.0, 2.0) == "3P"
    assert term_symbol(0.0, 6.0) == "1D"
