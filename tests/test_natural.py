"""Natural-orbital windows (ticket 22): the 1-RDM helper and the m-set
claim.

The closure theorem licenses any radial subspace replicated over ALL
members. These tests pin (a) the CI 1-RDM construction (projector
property, brute-force pair-sum agreement, root averaging), (b) the
structural claim of the ticket: an m-SCALAR density yields NOs in
complete degenerate m-sets (occupation groups of 2l+1), while a SINGLE
|L, M_L> root carries an m-diagonal but m-UNEQUAL density whose top
occupation picks ONE member -- the AVAS mistake pattern -- so root
AVERAGING is what makes NO selection symmetry-safe, and (c) the
certificate verdict on real AO spans (complete vs m-split window).
Pure numpy/libcint: no SCF anywhere (fast tier).
"""

import numpy as np
import pytest

from atomterm.detbasis import determinant_basis
from atomterm.natural import (
    ci_rdm1,
    natural_orbitals,
    occupation_groups,
    select_no_window,
)
from atomterm.windows import (
    Window,
    window_no,
    window_primitive,
    window_scf_mo,
    window_shell,
)
from atomterm.secondquant import annihilate, create
from atomterm.system import AtomConfig, build_mol

pytestmark = pytest.mark.fast


# ------------------------------------------------------------ 1-RDM

def test_rdm1_single_determinant_is_projector():
    dets = determinant_basis(4, 2, 1)
    c = np.zeros(len(dets))
    c[7] = 1.0
    gamma = ci_rdm1(dets, c, 4)
    assert np.allclose(gamma, gamma.T)
    assert np.isclose(np.trace(gamma), 3.0)
    # projector onto the occupied spatial orbitals of that determinant
    det = dets[7]
    expect = np.zeros((4, 4))
    for p in det:
        expect[p % 4, p % 4] += 1.0
    assert np.allclose(gamma @ gamma, gamma)
    assert np.allclose(gamma, expect)


def brute_rdm1(dets, c, norb):
    """Direct det-pair evaluation: gamma[p,q] = sum_ij c_i c_j
    <i| a+_p a_q |j> over both spins (independent of ci_rdm1's loop)."""
    gamma = np.zeros((norb, norb))
    index = {det: i for i, det in enumerate(dets)}
    for j, detj in enumerate(dets):
        for q in range(norb):
            for sq in (0, norb):
                killed = annihilate(detj, q + sq)
                if killed is None:
                    continue
                sg_q, o1 = killed
                for p in range(norb):
                    for sp in (0, norb):
                        made = create(o1, p + sp)
                        if made is None:
                            continue
                        sg_p, target = made
                        i = index.get(target)
                        if i is not None:
                            gamma[p, q] += sg_q * sg_p * c[i] * c[j]
    return gamma


def test_rdm1_matches_brute_force_on_random_ci():
    rng = np.random.default_rng(24)
    dets = determinant_basis(4, 2, 1)
    c = rng.normal(size=len(dets))
    c /= np.linalg.norm(c)
    gamma = ci_rdm1(dets, c, 4)
    assert np.allclose(gamma, brute_rdm1(dets, c, 4))
    assert np.isclose(np.trace(gamma), 3.0)
    assert np.allclose(gamma, gamma.T)


def test_rdm1_root_average_is_linear():
    rng = np.random.default_rng(4)
    dets = determinant_basis(4, 2, 0)
    vec = rng.normal(size=(len(dets), 2))
    g_avg = ci_rdm1(dets, vec, 4, weights=[0.3, 0.7])
    g1 = ci_rdm1(dets, vec[:, 0], 4)
    g2 = ci_rdm1(dets, vec[:, 1], 4)
    assert np.allclose(g_avg, 0.3 * g1 + 0.7 * g2)


# ------------------------------------------------- the m-set claim

def p_scalar_density():
    """Two p radials x three members, m-scalar: gamma = R (x) I_3."""
    return np.kron(np.diag([0.9, 0.5]), np.eye(3))


def test_m_scalar_density_yields_whole_m_sets():
    c_no, occ = natural_orbitals(np.eye(6), p_scalar_density())
    # occupations: two 3-fold degenerate levels, descending
    groups = occupation_groups(occ)
    assert [len(g) for g in groups] == [3, 3]
    assert np.allclose(sorted(occ, reverse=True)[::3], [0.9, 0.5])
    # whole-group selection: sizes 3 and 6, never 1/2/4
    for n_max, want in ((3, 3), (4, 3), (5, 3), (6, 6)):
        c_win, kept = select_no_window(c_no, occ, n_max)
        assert c_win.shape[1] == want


def test_single_ml_root_breaks_the_m_set():
    """The single determinant (p m=0 alpha, p m=0 beta) -- an M_S = 0,
    M_L = 0 state of p^2 (not a term eigenstate; any single det makes
    the point). Its density is m-DIAGONAL but m-UNEQUAL: the top
    occupation selects ONE member -- occupation selection on a single
    root is NOT symmetry-safe (the AVAS mistake, fact 19)."""
    dets = determinant_basis(3, 1, 1)
    det = (0, 3)  # m=0 alpha, m=0 beta
    c = np.zeros(len(dets))
    c[dets.index(det)] = 1.0
    gamma = ci_rdm1(dets, c, 3)
    assert np.allclose(np.diag(gamma), [2.0, 0.0, 0.0])
    c_no, occ = natural_orbitals(np.eye(3), gamma)
    c_win, kept = select_no_window(c_no, occ, 1)
    assert c_win.shape[1] == 1          # one member, not the m-set


def test_no_window_via_registry():
    """The uniform registry entry: toy basis, two roots averaged -- the
    m-scalar route -- gives the whole m-set even at n_max = 1."""
    dets = determinant_basis(3, 1, 1)
    rng = np.random.default_rng(5)
    # two roots whose densities are m-scalar in a shared eigenbasis:
    # root 1 = the (m0a, m0b) determinant, root 2 = (m1a, m1b) -- the
    # average is m-scalar over {m0, m1} but not m2 (honest: the
    # averaged SET must be rotation-invariant; a missing member shows).
    vec = np.zeros((len(dets), 2))
    vec[dets.index((0, 3)), 0] = 1.0
    vec[dets.index((1, 4)), 1] = 1.0
    w = window_no(None, None, None, c_act=np.eye(3), dets=dets,
                  vec=vec, n_max=3, weights=[0.5, 0.5])
    assert w.c_win.shape[1] == 3        # everything fits: all members
    # the degenerate pair is never split, even when only one is asked for
    w1 = window_no(None, None, None, c_act=np.eye(3), dets=dets,
                   vec=vec, n_max=1, weights=[0.5, 0.5])
    assert w1.c_win.shape[1] == 2


# --------------------------------------- certificate on real AO spans

def test_certificate_on_real_ao_windows():
    """Raw 2p AOs of the bare carbon mol: complete m-set closes at
    machine precision; dropping one member is rotation-defective."""
    from atomterm.drivers import select_valence_ao
    from atomterm.operators import window_closure

    mol = build_mol(AtomConfig(verbose=0))
    p_rows = select_valence_ao(mol, ["2p"])
    eye = np.eye(mol.nao_nr())
    c = window_closure(mol, eye[:, p_rows])
    assert c["lz"] < 1e-12 and c["lperp"] < 1e-12
    c = window_closure(mol, eye[:, p_rows[:2]])   # m-split p window
    assert c["lz"] > 1e-2 or c["lperp"] > 1e-2


def test_registry_ao_selectors_uniform_signature():
    mol = build_mol(AtomConfig(verbose=0))
    w = window_shell(mol, None, 2, patterns=["2p"])
    assert w.m_complete and w.ao_rows is not None and w.c_win is None
    assert len(w.ao_rows) == 3
    # the fact-1 guard: exponent-ranked primitive selection is undefined
    # for contracted shells (ccpvdz carbon p is one 3-exponent
    # contraction) -- the primitive selector is for uncontracted bases
    with pytest.raises(ValueError, match="contraction slot"):
        window_primitive(mol, None, 2, bands={1: 1})


def test_registry_scf_mo_selector_slices_mo_columns():
    class FakeMF:
        mo_coeff = np.eye(6)
    w = window_scf_mo(None, FakeMF(), 1, open_cols=[0], rank=[5, 4, 3], k=2)
    assert not w.m_complete and w.ao_rows is None
    assert w.c_win.shape == (6, 3)
    assert np.allclose(w.c_win[:, 1:], np.eye(6)[:, [5, 4]])
