"""The nested-basis seam: the contracted basis IS a declared
realization of its own primitive split, so every window construction
restricted to the contraction reproduces the contracted calculation
-- the shell-layer algebra certified against libcint's own
contraction (the carbon fixture at unit sizes: no SCF, integrals
only).  The restriction columns carry the contraction
normalization: libcint normalizes each contracted function, so the
published coefficients need the per-column factor measured in the
split metric."""

import numpy as np
import pytest
from pyscf import gto

from atomterm.detbasis import determinant_basis, slater_condon_matrix
from atomterm.operators import (l2_operator_parts, spin_block_one_body,
                                spin_block_two_body)
from atomterm.system import AtomConfig, build_mol
from atomterm.windows import lowdin

pytestmark = pytest.mark.fast

_UNC_BASIS = gto.parse(
    (__import__("pathlib").Path(__file__).resolve().parents[2]
     / "data" / "c_ccpvdz_unc.txt").read_text())


def _mols():
    """The contracted fixture mol and the split-primitive mol (same
    atom, geometry, spin; the split preserves block and primitive
    order)."""
    mol_c = build_mol(AtomConfig())
    mol_u = gto.M(atom="C 0 0 0", basis=_UNC_BASIS, spin=2, verbose=0)
    return mol_c, mol_u


def _restriction(mol_c, mol_u):
    """The contracted AOs as explicit columns of the split AO space:
    contracted AO (block, contraction j, m) = the published
    coefficients on the split primitives' rows, each column
    normalized in the split metric (libcint normalizes the
    contraction; the published coefficients alone do not)."""
    bs_c = mol_c._basis["C"]
    n_u, n_c = mol_u.nao_nr(), mol_c.nao_nr()
    cols = np.zeros((n_u, n_c))
    s_u = mol_u.intor("int1e_ovlp")
    cum_u = k = 0
    for block in bs_c:
        l = block[0]
        prims = block[1:]          # [exp, c1..cK] rows
        nprim = len(prims)
        nct = len(prims[0]) - 1
        dim = 2 * l + 1
        for j in range(nct):
            for m in range(dim):
                for p in range(nprim):
                    cols[cum_u + p * dim + m, k + j * dim + m] = \
                        prims[p][1 + j]
        cum_u += nprim * dim
        k += nct * dim
    assert cum_u == n_u and k == n_c
    for j in range(n_c):           # the contraction normalization
        cols[:, j] /= np.sqrt(cols[:, j] @ s_u @ cols[:, j])
    return cols


def test_restricted_overlap_is_the_contracted_overlap():
    """The split-metric Gram of the restriction columns IS the
    contracted overlap -- the split basis spans the contracted basis
    exactly, convention-for-convention."""
    mol_c, mol_u = _mols()
    c = _restriction(mol_c, mol_u)
    s_res = c.conj().T @ mol_u.intor("int1e_ovlp") @ c
    assert np.abs(s_res - mol_c.intor("int1e_ovlp")).max() < 1e-12


def test_restricted_one_body_is_the_contracted_one_body():
    """The restricted one-body Hamiltonian IS the contracted one."""
    mol_c, mol_u = _mols()
    c = _restriction(mol_c, mol_u)
    h_u = (mol_u.intor("int1e_kin") + mol_u.intor("int1e_nuc"))
    h_c = (mol_c.intor("int1e_kin") + mol_c.intor("int1e_nuc"))
    assert np.abs(c.conj().T @ h_u @ c - h_c).max() < 1e-11


def test_restricted_two_body_is_the_contracted_two_body():
    """The restricted two-body tensor IS the contracted chemist eri --
    libcint's own contraction reproduced by the explicit column
    contraction (stepwise tensordot, seconds at carbon sizes)."""
    mol_c, mol_u = _mols()
    c = _restriction(mol_c, mol_u)
    eri_u = mol_u.intor("int2e")
    eri_c = mol_c.intor("int2e")
    t = np.tensordot(c, eri_u, axes=(0, 0))   # t[j1, p2, p3, p4]
    t = np.tensordot(c, t, axes=(0, 1))       # t[j1, j2, p3, p4]
    t = np.tensordot(c, t, axes=(0, 2))       # t[j1, j2, j3, p4]
    t = np.tensordot(c, t, axes=(0, 3))       # t[j1, j2, j3, j4]
    assert np.abs(t - eri_c).max() < 1e-10


def test_orthonormalized_restriction_is_the_contracted_orthonormalization():
    """The campaign form: every construction consumes S-orthonormal
    columns, and the S-orthonormalized restriction IS the contracted
    AOs through the same symmetric orthonormalization."""
    mol_c, mol_u = _mols()
    c = _restriction(mol_c, mol_u)
    left = lowdin(c, mol_u.intor("int1e_ovlp"))
    eye_c = np.eye(mol_c.nao_nr())
    right = lowdin(eye_c, mol_c.intor("int1e_ovlp"))
    # the SAME canonical function on isometric spaces: the
    # orthonormalized restriction IS the contracted orthonormalization
    # mapped through the restriction columns
    assert np.abs(left - c @ right).max() < 1e-9


def test_restricted_window_reproduces_the_contracted_term_spectrum():
    """The payload invariant the campaign consumes: a window cut on
    the restriction (the 2p band, the triplet p2 sector) carries the
    SAME Slater-Condon H, the SAME L2 operator, and the SAME term
    spectrum as the contracted mol's own window -- the nested-basis
    anchor for every shell-layer construction (T-B's quantifier)."""
    mol_c, mol_u = _mols()
    c_all = _restriction(mol_c, mol_u)
    # the contracted 2p band: the p block's first contraction (3 AOs)
    p_start = 3          # 2 contracted s + 1 = the first p AO index
    c_win_u = c_all[:, p_start:p_start + 3]
    h_u = (mol_u.intor("int1e_kin") + mol_u.intor("int1e_nuc"))
    h_c = (mol_c.intor("int1e_kin") + mol_c.intor("int1e_nuc"))
    s_u = mol_u.intor("int1e_ovlp")
    s_c = mol_c.intor("int1e_ovlp")
    # canonical window columns on each side: the restriction on the
    # split space, the raw AO rows on the contracted space
    win_c = lowdin(np.eye(mol_c.nao_nr())[:, p_start:p_start + 3], s_c)
    win_u = lowdin(c_win_u, s_u)
    h1_c = win_c.conj().T @ h_c @ win_c
    h1_u = win_u.conj().T @ h_u @ win_u
    eri_c = mol_c.intor("int2e")
    eri_u = mol_u.intor("int2e")

    def phys_eri(win, eri):
        t = np.tensordot(win, eri, axes=(0, 0))
        t = np.tensordot(win, t, axes=(0, 1))
        t = np.tensordot(win, t, axes=(0, 2))
        t = np.tensordot(win, t, axes=(0, 3))
        # the boundary: the project tensor is the physicist
        # <pq|rs> = chemist (pr|qs) -- the ONE transpose
        return t.transpose(0, 2, 1, 3)

    g_c = phys_eri(win_c, eri_c)
    g_u = phys_eri(win_u, eri_u)
    dets = determinant_basis(3, 2, 0)
    h_m_c = (slater_condon_matrix(dets, one_body=spin_block_one_body(h1_c))
             + slater_condon_matrix(dets, two_body=spin_block_two_body(g_c)))
    h_m_u = (slater_condon_matrix(dets, one_body=spin_block_one_body(h1_u))
             + slater_condon_matrix(dets, two_body=spin_block_two_body(g_u)))
    l2_c = slater_condon_matrix(
        dets, *l2_operator_parts(mol_c, win_c))
    l2_u = slater_condon_matrix(
        dets, *l2_operator_parts(mol_u, win_u))
    assert np.abs(h_m_c - h_m_u).max() < 1e-9
    assert np.abs(l2_c - l2_u).max() < 1e-9
    evals_c = np.linalg.eigvalsh(h_m_c)
    evals_u = np.linalg.eigvalsh(h_m_u)
    assert np.abs(evals_c - evals_u).max() < 1e-9
    # the known answer: the triplet p2 ground term is 3P, L2 = 2
    _, evecs = np.linalg.eigh(h_m_c)
    v0 = evecs[:, 0]
    assert abs(float(np.real(v0.conj() @ l2_c @ v0)) - 2.0) < 1e-6
