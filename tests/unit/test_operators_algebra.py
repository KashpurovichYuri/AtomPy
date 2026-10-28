"""Operator algebra without an SCF reference: one-body symbols, the
ladder twin constructions, the parity/population labels -- everything
buildable on a fresh mol (or pure numpy) without an SCF. The SCF-bearing
shell-space and closure-certificate tests live in
tests/integration/test_operators.py."""

import numpy as np
import pytest

from atomterm.detbasis import determinant_basis, frobenius_norm, slater_condon_matrix
from atomterm.windows import select_valence_ao

pytestmark = pytest.mark.fast


def test_parity_operator_diagonal_and_product():
    """Parity: raw AO columns give (-1)^l per orbital; the determinant
    operator is diagonal with eigenvalue = the PRODUCT over occupied
    orbitals (multiplicative -- outside the additive parts protocol). A
    p^2 sector is all-even; an {s,p} sector carries both parities."""
    from atomterm.operators import parity_one_body, parity_operator
    from atomterm.system import AtomConfig, build_mol

    mol = build_mol(AtomConfig(verbose=0))
    p_rows = select_valence_ao(mol, ["2p"])
    c_p = np.eye(mol.nao_nr())[:, p_rows]
    par_p = parity_one_body(mol, c_p)
    assert np.allclose(np.diag(par_p), -1.0)            # l = 1 orbitals

    op = parity_operator(mol, c_p, 2, 0)                # p^2, Ms = 1
    assert np.allclose(op, np.eye(op.shape[0]))         # two p -> even

    sp = select_valence_ao(mol, ["2s", "2p"])
    c_sp = np.eye(mol.nao_nr())[:, sp]
    op_sp = parity_operator(mol, c_sp, 1, 1)            # {2s,2p}, Ms = 0
    assert np.allclose(op_sp, np.diag(np.diag(op_sp)))  # det-diagonal
    assert set(np.round(np.diag(op_sp), 9)) == {-1.0, 1.0}

    mixed = 0.5 * c_sp[:, :1] + 0.5 * c_sp[:, 3:4]  # s/p l-MIX: no
    with pytest.raises(ValueError):                 # definite parity
        parity_one_body(mol, np.hstack([c_sp, mixed]))
    # same-l rotations (Loewdin of overlapping primitives) stay pure
    rot = 0.6 * c_p[:, :1] + 0.8 * c_p[:, 1:2]
    assert np.allclose(np.diag(parity_one_body(mol, rot)), -1.0)


def test_population_operator_partitions_the_window():
    """Population (per-l electron counting for configuration tables): the
    S-metric row projector makes l-disjoint row sets PARTITION the number
    operator elementwise (exact angular orthogonality on one center), a
    raw-AO determinant's diagonal entry counts its occupied spin-orbitals
    of that l exactly, and the one-particle spectrum is {0, 1}."""
    from atomterm.operators import population_operator
    from atomterm.system import AtomConfig, build_mol

    mol = build_mol(AtomConfig(verbose=0))
    iao = np.eye(mol.nao_nr())
    s_rows = select_valence_ao(mol, ["2s"])
    p_rows = select_valence_ao(mol, ["2p"])
    sp_rows = select_valence_ao(mol, ["2s", "2p"])
    c_sp = iao[:, sp_rows]

    na, nb = 1, 1
    n_s = population_operator(mol, c_sp, s_rows, na, nb)
    n_p = population_operator(mol, c_sp, p_rows, na, nb)
    n_all = population_operator(mol, c_sp, sp_rows, na, nb)

    # the partition: different-l row sets are exactly S-orthogonal, so
    # the two projectors resolve the identity on the window span
    assert frobenius_norm(n_s + n_p - n_all) < 1e-12

    # raw AO columns: each determinant's diagonal IS its occupation count
    dets = determinant_basis(len(sp_rows), na, nb)
    expected = np.array([
        sum(1 for so in det if so % len(sp_rows) == 0)  # 2s is column 0
        for det in dets
    ], dtype=float)
    assert np.allclose(np.diag(n_s), expected)

    # many-body spectrum: the possible occupation counts, integers only
    eigs = np.linalg.eigvalsh(n_p)
    assert set(np.round(eigs, 9)) <= {0.0, 1.0, 2.0}


def test_lz2_operator_product_identity_and_loewner(shell_ops):
    """Lz2 == Lz@Lz elementwise, and L^2 - Lz^2 is psd on every sector
    (the Casimir bound <Lz^2> <= <L^2> at the operator level)."""
    _, _, ops1, ops0 = shell_ops
    for ops in (ops1, ops0):
        assert np.allclose(ops["Lz2"], ops["Lz"] @ ops["Lz"], atol=1e-12)
        eig = np.linalg.eigvalsh(np.real(ops["L2"] - ops["Lz2"]))
        assert eig.min() > -1e-10


def test_sz_matrix_is_the_exact_sector_label():
    """Sz (+1/2 alpha, -1/2 beta): its many-body matrix is exactly
    M_S * identity inside every (na, nb) sector -- the sector label.
    Nothing further to certify: H is spin-free, so [H, Sz] vanishes
    identically (that is why the certificates are Lz/L^2 only)."""
    from atomterm.operators import sz_operator

    for norb, na, nb in [(3, 2, 0), (3, 1, 1), (4, 3, 1), (2, 0, 2)]:
        dets = determinant_basis(norb, na, nb)
        sz = slater_condon_matrix(dets, one_body=sz_operator(norb))
        assert np.allclose(sz, ((na - nb) / 2) * np.eye(len(dets))), (na, nb)


def test_spin_ladders_one_body_algebra():
    """The one-electron ladder algebra the S^2 pair terms rest on:
    {s+, s-} = 1, [sz, s±] = ±s±, s±² = 0 (spin-1/2)."""
    from atomterm.operators import spin_ladders, sz_operator

    sp, sm = spin_ladders(3)
    sz = sz_operator(3)
    eye = np.eye(2 * 3)
    assert np.allclose(sp @ sm + sm @ sp, eye)
    assert np.allclose(sp @ sp, 0) and np.allclose(sm @ sm, 0)
    assert np.allclose(sz @ sp - sp @ sz, sp)
    assert np.allclose(sz @ sm - sm @ sz, -sm)


def test_orbital_ladders_twin_construct_l2():
    """L± = Lx ± iLy as one-body window operators. The independent
    ladder construction L^2 = Lz^2 + (L+ L- + L- L+)/2 is elementwise
    identical to the delivered L^2 parts, and [Lz, L±] = ±L± pins the
    m-ladder structure (p², Ms=0 sector)."""
    from atomterm.operators import (
        l2_operator_parts,
        l_ladder_operator,
        lz_operator,
    )
    from atomterm.system import AtomConfig, build_mol

    mol = build_mol(AtomConfig(verbose=0))
    c = np.eye(mol.nao_nr())[:, select_valence_ao(mol, ["2p"])]
    dets = determinant_basis(3, 1, 1)

    ladders = l_ladder_operator(mol, c)
    lz = slater_condon_matrix(dets, one_body=lz_operator(mol, c))
    lp = slater_condon_matrix(dets, one_body=ladders["+"])
    lm = slater_condon_matrix(dets, one_body=ladders["-"])
    assert np.allclose(lz @ lp - lp @ lz, lp)
    assert np.allclose(lz @ lm - lm @ lz, -lm)

    one, two = l2_operator_parts(mol, c)
    l2 = slater_condon_matrix(dets, one_body=one, two_body=two)
    twin = lz @ lz + 0.5 * (lp @ lm + lm @ lp)
    assert np.allclose(l2, twin)


def test_spin_ladders_twin_construct_s2_across_sectors():
    """S± as explicit inter-sector maps (the one operator class that
    does NOT preserve a sector): S² = Sz² + (S- S+ + S+ S-)/2 holds
    elementwise on EVERY sector of the p² multiplet -- the pair
    products leave and return home. Strongest check of the Ŝ²
    flip-flop terms and of the inter-sector signs."""
    from atomterm.operators import (
        s2_operator_parts,
        spin_ladder_matrix,
        sz_operator,
    )

    norb = 3
    one, two = s2_operator_parts(norb)
    sz1 = sz_operator(norb)

    def sector(dets):
        s2 = slater_condon_matrix(dets, one_body=one, two_body=two)
        sz = slater_condon_matrix(dets, one_body=sz1)
        return s2, sz

    dets_0 = determinant_basis(norb, 1, 1)
    dets_up = determinant_basis(norb, 2, 0)
    dets_dn = determinant_basis(norb, 0, 2)
    s2_0, sz_0 = sector(dets_0)
    s2_up, sz_up = sector(dets_up)
    s2_dn, sz_dn = sector(dets_dn)

    sp_0u = spin_ladder_matrix(norb, dets_0, dets_up, plus=True)   # (1,1)->(2,0)
    sm_0d = spin_ladder_matrix(norb, dets_0, dets_dn, plus=False)  # (1,1)->(0,2)

    # Ms=0: both pair products leave and come home
    twin_0 = sz_0 @ sz_0 + 0.5 * (
        sp_0u.conj().T @ sp_0u + sm_0d.conj().T @ sm_0d
    )
    assert np.allclose(s2_0, twin_0)
    # the flipped sectors see only the product through Ms=0
    twin_up = sz_up @ sz_up + 0.5 * (sp_0u @ sp_0u.conj().T)
    assert np.allclose(s2_up, twin_up)
    twin_dn = sz_dn @ sz_dn + 0.5 * (sm_0d @ sm_0d.conj().T)
    assert np.allclose(s2_dn, twin_dn)
