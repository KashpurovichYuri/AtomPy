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


def test_parity_one_body_exact_on_pure_windows():
    """(-1)^l per column: the p window and an f column read -1, the
    {s,p} window diag(+1, -1, -1, -1), off-diagonal exactly zero;
    same-l rotations stay pure."""
    from pyscf import gto

    from atomterm.operators import parity_one_body
    from atomterm.system import AtomConfig, build_mol

    mol = build_mol(AtomConfig(verbose=0))
    c_p = np.eye(mol.nao_nr())[:, select_valence_ao(mol, ["2p"])]
    assert np.allclose(parity_one_body(mol, c_p), -np.eye(3), atol=1e-12)

    c_sp = np.eye(mol.nao_nr())[:, select_valence_ao(mol, ["2s", "2p"])]
    assert np.allclose(parity_one_body(mol, c_sp),
                       np.diag([1.0, -1.0, -1.0, -1.0]), atol=1e-12)

    # same-l rotations (Loewdin of overlapping primitives) stay pure
    rot = 0.6 * c_p[:, :1] + 0.8 * c_p[:, 1:2]
    assert np.allclose(np.diag(parity_one_body(mol, rot)), -1.0)

    # an f shell: odd l reads -1 (the shell walker, no label parsing)
    mol_f = gto.M(atom="C 0 0 0",
                  basis=[[0, [3.0, 1.0]], [3, [1.0, 1.0]]], verbose=0)
    c_f = np.eye(mol_f.nao_nr())[:, 1:2]
    assert np.allclose(parity_one_body(mol_f, c_f), -np.eye(1), atol=1e-12)


def test_parity_one_body_mixed_columns_are_data():
    """An l-mixed column is DATA, not an error: the window
    [2s, (2s + 2p0)/sqrt(2)] gives the hand-derived
    Pi = [[1, 1/sqrt(2)], [1/sqrt(2), 0]] with eigenvalues
    (1 +- sqrt(3))/2."""
    from atomterm.operators import parity_one_body
    from atomterm.system import AtomConfig, build_mol

    mol = build_mol(AtomConfig(verbose=0))
    sp = select_valence_ao(mol, ["2s", "2p"])
    eye = np.eye(mol.nao_nr())
    w = np.hstack([eye[:, sp[:1]],
                   (eye[:, sp[:1]] + eye[:, sp[3:4]]) / np.sqrt(2.0)])

    pi = parity_one_body(mol, w)
    expect = np.array([[1.0, 1 / np.sqrt(2)], [1 / np.sqrt(2), 0.0]])
    assert np.allclose(pi, expect, atol=1e-12)
    vals = np.sort(np.linalg.eigvalsh(pi))
    assert np.allclose(vals, [(1 - np.sqrt(3)) / 2, (1 + np.sqrt(3)) / 2],
                       atol=1e-12)


def test_parity_operator_exterior_power_exact():
    """Lambda^N(Pi): on l-pure windows diagonal with the PRODUCT of
    occupied signs (the old construction, elementwise) and an exact
    involution -- parity-closed is where +-1 is sharp. On the mixed
    window [2s, (2s+2p0)/sqrt(2)] the full 4x4 sector matrix is
    hand-derived from the minors: first row [1, r, r, 1/2],
    r = 1/sqrt(2) -- the honest dense object."""
    from atomterm.operators import parity_operator
    from atomterm.system import AtomConfig, build_mol

    mol = build_mol(AtomConfig(verbose=0))
    eye = np.eye(mol.nao_nr())
    c_p = eye[:, select_valence_ao(mol, ["2p"])]

    op = parity_operator(mol, c_p, 2, 0)                # p^2: all even
    assert np.allclose(op, np.eye(3), atol=1e-12)
    assert np.allclose(op @ op, np.eye(3), atol=1e-12)  # closed: sharp

    sp = select_valence_ao(mol, ["2s", "2p"])
    c_sp = eye[:, sp]
    op_sp = parity_operator(mol, c_sp, 1, 1)            # {2s,2p}, Ms = 0
    signs = np.array([1.0, -1.0, -1.0, -1.0])           # the product form
    dets = determinant_basis(4, 1, 1)
    product = np.diag([float(np.prod([signs[d % 4] for d in det]))
                       for det in dets])
    assert np.allclose(op_sp, product, atol=1e-12)
    assert np.allclose(op_sp @ op_sp, np.eye(len(op_sp)), atol=1e-12)

    mixed_w = np.hstack([eye[:, sp[:1]],
                         (eye[:, sp[:1]] + eye[:, sp[3:4]]) / np.sqrt(2.0)])
    op_m = parity_operator(mol, mixed_w, 1, 1)
    r = 1 / np.sqrt(2.0)
    expect = np.array([
        [1.0, r, r, 0.5],
        [r, 0.0, 0.5, 0.0],
        [r, 0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0, 0.0],
    ])
    assert np.allclose(op_m, expect, atol=1e-12)        # minors, by hand


def test_one_electron_sector_is_exact_one_body():
    """One electron -- the degenerate anchor: the many-body operators
    ARE the one-body ones. p^1: L^2 = 2 I, Lz spectrum {-1, 0, 1},
    S^2 = 3/4 I, N = I; s^1: L^2 = 0 identically."""
    from atomterm.operators import (
        l2_operator_parts,
        lz_operator,
        number_operator,
        s2_operator_parts,
    )
    from atomterm.system import AtomConfig, build_mol

    mol = build_mol(AtomConfig(verbose=0))
    eye = np.eye(mol.nao_nr())

    for cols, l2_expect, lz_vals in [
        (select_valence_ao(mol, ["2p"]), 2.0, [-1, 0, 1]),
        (select_valence_ao(mol, ["2s"]), 0.0, [0]),
    ]:
        c = eye[:, cols]
        dets = determinant_basis(c.shape[1], 1, 0)
        one, two = l2_operator_parts(mol, c)
        l2 = slater_condon_matrix(dets, one_body=one, two_body=two)
        s2_one, s2_two = s2_operator_parts(c.shape[1])
        s2 = slater_condon_matrix(dets, one_body=s2_one, two_body=s2_two)
        lz = slater_condon_matrix(dets, one_body=lz_operator(mol, c))
        n = slater_condon_matrix(dets, one_body=number_operator(c.shape[1]))
        eye_n = np.eye(len(dets))
        assert np.allclose(l2, l2_expect * eye_n, atol=1e-12)
        assert np.allclose(s2, 0.75 * eye_n, atol=1e-12)
        assert np.allclose(n, eye_n, atol=1e-12)
        assert np.allclose(np.sort(np.real(np.linalg.eigvals(lz))),
                           lz_vals, atol=1e-12)


def test_closure_of_the_whole_ao_space_is_exact(p_shell_mol):
    """P = I on the full AO space: both closure defects exactly zero
    -- the degenerate anchor of the meter."""
    from atomterm.operators import window_closure

    c = window_closure(p_shell_mol, np.eye(p_shell_mol.nao_nr()))
    assert c["lz"] < 1e-12 and c["lperp"] < 1e-12


def test_build_operators_method_seam_elementwise_identical(carbon):
    """The public dispatch: both evaluators over the same parts give
    elementwise-identical matrices (the evaluators are pinned on
    random integrals; this pins the seam that feeds them)."""
    from atomterm.operators import build_operators

    _, mol, mf, act_idx = carbon
    c_act = np.eye(mol.nao_nr())[:, act_idx]
    a = build_operators(mol, mf, c_act, 2, 1, 1)
    b = build_operators(mol, mf, c_act, 2, 1, 1, method="secondquant")
    for k in ("H", "Lz", "Lz2", "L2", "S2", "Sz", "N"):
        assert np.allclose(a[k], b[k], atol=1e-12), k


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


def test_pyscf_apply_is_the_slater_condon_action(carbon):
    """The C-contraction boundary: pyscf_apply over the kernel feed's
    own (h1e, eri) reproduces our Slater-Condon Hamiltonian's action
    on random vectors (the vector-level twin of the lane's root
    anchor), real and complex; a kernel root is an eigenvector of
    ours; the state residual [H, -i Lz] applied through the boundary
    + the numpy apply vanishes on the term root."""
    from pyscf import ao2mo
    from pyscf.fci import direct_spin1

    from atomterm.detbasis import operator_apply
    from atomterm.operators import (active_space_integrals,
                                    build_operator_parts,
                                    ci_in_detbasis_order, pyscf_apply)
    from atomterm.windows import build_full_mo, occupied_core_idx

    _, mol, mf, _ = carbon
    ncore = 2
    mo = build_full_mo(mf, occupied_core_idx(mf, ncore),
                       mf.mo_coeff[:, ncore:])
    c = mo[:, ncore:]
    na, nb = 2, 1
    norb = c.shape[1]
    h1e, _ = active_space_integrals(mol, mf, c, ncore)
    eri = ao2mo.kernel(mol, c)
    parts = build_operator_parts(mol, mf, c, ncore, na, nb)
    dets = parts["dets"]
    ours = operator_apply(dets, one_body=parts["H"][0],
                          two_body=parts["H"][1])

    f = pyscf_apply(h1e, eri, norb, na, nb)
    rng = np.random.default_rng(19)
    for _ in range(2):
        v = rng.standard_normal(len(dets)) + 1j * rng.standard_normal(len(dets))
        assert np.allclose(f(v), ours(v), atol=1e-10)

    # the anchor's vector form: a kernel root is an eigenvector of ours
    e, ci = direct_spin1.kernel(np.real(h1e).astype(float),
                                np.real(eri).astype(float), norb,
                                (na, nb), nroots=1)
    w = ci_in_detbasis_order(np.asarray(ci).ravel(), norb, na, nb)
    assert abs(w @ ours(w) - float(np.atleast_1d(e)[0])) < 1e-10

    # the state residual through the pair of applies: M = -i Lz real,
    # [H, Lz] = i[H, M] -- the norm is M's; on the closed 2p window
    # the root is an Lz eigenstate-ish (the term's |M_L| member) so
    # the residual is small but nonzero at the window level -- the
    # INVARIANT is agreement with the dense commutator's action
    m_apply = operator_apply(dets, one_body=-1j * parts["Lz"][0])
    lz_m = slater_condon_matrix(dets, one_body=-1j * parts["Lz"][0])
    h_m = slater_condon_matrix(dets, one_body=parts["H"][0],
                               two_body=parts["H"][1])
    r_apply = f(m_apply(w)) - m_apply(f(w))
    r_dense = h_m @ (lz_m @ w) - lz_m @ (h_m @ w)
    assert np.allclose(r_apply, r_dense, atol=1e-9)
