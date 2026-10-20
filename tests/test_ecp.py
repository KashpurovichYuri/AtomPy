"""ECP mode (tickets 18/19): ccECP pairing, valence bookkeeping, exact p^2
multiplets, the ncore=0 full-valence fold, the louder raw-core detector,
the gated TM references, TM free-ion tables, and the fact-20 closure
certificate at ncore=4. All bands carry the measured fresh-draw values in
comments."""

import numpy as np
import pytest

from atomterm import analysis, drivers, operators
from atomterm.drivers import ReferenceSpec
from atomterm.system import AtomConfig, build_mol, ncore_of, sector_occupations

C_ECP = AtomConfig(symbol="C", basis="ccecpccpvdz", ecp="ccecp", spin=2)
SI_ECP = AtomConfig(symbol="Si", basis="ccecpccpvdz", ecp="ccecp", spin=2)
SC_ECP = AtomConfig(symbol="Sc", basis="ccecpccpvdz", ecp="ccecp", spin=1)
TI_ECP = AtomConfig(symbol="Ti", basis="ccecpccpvdz", ecp="ccecp", spin=2)
SC_REF = ReferenceSpec(
    occ_closed=("3s", "3p", "4s"), occ_open=("3d",),
    occ_e_ref=-46.12052897)  # the forced basin: 4/4 guesses (minao/atom/
# huckel/vsap) land it -- the REPRODUCIBLE choice; the unforced minao/vsap
# basin (-46.12083302, 0.30 mEh lower) exists but is guess-dependent, and
# unforced atom/huckel run away to 4s^2 4p^1 (+104.7 mEh). Ticket 19.
TI_REF = ReferenceSpec(
    occ_closed=("3s", "3p", "4s"), occ_open=("3d",),
    occ_e_ref=-57.60649802)  # good basin (window defect 1.4e-8); huckel
# and unforced-default land the +11.1 mEh member-gauge basin (defect 1.6e-2)
# -- the fact-21a phenomenology, reversed: under ECP the DEFAULT is bad.


@pytest.mark.fast
def test_ecp_wiring_valence_bookkeeping():
    """mol.nelectron is ECP-honest (valence only) and the sector geometry
    derives from it: C/Si -> 4 e (ncore=1, the explicit ns^2), Sc -> 11 e
    (ncore=4, the explicit 3s^2 3p^6)."""
    for cfg, pat, nelec, ncore in [(C_ECP, "2p", 4, 1), (SI_ECP, "3p", 4, 1)]:
        mol = build_mol(cfg)
        assert mol.nelectron == nelec
        assert not mol.cart
        assert len(drivers.select_shell_ao(mol, pat)) == 3  # complete m-set
        assert ncore_of(mol.nelectron, 2) == ncore
        assert sector_occupations(2, cfg.spin) == (2, 0)
    mol = build_mol(SC_ECP)
    assert mol.nelectron == 11  # Ne-core ECP: 3s 3p 3d 4s explicit
    assert ncore_of(mol.nelectron, 3) == 4


@pytest.mark.slow
def test_ecp_c_multiplet_exact():
    """C/ccECP p^2: the textbook table at machine precision -- 3P x3
    (<L^2> = 2) < 1D x5 (6) < 1S (0), degenerate within each term, and
    the frozen-configuration spacing ratio (1S-3P)/(1D-3P) = 5/2 (the
    one-Slater-parameter structure of a pure p^2 window)."""
    mol, mf = drivers.run_rohf(C_ECP, "results/C_ccecpccpvdz_ccecp_s2_rohf.chk")
    act = drivers.select_valence_ao(mol, ["2p"])
    c_shell = np.eye(mol.nao_nr())[:, act]
    res = drivers.run_casci(mf, 3, (1, 1), drivers.build_full_mo(mf, [0], act),
                            "results/C_ccecpccpvdz_ccecp_s2_rohf.chk",
                            nroots=9, name="casci_ms0")
    ops = operators.build_operators(mol, mf, c_shell, 1, 1, 1)
    h, lz, l2 = ops["H"], ops["Lz"], ops["L2"]
    assert np.linalg.norm(h @ lz - lz @ h) < 1e-13  # measured 2.2e-15
    terms = []
    for k in range(9):
        ci = np.asarray(res.ci[k]).ravel()
        terms.append(analysis.term_symbol(
            float(np.real(ci.conj() @ ops["S2"] @ ci)),
            float(np.real(ci.conj() @ l2 @ ci))))
    assert terms == ["3P"] * 3 + ["1D"] * 5 + ["1S"]
    ev = np.linalg.eigvalsh(h)  # OUR averaged fold: the exact one
    assert ev[2] - ev[0] < 1e-12  # measured 0.0
    assert ev[7] - ev[3] < 1e-12  # 1D degenerate, measured 0.0
    d_1d, d_1s = ev[3] - ev[0], ev[8] - ev[0]
    assert 0.04 < d_1d < 0.08    # measured 0.0586 Eh (all-electron 0.0881)
    assert abs(d_1s / d_1d - 2.5) < 1e-3  # measured 2.5000


@pytest.mark.slow
def test_ecp_raw_detector_louder():
    """The polarized-core detector fires LOUDER under ECP: with only the
    ns^2 explicit, the core sits closer to the open shell. Measured:
    ||[H,Lz]||_raw = 1.09e-4 (all-electron C: 6.1e-6) and the 3P split
    8.0e-5 Eh (AE: 4.4e-5)."""
    mol, mf = drivers.run_rohf(C_ECP, "results/C_ccecpccpvdz_ccecp_s2_rohf.chk")
    act = drivers.select_valence_ao(mol, ["2p"])
    c_shell = np.eye(mol.nao_nr())[:, act]
    ops = operators.build_operators(mol, mf, c_shell, 1, 2, 0,
                                    spherical_core=False)
    h, lz = ops["H"], ops["Lz"]
    comm = np.linalg.norm(h @ lz - lz @ h)
    assert 5e-5 < comm < 5e-3  # measured 1.09e-4
    ev = np.linalg.eigvalsh(h)
    assert 1e-5 < ev[2] - ev[0] < 1e-3  # measured 8.0e-5 Eh


@pytest.mark.slow
def test_ecp_si_multiplet_and_ncore0_fold():
    """Si/ccECP: same exact table; and the FULL-VALENCE CAS(4,4) at
    ncore=0 -- the empty-core fold must be an exact no-op (veff == 0)
    with the closure certificate intact."""
    mol, mf = drivers.run_rohf(SI_ECP, "results/Si_ccecpccpvdz_ccecp_s2_rohf.chk")
    act = drivers.select_valence_ao(mol, ["3p"])
    c_shell = np.eye(mol.nao_nr())[:, act]
    res = drivers.run_casci(mf, 3, (1, 1), drivers.build_full_mo(mf, [0], act),
                            "results/Si_ccecpccpvdz_ccecp_s2_rohf.chk",
                            nroots=9, name="casci_ms0")
    ops = operators.build_operators(mol, mf, c_shell, 1, 1, 1)
    h, lz, l2 = ops["H"], ops["Lz"], ops["L2"]
    assert np.linalg.norm(h @ lz - lz @ h) < 1e-13  # measured 5.8e-15
    ev = np.linalg.eigvalsh(h)
    assert ev[7] - ev[3] < 1e-12  # 1D degenerate (measured 0.0)
    assert 0.03 < ev[3] - ev[0] < 0.06  # measured 0.0406 Eh (AE 0.0594)

    # ncore=0: window {3s, 3p}, CAS(4,4), (3, 1)
    fv = drivers.select_valence_ao(mol, ["3s", "3p"])
    c_fv = np.eye(mol.nao_nr())[:, fv]
    h1, _g = operators.active_space_integrals(mol, mf, c_fv, 0)
    h1_direct = c_fv.conj().T @ mf.get_hcore() @ c_fv
    assert np.abs(h1 - h1_direct).max() == 0.0  # veff exactly empty
    ops0 = operators.build_operators(mol, mf, c_fv, 0, 3, 1)
    h0, lz0 = ops0["H"], ops0["Lz"]
    assert np.linalg.norm(h0 @ lz0 - lz0 @ h0) < 1e-13  # measured 1.8e-15


@pytest.mark.slow
def test_ecp_tm_references_and_basins():
    """The gated references (ticket 19 census). Ti: the good basin
    (window defect 1.4e-8). Sc: the forced basin is the reproducible one
    from EVERY guess but carries the AMPLIFIED fact-17 baseline (~1e-2:
    the explicit 3s3p core sits adjacent to the {3d,4s} window)."""
    from atomterm.operators import window_closure

    for cfg, ref, e_pinned, defect_lo, defect_hi in [
        (SC_ECP, SC_REF, -46.12052897, 1e-3, 1e-1),  # measured 1.4e-2
        (TI_ECP, TI_REF, -57.60649802, 0.0, 1e-6),   # measured 1.4e-8
    ]:
        mol, mf = drivers.run_rohf(
            cfg, f"results/{cfg.label}_rohf.chk", ref=ref)
        assert abs(mf.e_tot - e_pinned) < 1e-6
        act = drivers.select_valence_ao(mol, ["3d", "4s"])
        ncore = ncore_of(mol.nelectron, len(act))
        mo = drivers.build_full_mo(mf, list(range(ncore)), act)
        cl = window_closure(mol, mo[:, ncore:ncore + len(act)])
        assert defect_lo < cl["lz"] < defect_hi


@pytest.mark.slow
def test_ecp_tm_multiplet_tables():
    """The free-ion tables under ECP, from OUR folded H (raw shell window,
    averaged fold): Sc 2D ground (<L^2> = 6); Ti 3F x7 ground, exactly
    degenerate, with 5F x7 next -- the 120-detector sector, all integer
    L(S+1)/S(S+1)."""
    for cfg, ref, nelec, spin, ground_l2, first_ev, last_ev, l2_next in [
        (SC_ECP, SC_REF, 3, 1, 6.0, 0, 4, None),      # 2D x5 ground
        (TI_ECP, TI_REF, 4, 2, 12.0, 0, 6, 12.0),     # 3F x7, then 5F x7
    ]:
        mol, mf = drivers.run_rohf(
            cfg, f"results/{cfg.label}_rohf.chk", ref=ref)
        act = drivers.select_valence_ao(mol, ["3d", "4s"])
        c_shell = np.eye(mol.nao_nr())[:, act]
        na, nb = sector_occupations(nelec, spin)
        ncore = ncore_of(mol.nelectron, nelec)
        ops = operators.build_operators(mol, mf, c_shell, ncore, na, nb)
        h, lz, l2, s2 = ops["H"], ops["Lz"], ops["L2"], ops["S2"]
        assert np.linalg.norm(h @ lz - lz @ h) < 1e-11  # Sc 1.5e-12, Ti 5.7e-13
        # per-root labels via the library path (V^dag O V diagonal) -- an
        # einsum with the root index on the WRONG axis silently returns
        # garbage on exactly-degenerate multiplets (caught on Sc's 2D x5)
        res = analysis.analyze(h, lz, l2, s2, "sector")
        st = res.states
        assert abs(st[first_ev].l2 - ground_l2) < 1e-8
        assert abs(st[first_ev].s2 - spin * (spin + 2) / 4) < 1e-8
        assert st[last_ev].energy - st[first_ev].energy < 1e-10  # term-exact
        if l2_next is not None:  # Ti: roots 7..13 = 5F (S = 2)
            assert abs(st[7].l2 - 12.0) < 1e-8
            assert abs(st[7].s2 - 6.0) < 1e-8
            assert st[13].energy - st[7].energy < 1e-10


@pytest.mark.slow
def test_ecp_tm_purified_core_closes_the_input_windows():
    """Fact 20 under ECP at ncore=4: the purified core closes both input
    windows at machine precision (measured 4.9e-15 Sc / 7.2e-15 Ti vs the
    canonical-core defects 1.4e-2 / 9.5e-8); contamination at 1e-5-class
    foreign-l tails (measured 5.7e-5 / 7.0e-5)."""
    from atomterm.operators import window_closure

    for cfg, ref, nelec in [(SC_ECP, SC_REF, 3), (TI_ECP, TI_REF, 4)]:
        mol, mf = drivers.run_rohf(
            cfg, f"results/{cfg.label}_rohf.chk", ref=ref)
        act = drivers.select_valence_ao(mol, ["3d", "4s"])
        ncore = ncore_of(mol.nelectron, nelec)
        c_pur, contam = drivers.purified_core_orbitals(mol, mf, ncore)
        assert contam.max() < 1e-3
        for mo in [
            drivers.build_full_mo(mf, list(range(ncore)), act, c_core=c_pur),
            drivers.radial_basis(mf, list(range(ncore)), act, c_core=c_pur)[0],
        ]:
            cl = window_closure(mol, mo[:, ncore:ncore + len(act)])
            assert cl["lz"] < 1e-12 and cl["lperp"] < 1e-12
