"""Tests of atomterm.drivers: orbital selection, block-safe window
construction, S-metric weights, and the fractional-occupation reference.

The build_full_mo tests are the regression guard for the ticket-09 bug:
a naive complement projector silently returned a NON-orthonormal MO matrix
whenever the raw window AOs overlapped the canonical core (Si: ~0.3) --
and everything fed that matrix (CASCI tolerated it; the CASSCF optimizer
diverged to sub-variational energies).
"""

import dataclasses
from pathlib import Path

import numpy as np
import pytest

from atomterm.system import AtomConfig
from atomterm.drivers import (
    ReferenceSpec,
    build_full_mo,
    lm_sectors,
    radial_basis,
    radial_mask,
    run_casci,
    run_frac_occ_scf,
    run_rohf,
    select_shell_ao,
    select_valence_ao,
    shell_weights,
)
pytestmark = pytest.mark.slow


def _s_ortho_err(mo, s):
    """max |C^dag S C - I| of an orbital block (S-orthonormality error)."""
    mo = np.asarray(mo)
    return float(np.abs(mo.conj().T @ s @ mo - np.eye(mo.shape[1])).max())


def test_select_valence_ao_concatenates_in_pattern_order(carbon):
    _, mol, _, act_idx = carbon
    assert select_valence_ao(mol, ["2p"]) == act_idx
    assert select_valence_ao(mol, ["2s", "2p"]) == \
        select_shell_ao(mol, "2s") + select_shell_ao(mol, "2p")


def test_select_valence_ao_rejects_unknown_pattern(carbon):
    _, mol, _, _ = carbon
    with pytest.raises(ValueError, match="9f"):
        select_valence_ao(mol, ["9f"])


def test_shell_weights_partition_unity_and_group_by_l(carbon):
    _, mol, mf, act_idx = carbon
    s = mf.get_ovlp()
    c = np.eye(mol.nao_nr())[:, act_idx]  # raw 2p AOs: exactly l = 1
    names, w = shell_weights(mol, c, s)
    assert names == ["s", "p", "d"]
    assert np.allclose(w.sum(axis=1), 1.0, atol=1e-12)
    assert np.allclose(w[:, names.index("p")], 1.0, atol=1e-12)


@pytest.fixture(scope="session")
def build_full_mo_refs(tmp_path_factory):
    """ROHF per atom, once per session (chkfiles under pytest tmp, never
    leaked across sessions -- stale caches would mask physics changes)."""
    out = []
    for cfg, patterns in [
        (AtomConfig(), ["2p"]),
        (AtomConfig("Si", "ccpvdz", spin=2), ["3p"]),
    ]:
        chk = str(tmp_path_factory.mktemp(f"bfm_{cfg.symbol}") / "rohf.chk")
        mol, mf = run_rohf(cfg, chk)
        out.append((cfg, mol, mf, patterns))
    return out


def test_reference_spec_validation():
    """ReferenceSpec owns the occ-pattern contract (fact 21a): the two
    pattern lists come as a pair, and the basin gate needs them."""
    ReferenceSpec()  # plain reference: all-None is the default run
    ReferenceSpec(occ_closed=("4s",), occ_open=("3d",))
    with pytest.raises(ValueError, match="together"):
        ReferenceSpec(occ_closed=("4s",))
    with pytest.raises(ValueError, match="needs the"):
        ReferenceSpec(occ_e_ref=-848.40648875)


def test_build_full_mo_blocks_orthonormal_and_span_preserving(build_full_mo_refs):
    """The ticket-09 regression guard, on both atoms.

    Si is the case that broke: raw 3p AOs overlap the canonical 2p core
    MOs by ~0.3, so the block-safe construction has real work to do.
    """
    for cfg, mol, mf, patterns in build_full_mo_refs:
        s = mf.get_ovlp()
        act_idx = select_valence_ao(mol, patterns)
        ncore = (mol.nelectron - 2) // 2

        mo = build_full_mo(mf, list(range(ncore)), act_idx)
        # (nao, nao) and fully S-orthonormal -- the old bug returned
        # (nao, nao+6) with err ~ 1 on Si
        assert mo.shape == (mol.nao_nr(), mol.nao_nr()), cfg.symbol
        assert _s_ortho_err(mo, s) < 1e-10, cfg.symbol

        # core SPAN preserved: the block projectors agree with the plain
        # canonical core projector (Löwdin within the span must not move it)
        c_can = mf.mo_coeff[:, :ncore]
        c_new = mo[:, :ncore]
        p_can = c_can @ c_can.conj().T @ s
        p_new = c_new @ c_new.conj().T @ s
        assert np.abs(p_new - p_can).max() < 1e-8, cfg.symbol

        # window l-purity survives the core-overlap correction
        c_win = mo[:, ncore : ncore + len(act_idx)]
        names, w = shell_weights(mol, c_win, s)
        assert w[:, names.index("p")].min() > 0.99, cfg.symbol


@pytest.fixture(scope="session")
def si_frac_chk(tmp_path_factory):
    return tmp_path_factory.mktemp("frac_si") / "frac.chk"


def test_frac_occ_scf_open_block_ordering(si_frac_chk):
    """Average-of-configuration reference for Si: the fractional block must
    land on the 3p shell -- closed [Ne]+3s below (the ordering assumption
    of run_frac_occ_scf, first nontrivial atom after carbon).

    Runs twice deliberately: first call = fresh SCF, second = chkfile cache
    hit. Both returned objects must carry the ensemble occupation hook --
    re-deriving occupations from energies must reproduce the stored mo_occ
    (the ticket-10 finding: the hook used to be missing on cache hits).
    """
    cfg = AtomConfig("Si", "ccpvdz", spin=2)
    Path(si_frac_chk).unlink(missing_ok=True)
    for k in range(2):  # 1st pass: fresh; 2nd pass: cache hit
        mol_f, mf_f = run_frac_occ_scf(
            dataclasses.replace(cfg, spin=0), str(si_frac_chk),
            nelec_open=2, n_open=3)
        assert mf_f.converged, k
        occ = np.asarray(mf_f.mo_occ)  # the occupations the SCF converged with
        n_closed = (mol_f.nelectron - 2) // 2
        assert np.allclose(occ[:n_closed], 2.0), k
        assert np.allclose(occ[n_closed : n_closed + 3], 2.0 / 3.0), k
        assert np.allclose(occ[n_closed + 3 :], 0.0), k
        # the hook is part of the object's identity on BOTH paths
        assert np.allclose(mf_f.get_occ(mf_f.mo_energy), occ), k
        # the fractional block sits on the 3p AOs (l-purity of the open shell)
        s = mf_f.get_ovlp()
        names, w = shell_weights(mol_f, mf_f.mo_coeff[:, n_closed : n_closed + 3], s)
        assert w[:, names.index("p")].min() > 0.99, k


def test_sc_rohf_reference_is_3d1_4s2_with_argon_core(sc_rohf):
    """Ticket 11's TM reference risks, pinned: the ROHF converges to the
    3d1 4s2 occupation with the canonical [:9] core = [Ar] exactly -- no
    occupation flipping, no 3d charge in the folded core."""
    cfg, mol, mf, act_idx = sc_rohf
    assert mf.converged
    occ = np.asarray(mf.mo_occ)
    ncore = 9
    assert np.allclose(occ[:10], 2.0) and occ[10] == 1.0
    # the 9 canonical core orbitals are the 9 lowest, and carry no d charge
    assert mf.mo_energy[:ncore].max() < mf.mo_energy[ncore]
    s = mf.get_ovlp()
    names, w = shell_weights(mol, mf.mo_coeff[:, :ncore], s)
    # measured ~5e-5; flaked once at the 1e-4 bound in a full-suite draw
    # (2026-10-19, 2x-headroom class per fact 10) -- 5e-4 still certifies
    # "no d charge in the core" with 3 orders of margin
    assert w[:, names.index("d")].max() < 5e-4


def test_ti_rohf_reference_needs_character_forced_occupations(ti_rohf, tmp_path):
    """Ticket 12's blocker, pinned: the DEFAULT energy-ordered ROHF on Ti is
    an aufbau runaway (a doubly-occupied 3d inside canonical [:10] and a
    {4s, 4p} open block, +123.07 mEh above the physical 3d2 4s2 3F state),
    while the character-forced hook lands the physical reference
    deterministically (E = -848.40648875, core d ~ 2e-13, open = 3d2)."""
    cfg, mol, mf, act_idx = ti_rohf
    assert mf.converged
    occ = np.asarray(mf.mo_occ)
    assert np.allclose(occ[:10], 2.0) and occ[10] == 1.0 and occ[11] == 1.0
    assert -848.41 < mf.e_tot < -848.40  # trap converges to -848.2834
    s = mf.get_ovlp()
    names, w = shell_weights(mol, mf.mo_coeff[:, :10], s)
    assert w[:, names.index("d")].max() < 1e-6  # measured ~2e-13
    names_o, w_o = shell_weights(mol, mf.mo_coeff[:, 10:12], s)
    assert w_o[:, names_o.index("d")].min() > 0.99  # two pure 3d open MOs

    # the documented failure face: default (energy-ordered) filling is
    # unusable as the reference. Like the Sc frac-occ runaway it has two
    # faces: (a) converges to the WRONG configuration (a 3d in the core),
    # (b) never converges. Either way the hook is mandatory.
    chk = tmp_path / "ti_default_rohf.chk"
    _, mf_d = run_rohf(cfg, str(chk))
    if mf_d.converged:
        names_d, w_d = shell_weights(mol, mf_d.mo_coeff[:, :10], s)
        assert w_d[:, names_d.index("d")].max() > 0.9  # 3d invaded the core
        assert mf_d.e_tot > -848.30  # +123.07 mEh above the reference
    # else: the non-converging face -- nothing to pin (the runaway
    # surface simply failed to converge); either face is unusable.


def test_sc_frac_occ_selects_open_shell_by_projection(sc_rohf, tmp_path):
    """The fact-15(c) invariant: on a TM the fractional block MUST be
    selected by character. Energy-ordered filling converges -- to a
    p-containing ensemble (4p intrudes into the 3d/4s block); the
    open_idx projection picks {3d, 4s} exactly."""
    cfg, mol, mf, act_idx = sc_rohf
    n_open = len(act_idx)

    chk = tmp_path / "frac_proj.chk"
    mol_p, mf_p = run_frac_occ_scf(
        cfg, str(chk), nelec_open=3, n_open=n_open, open_idx=act_idx)
    assert mf_p.converged
    frac_cols = np.where(np.isclose(np.asarray(mf_p.mo_occ), 0.5))[0]
    s = mf_p.get_ovlp()
    names, w = shell_weights(mol_p, mf_p.mo_coeff[:, frac_cols], s)
    wd = dict(zip(names, w.mean(axis=0)))
    assert abs(wd.get("d", 0.0) - 5 / 6) < 1e-3
    assert abs(wd.get("s", 0.0) - 1 / 6) < 1e-3
    assert wd.get("p", 0.0) < 1e-3

    chk_e = tmp_path / "frac_energy_ordered.chk"
    mol_e, mf_e = run_frac_occ_scf(
        cfg, str(chk_e), nelec_open=3, n_open=n_open)  # legacy energy order
    if mf_e.converged:  # the documented failure mode has two faces:
        # (a) it converges to the WRONG ensemble (4p intrudes into the
        # 3d/4s block), (b) the runaway surface sometimes never converges
        # (chaotic across SCF starts -- fact 10/15c). Either way the
        # legacy path is unusable as the reference.
        frac_cols = np.where(np.isclose(np.asarray(mf_e.mo_occ), 0.5))[0]
        s = mf_e.get_ovlp()
        names, w = shell_weights(mol_e, mf_e.mo_coeff[:, frac_cols], s)
        assert w[:, names.index("p")].mean() > 0.3  # p-contaminated block
    else:
        pass  # face (b): non-convergence is itself the failure mode


# --------------------------------------------------------------------------
# The radial-constraint substrate (ticket 14): (l, member)-labeled bases.
# --------------------------------------------------------------------------


def test_lm_sectors_partition_fused_shells(sc_rohf):
    """(l, member) sectors: exactly S-orthogonal, and the fused shells
    (nctr > 1, contraction-outer/m-inner) are unpacked so every AO lands in
    exactly one sector with its radial chain."""
    _, mol, mf, _ = sc_rohf
    sec = lm_sectors(mol)
    rows = [r for rs in sec.values() for r in rs]
    assert sorted(rows) == list(range(mol.nao_nr()))
    s = mf.get_ovlp()
    keys = sorted(sec)
    for i, ki in enumerate(keys):
        for kj in keys[i + 1:]:
            assert np.abs(s[np.ix_(sec[ki], sec[kj])]).max() < 1e-10
    # Sc: 6 s + 3x5 p + 3x5 d + 1x7 f radial chains per member sector
    assert len(sec[(0, 0)]) == 6 and len(sec[(1, 0)]) == 5
    assert len(sec[(2, 2)]) == 3 and len(sec[(3, 6)]) == 1


def test_radial_basis_labels_and_partners(sc_rohf):
    """radial_basis: square S-orthonormal matrix whose window slots carry
    their AO's (l, member) sector; every window slot has at least one
    same-sector rotation partner; the mask never pairs core or cross-(l,m)
    slots."""
    cfg, mol, mf, act_idx = sc_rohf
    ncore = (mol.nelectron - 3) // 2
    mo, labels = radial_basis(mf, list(range(ncore)), act_idx)
    s = mf.get_ovlp()
    assert mo.shape == (mol.nao_nr(), mol.nao_nr())
    assert np.abs(mo.conj().T @ s @ mo - np.eye(mo.shape[1])).max() < 1e-10
    row_sector = {r: k for k, rs in lm_sectors(mol).items() for r in rs}
    got = [labels[ncore + i] for i in range(len(act_idx))]
    assert got == [row_sector[a] for a in act_idx]
    mask = radial_mask(labels, ncore, len(act_idx))(
        None, mo.shape[1], ncore, len(act_idx), None)
    assert mask.sum() >= len(act_idx)  # every slot: >= 1 radial partner
    assert not mask[:ncore, :].any()   # core slots never pair
    for i, j in np.argwhere(mask):
        assert labels[i] == labels[j] is not None


def test_casscf_input_window_orthogonalization_defect(sc_rohf):
    """Fact 17: the raw {3d, 4s} shell window is exactly rotation-closed,
    but the core-ORTHOGONALIZED window every CASSCF must start from carries
    a ~1e-3 closure defect (the projection subtracts non-closed canonical
    core directions; the Sc window-core overlap is ~0.3). This -- not core
    polarization, not the optimizer -- is the relax rung's commutator
    baseline, and it is present BEFORE any optimization runs."""
    from atomterm.analysis import analyze
    from atomterm.operators import build_operators

    cfg, mol, mf, act_idx = sc_rohf
    ncore = (mol.nelectron - 3) // 2
    c_shell = np.eye(mol.nao_nr())[:, act_idx]
    mo_shell = build_full_mo(mf, list(range(ncore)), act_idx)
    mo_rad, _ = radial_basis(mf, list(range(ncore)), act_idx)
    for name, w in [("raw shell AOs", c_shell),
                    ("build_full_mo window", mo_shell[:, ncore:ncore + 6]),
                    ("radial_basis window", mo_rad[:, ncore:ncore + 6])]:
        ops = build_operators(mol, mf, w, ncore, 2, 1)
        res = analyze(ops["H"], ops["Lz"], ops["L2"], ops["S2"], name)
        if name == "raw shell AOs":
            assert res.comm_norm_h_lz < 1e-12
        else:
            # the fact-17 input-window defect; SCF-chaos-exposed tails:
            # observed 6.8e-4 .. 1.7e-2 across cold reruns (band widened
            # 2026-10-15)
            assert 1e-4 < res.comm_norm_h_lz < 3e-2, name


# --------------------------------------------------------------------------
# The purified core (ticket 13 / fact 20): the exactly L-closed core span.
# --------------------------------------------------------------------------


def _pur_cases(build_full_mo_refs, sc_rohf):
    """(symbol, mol, mf, patterns, ncore_add) rows for C / Si / Sc."""
    out = []
    for cfg, mol, mf, patterns in build_full_mo_refs:
        out.append((cfg.symbol, mol, mf, patterns, 2))
    cfg, mol, mf, act_idx = sc_rohf
    out.append((cfg.symbol, mol, mf, ["3d", "4s"], 3))
    return out


def test_purified_core_orbitals_exact_closed_span(build_full_mo_refs, sc_rohf):
    """Fact 20: the purified core span is EXACTLY rotation-closed on every
    atom class (s core C, p core Si, p+d core Sc) -- the fix for the
    fact-17 input-window baseline. Pinned: S-orthonormality, one-l purity
    per column, the one-particle closure certificate (lz AND lperp, no
    blind spots), the m-scalar density (shared radial: equal member
    blocks, zero cross-member coherences, zero anisotropy) and the tiny
    contamination / averaged-density-change scales."""
    from atomterm.drivers import purified_core_orbitals
    from atomterm.operators import core_anisotropy, window_closure

    # measured 2026-10-15 (cold chkfiles): max contamination
    # C 6.3e-6, Si 9.8e-6, Sc 5.4e-5; spherically averaged density change
    # C 4.2e-6, Si 3.6e-6, Sc 1.05e-4 (tail^2 of the core polarization)
    dens_band = {"C": 5e-5, "Si": 5e-5, "Sc": 5e-4}
    for sym, mol, mf, patterns, nelec in _pur_cases(build_full_mo_refs, sc_rohf):
        s = mf.get_ovlp()
        ncore = (mol.nelectron - nelec) // 2
        c_pur, contam = purified_core_orbitals(mol, mf, ncore)

        assert np.abs(c_pur.conj().T @ s @ c_pur - np.eye(ncore)).max() < 1e-12, sym
        names, w = shell_weights(mol, c_pur, s)
        assert (w.max(axis=1) > 1 - 1e-12).all(), sym        # one l per column
        clo = window_closure(mol, c_pur)
        assert clo["lz"] < 1e-12 and clo["lperp"] < 1e-12, (sym, clo)
        assert contam.max() < 3e-4, (sym, contam)

        # m-scalar density: shared radial made visible -- the member
        # diagonal blocks of one l are IDENTICAL, cross-member coherences
        # vanish, so the core needs no spherical averaging at all
        sec = lm_sectors(mol)
        d = c_pur @ c_pur.conj().T
        for l in range(4):
            members = sorted(m for (ll, m) in sec if ll == l)
            if len(members) < 2:
                continue
            base = d[np.ix_(sec[(l, members[0])], sec[(l, members[0])])]
            for m in members:
                blk = d[np.ix_(sec[(l, m)], sec[(l, m)])]
                assert np.abs(blk - base).max() < 1e-12, (sym, l)
            for i, mi in enumerate(members):
                for mj in members[i + 1:]:
                    cross = d[np.ix_(sec[(l, mi)], sec[(l, mj)])]
                    assert np.abs(cross).max() < 1e-12, (sym, l, mi, mj)
        assert core_anisotropy(mol, mf, ncore, c_core=c_pur) < 1e-12, sym

        # the reference moved only at the tail^2 level of the fold's own
        # scalar core (spherically averaged densities compared)
        from atomterm.operators import _spherical_average_density
        d_can = _spherical_average_density(
            mol, 2.0 * mf.mo_coeff[:, :ncore] @ mf.mo_coeff[:, :ncore].conj().T)
        d_new = _spherical_average_density(mol, 2.0 * d)
        rel = (np.linalg.norm(d_new - d_can) / np.linalg.norm(d_can))
        assert rel < dens_band[sym], (sym, rel)


def test_purified_substrate_closes_the_input_windows(build_full_mo_refs, sc_rohf):
    """Gate (b), the headline: windows orthogonalized against the purified
    core span are exactly rotation-closed (canonical core: Si 3.4e-3,
    Sc 5.7e-3 / 2.1e-2 -- the fact-17 baseline every CASSCF rung starts
    from; purified: machine precision, because the projector onto
    E (x) members commutes with L for ANY radial subspace E)."""
    from atomterm.drivers import purified_core_orbitals
    from atomterm.operators import window_closure

    for sym, mol, mf, patterns, nelec in _pur_cases(build_full_mo_refs, sc_rohf):
        ncore = (mol.nelectron - nelec) // 2
        act_idx = select_valence_ao(mol, patterns)
        ncas = len(act_idx)
        c_pur, _ = purified_core_orbitals(mol, mf, ncore)
        core_idx = list(range(ncore))
        for mo in (build_full_mo(mf, core_idx, act_idx, c_core=c_pur),
                   radial_basis(mf, core_idx, act_idx, c_core=c_pur)[0]):
            clo = window_closure(mol, mo[:, ncore:ncore + ncas])
            assert clo["lz"] < 1e-12 and clo["lperp"] < 1e-12, (sym, clo)


def test_avas_cross_validates_the_shell_window(sc_rohf, carbon):
    """AVAS (an independent ecosystem window selector) picks exactly our
    shell window's size and electron count, with a span overlapping the raw
    shell AOs (principal S-overlaps >= 0.9): the window choice is not an
    artifact of our selection machinery."""
    from pyscf.mcscf import avas

    for (cfg, mol, mf, act_idx), labels, nelec in [
        (sc_rohf, ["Sc 3d", "Sc 4s"], 3),
        (carbon, ["C 2p"], 2),
    ]:
        ncore = (mol.nelectron - nelec) // 2
        ncas, nelecas, mo = avas.avas(mf, labels, ncore=ncore)
        assert ncas == len(act_idx), labels
        assert int(np.sum(nelecas)) == nelec, labels
        s = mf.get_ovlp()
        c_raw = np.eye(mol.nao_nr())[:, act_idx]
        sv = np.linalg.svd(
            c_raw.conj().T @ s @ mo[:, ncore:ncore + ncas], compute_uv=False)
        assert sv.min() > 0.9, (labels, sv)


def test_avas_window_reproduces_the_shell_cascade(sc_rohf, carbon, tmp_path):
    """Fact 19 (AVAS verdict, VERIFIED not assumed): AVAS selects the same
    window SIZE and a span-CLOSE overlap, but the span is NOT the closed
    shell span -- and the two axes (overlap vs closure) are independent:
    on Sc a 0.9855 principal overlap hides a ~0.17-rad member tilt
    (closure defects O(0.5), CASCI commutators O(0.1..1), energy +80 mEh
    WORSE: the mean-field-inherited virtual d set is diffuse). On C the
    inherited relaxation is exactly the fact-18 theta* radial mixing
    (E = -250.8 mEh vs raw shell = the full frozen-scan gain) with
    near-m-symmetry. Raw shell AOs: overlap 1, closure 1e-15."""
    from atomterm.drivers import avas_window, span_overlap
    from atomterm.operators import window_closure
    from atomterm.stages import stage_facts

    for (cfg, mol, mf, act_idx), labels, nelec, spin in [
        (carbon, ["C 2p"], 2, 2),
        (sc_rohf, ["Sc 3d", "Sc 4s"], 3, 1),
    ]:
        ncore = (mol.nelectron - nelec) // 2
        na, nb = (nelec + spin) // 2, (nelec - spin) // 2
        eye = np.eye(mol.nao_nr())
        c_shell = eye[:, act_idx]
        s = mf.get_ovlp()

        mo, ncas, _ = avas_window(mf, labels, ncore)
        assert ncas == len(act_idx), labels
        w = mo[:, ncore:ncore + ncas]
        sv = span_overlap(c_shell, w, s)
        assert 0.9 < sv < 1 - 1e-6, (labels, sv)   # span-close, NOT identical
        clo = window_closure(mol, w)
        assert clo["lperp"] > 1e-5, (labels, clo)  # ... and rotation-defective

        cas = run_casci(mf, ncas, (na, nb), mo,
                        str(tmp_path / f"avas_{cfg.label}.chk"),
                        nroots=3, name=f"avas_{cfg.label}")
        facts = stage_facts(mol, mf, cas, w, ncore, na, nb, tag=f"avas_{cfg.label}")
        ref = run_casci(mf, ncas, (na, nb),
                        build_full_mo(mf, list(range(ncore)), act_idx),
                        str(tmp_path / f"shell_{cfg.label}.chk"),
                        nroots=3, name=f"shell_{cfg.label}")
        ref_facts = stage_facts(mol, mf, ref, c_shell, ncore, na, nb,
                                tag=f"shell_{cfg.label}")
        d_e = facts.e_ground - ref_facts.e_ground
        if cfg.label.startswith("C"):
            # the inherited ROHF radial relaxation = the fact-18 gain
            assert -0.26 < d_e < -0.24, d_e
            assert facts.analysis.comm_norm_h_lz < 1e-5
        else:
            assert 0.05 < d_e < 0.12, d_e          # diffuse virtual d set
            assert facts.analysis.comm_norm_h_lz > 0.1
            assert facts.analysis.comm_norm_h_l2 > 1.0


@pytest.mark.fast
def test_primitive_chain_and_select_primitives():
    """The exponent-walk selector (kin proxy retired): the chain IS the
    basis's radial coordinate (ascending exponent, one m-complete slot
    per primitive), bands/counts select whole primitives, and contractions
    -- which have no single exponent -- are rejected with the pointer to
    label/MO selection."""
    from pyscf import gto

    from atomterm.drivers import primitive_chain, select_primitives
    from atomterm.operators import window_closure

    mol = gto.Mole()
    mol.atom = "C 0 0 0"
    mol.basis = gto.parse("""
        C S
             8.0  1.0
        C S
             1.0  1.0
        C S
             0.1  1.0
        C P
             2.0  1.0
        C P
             0.2  1.0
        C D
             0.5  1.0
    """)
    mol.spin = 2
    mol.verbose = 0
    mol.build()

    # the chain: exponent ordering DIRECT, m-complete slots, row identity
    chain = primitive_chain(mol, 1)
    assert [e for e, _ in chain] == [0.2, 2.0]
    assert chain[0][1] == [6, 7, 8]          # the exp-0.2 p shell
    assert chain[1][1] == [3, 4, 5]          # the exp-2.0 p shell

    # band semantics: exactly the exponents inside the interval
    assert select_primitives(mol, 1, exponent_band=(0.5, 5.0)) == [3, 4, 5]
    assert select_primitives(mol, 1, exponent_band=(0.1, 0.3)) == [6, 7, 8]
    # count = the most diffuse end
    assert select_primitives(mol, 1, count=1) == [6, 7, 8]
    assert select_primitives(mol, 1, count=2) == [3, 4, 5, 6, 7, 8]

    # the error surface: both/neither selector, over-count, empty band
    with pytest.raises(ValueError, match="exactly one"):
        select_primitives(mol, 1)
    with pytest.raises(ValueError, match="exactly one"):
        select_primitives(mol, 1, count=1, exponent_band=(0.0, 9.0))
    with pytest.raises(ValueError, match="chain has"):
        select_primitives(mol, 1, count=3)
    with pytest.raises(ValueError, match="no l=1 primitive"):
        select_primitives(mol, 1, exponent_band=(100.0, 200.0))

    # the physics pin: a selected primitive m-set is exactly closed
    cl = window_closure(mol, np.eye(mol.nao_nr())[
        :, select_primitives(mol, 1, count=1)])
    assert cl["lz"] < 1e-12 and cl["lperp"] < 1e-12

    # contracted shells have no single exponent: clear error, and the
    # message points at the label/MO routes instead
    ccd = gto.Mole()
    ccd.atom = "C 0 0 0"
    ccd.basis = "cc-pvdz"
    ccd.verbose = 0
    ccd.build()
    with pytest.raises(ValueError, match="contractions"):
        primitive_chain(ccd, 1)


@pytest.mark.fast
def test_forced_occ_l_specifiers_and_open_premium():
    """The occupation hook's span vocabulary and ranking: 'l0'..'l9'
    select whole angular-momentum subspaces (the honest vocabulary on
    uncontracted bases), and the open block is ranked by 2 w_open +
    w_closed -- with whole-l spans a 1:1 sum lets closed-character
    leftovers tie the open MOs (measured on U: DIIS then closed s/p
    garbage, E > +1e3), while the premium keeps the early-cycle leakage
    case (Sc/Ti) findable."""
    from pyscf import gto

    from atomterm.reference import _forced_rohf_occ

    mol = gto.Mole()
    mol.atom = "C 0 0 0"
    mol.basis = gto.parse("""
        C S
             8.0  1.0
        C S
             1.0  1.0
        C S
             0.1  1.0
        C P
             2.0  1.0
        C P
             0.2  1.0
        C D
             0.5  1.0
    """)
    mol.spin = 2   # nelec = 6: closed = 2, open = 2
    mol.verbose = 0
    mol.build()

    from atomterm.drivers import lm_sectors
    rows_l = {l: sorted(r for (ll, m), rws in lm_sectors(mol).items()
                        if ll == l for r in rws) for l in range(4)}
    hook = _forced_rohf_occ(mol, ["l0", "l1"], ["l2"])
    mo = np.eye(mol.nao_nr())
    occ = hook(np.zeros(mol.nao_nr()), mo)
    closed = set(np.flatnonzero(occ == 2.0))
    open_ = set(np.flatnonzero(occ == 1.0))
    assert all(i in rows_l[0] + rows_l[1] for i in closed)  # s/p-pure
    assert open_ <= set(rows_l[2])                          # d wins the
    # 2:1 premium even though closed leftovers tie at 1.0
    with pytest.raises(ValueError, match="matched no AO labels"):
        _forced_rohf_occ(mol, ["4z"], ["l2"])


def test_occupied_core_idx_selects_the_closed_block_not_canonical_order():
    """The core fold takes the DOUBLY-OCCUPIED MOs; canonical MO order
    is meaningless on near-degenerate uncontracted bases (U: the first
    canonical MOs are virtual junk at eps ~ 0.02 while the closed block
    sits at 22/61/88/120/143 -- folding range(ncore) silently corrupted
    the fold at max|dE| = 14.5 Eh). The contracted-basis coincidence is
    a checked invariant, not an assumption."""
    import warnings

    from atomterm.drivers import occupied_core_idx

    class FakeMF:
        def __init__(self, mo_occ):
            self.mo_occ = np.asarray(mo_occ, dtype=float)

    # U-style: closed block scattered (canonical [0:5] are virtual junk)
    occ = [0.0] * 10
    for k in (2, 5, 6, 8, 9):
        occ[k] = 2.0
    mf = FakeMF(occ)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        core = occupied_core_idx(mf, 5)
    assert core == [2, 5, 6, 8, 9]
    assert any("NOT the first" in str(w.message) for w in rec)

    # contracted-style coincidence: closed block IS the canonical head
    mf = FakeMF([2.0, 2.0, 2.0, 1.0, 0.0])
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        assert occupied_core_idx(mf, 3) == [0, 1, 2]
    assert not any("NOT the first" in str(w.message) for w in rec)

    # asking for more closed orbitals than exist is refused loudly
    with pytest.raises(ValueError, match="closed block has"):
        occupied_core_idx(mf, 4)
