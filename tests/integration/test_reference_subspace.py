"""The subspace-ROHF pins (ticket 37 Stage 3): the confined-valence
runner's exactness at the degenerate limit, its closed-shell
cross-check against confined_core_orbitals, the runner's guard modes,
and the twins' per-contraction certificate values.  The substrates
build through the certified seeded product path (seconds); the solves
themselves are sub-second."""

import numpy as np
import pytest

from dataclasses import replace

from atomterm.cases import CAMPAIGN
from atomterm.core import closed_shell_core_energy, confined_core_orbitals
from atomterm.reference import (ReferenceSpec, _ao_family, run_rohf,
                                subspace_rohf)
from atomterm.windows import lowdin, select_valence_ao

pytestmark = pytest.mark.fast

CORE_SHELLS = ("1s", "2s", "2p", "3s", "3p")


@pytest.fixture(scope="module")
def carbon_references(tmp_path_factory):
    """(mol, mf_rohf) for the campaign's C/cc-pVDZ configuration: the
    free reference the whole-space leg must reproduce."""
    cfg = CAMPAIGN["C"].atom_config()
    root = tmp_path_factory.mktemp("c_subspace")
    return run_rohf(cfg, str(root / "rohf.chk"))


def test_subspace_whole_space_reproduces_the_free_reference_exactly(
        carbon_references):
    """The degenerate limit: subspace == the whole space == the free
    reference (exact-by-construction; the fresh-draw band: < 1e-8 Eh,
    2026-11-16)."""
    mol, mf = carbon_references
    n = mol.nao_nr()
    eye = np.eye(n)
    core1 = select_valence_ao(mol, ["1s"])
    others = [i for i in range(n) if i not in core1]
    _, mf_sub = subspace_rohf(mf, eye[:, core1], 1, eye[:, others],
                              ("l0",), ("l1",))
    assert abs(mf_sub.e_tot - mf.e_tot) < 1e-8
    assert int(round(np.sum(mf_sub.mo_occ))) == 6
    assert mf_sub.subspace["n_fam"] + mf_sub.subspace["n_win"] == n


def test_subspace_ion_matches_confined_core_orbitals():
    """The valence-empty solve == confined_core_orbitals's own
    Roothaan loop (two constructions, one confined problem)."""
    mol, mf_ion, fam_raw, e_loop, rep = _carbon_ion_problem()
    _, mf_sub = subspace_rohf(mf_ion, fam_raw, 1,
                              np.zeros((mol.nao_nr(), 0)),
                              ("l0", "l1"), ())
    assert abs(mf_sub.e_tot - e_loop) < 1e-6
    assert rep["comm"] < 1e-9


def test_subspace_ion_seedless_form_runs_under_the_newton_arm():
    """The seedless form composes with the certificate arm: the same
    confined ion under newton lands the Roothaan loop's energy (the
    fresh-draw band: < 1e-6 Eh, 2026-11-17)."""
    mol, mf_ion, fam_raw, e_loop, _rep = _carbon_ion_problem()
    _, mf_sub = subspace_rohf(mf_ion, fam_raw, 1,
                              np.zeros((mol.nao_nr(), 0)),
                              ("l0", "l1"), (), conv_accelerator="newton")
    assert abs(mf_sub.e_tot - e_loop) < 1e-6


def _carbon_ion_problem():
    """The C4+ ion's confined-one-pair problem: (mol, the bare ROHF
    substrate, the family's raw columns, the Roothaan loop's energy,
    its report)."""
    from pyscf import scf
    from atomterm.system import build_mol
    cfg = CAMPAIGN["C"].atom_config()
    # the C4+ ion (2e = the one confined pair): the nuclear field is
    # atom-blind to the charge; the counting follows mol.nelectron
    mol = build_mol(replace(cfg, charge=4, spin=0))
    n = mol.nao_nr()
    s = mol.intor_symmetric("int1e_ovlp")
    fam_rows = select_valence_ao(mol, ["1s"]) + select_valence_ao(mol, ["2s"])
    eye = np.eye(n)
    fam_raw = eye[:, fam_rows]
    fam = lowdin(fam_raw, s)
    mf_ion = scf.ROHF(mol)   # the ion-like substrate (reads no orbitals)
    c_conf, rep = confined_core_orbitals(mf_ion, fam, 1)
    e_loop = closed_shell_core_energy(mf_ion, c_conf)
    return mol, mf_ion, fam_raw, e_loop, rep


def test_subspace_refuses_a_family_narrower_than_the_pairs(
        carbon_references):
    """The family-capacity guard: the core's pairs need that many
    family columns."""
    mol, mf = carbon_references
    n = mol.nao_nr()
    eye = np.eye(n)
    core1 = select_valence_ao(mol, ["1s"])
    others = [i for i in range(n) if i not in core1]
    with pytest.raises(ValueError, match="too narrow"):
        subspace_rohf(mf, eye[:, core1], 2, eye[:, others],
                      ("l0",), ("l1",))


def test_subspace_refuses_a_core_that_overfreezes_the_atom(
        carbon_references):
    """The filling-arithmetic guard: three pairs on carbon leave a
    negative valence count against the 2-open sector."""
    mol, mf = carbon_references
    n = mol.nao_nr()
    eye = np.eye(n)
    core1 = select_valence_ao(mol, ["1s"])
    others = [i for i in range(n) if i not in core1]
    with pytest.raises(ValueError, match="over-freezes"):
        subspace_rohf(mf, eye[:, core1], 3, eye[:, others],
                      ("l0",), ("l1",))


def test_subspace_refuses_a_window_too_small_for_the_filling(
        carbon_references):
    """The window-capacity guard: the window block must hold the
    valence filling's columns."""
    mol, mf = carbon_references
    n = mol.nao_nr()
    eye = np.eye(n)
    core1 = select_valence_ao(mol, ["1s"])
    d_rows = select_valence_ao(mol, ["3d"])
    with pytest.raises(ValueError, match="window block holds"):
        subspace_rohf(mf, eye[:, core1], 1, eye[:, d_rows[:2]],
                      ("l0",), ("l1",))


def test_subspace_refuses_a_seed_the_subspace_cannot_host(
        carbon_references):
    """The capture guard: occupied columns with no carrier in the
    subspace (2s/2p over [1s family | 3d window]) refuse below 0.5."""
    mol, mf = carbon_references
    n = mol.nao_nr()
    eye = np.eye(n)
    core1 = select_valence_ao(mol, ["1s"])
    d_rows = select_valence_ao(mol, ["3d"])
    with pytest.raises(RuntimeError, match="cannot host"):
        subspace_rohf(mf, eye[:, core1], 1, eye[:, d_rows],
                      ("l0",), ("l1",))


def _twin_substrate(key, spec=None):
    """The twin's trusted reference through the seeded product path
    (the head, or one series member)."""
    cs = CAMPAIGN[key]
    cfg = cs.atom_config()
    if spec is None:
        cfg, ref = cfg, cs.ref
        tag = "head"
    else:
        cfg = replace(cfg, spin=spec.spin)
        ref = ReferenceSpec(occ_closed=spec.closed, occ_open=spec.open,
                            occ_e_ref=spec.e_ref,
                            seed_basis=cs.series_seed,
                            source_closed=spec.source_closed,
                            source_open=spec.source_open,
                            conv_accelerator=spec.accel)
        tag = spec.tag
    mol, mf = run_rohf(cfg, f"elements/{key}/rungs/subspace_verify_{tag}.chk",
                       ref=ref)
    return cs, mol, mf, tag


@pytest.mark.parametrize("key", ["Sc-unc", "Ti-unc"])
def test_twin_whole_space_legs_reproduce_their_anchors(key):
    """The twins' whole-space legs == their declared anchors (head +
    both members; the fresh-draw band: < 1e-7 Eh, 2026-11-16)."""
    cs = CAMPAIGN[key]
    for spec in (None, cs.series[1], cs.series[2]):
        cs_t, mol, mf, tag = _twin_substrate(key, spec)
        anchor = (cs_t.ref.occ_e_ref if spec is None else spec.e_ref)
        fam = _ao_family(mol, CORE_SHELLS, cs_t.shell_bands)
        n = mol.nao_nr()
        fam_rows = set(np.flatnonzero(np.any(np.abs(fam) > 1e-12,
                                             axis=1)).tolist())
        rest = [i for i in range(n) if i not in fam_rows]
        eye = np.eye(n)
        closed = cs_t.ref.occ_closed if spec is None else spec.closed
        open_ = cs_t.ref.occ_open if spec is None else spec.open
        _, mf_sub = subspace_rohf(mf, fam, cs_t.ncore, eye[:, rest],
                                  closed, open_, conv_accelerator="newton")
        assert abs(mf_sub.e_tot - anchor) < 1e-7, (key, tag)
        assert int(round(np.sum(mf_sub.mo_occ))) == mol.nelectron


@pytest.mark.parametrize(
    "key,d_s1,d_s2",
    [("Sc-unc", 21.233, 1.267), ("Ti-unc", 38.031, 2.197)])
def test_twin_head_s1_s2_monotone(key, d_s1, d_s2):
    """Anchor 7's monotonicity: dE(s2) <= dE(s1) (the measured values
    pinned, 2026-11-16)."""
    cs_t, mol, mf, _ = _twin_substrate(key)
    anchor = cs_t.ref.occ_e_ref
    fam = _ao_family(mol, CORE_SHELLS, cs_t.shell_bands)
    from atomterm.basis import series_window
    s = mf.get_ovlp()
    # the @s windows are the SERIES average (a single-configuration
    # block is rank-1 -- the construction refuses k=2 without it)
    series_mfs = [mf]
    for spec in cs_t.series[1:]:
        cfg_m = replace(cs_t.atom_config(), spin=spec.spin)
        ref_m = ReferenceSpec(occ_closed=spec.closed, occ_open=spec.open,
                              occ_e_ref=spec.e_ref,
                              seed_basis=cs_t.series_seed,
                              source_closed=spec.source_closed,
                              source_open=spec.source_open,
                              conv_accelerator=spec.accel)
        _, mf_m = run_rohf(
            cfg_m, f"elements/{key}/rungs/subspace_verify_{spec.tag}.chk",
            ref=ref_m)
        series_mfs.append(mf_m)
    de = {}
    for wtag, k in (("s1", 1), ("s2", 2)):
        win, _ = series_window(mol, series_mfs, s, ("3d", "4s"),
                               cs_t.prior_bands,
                               {"3d": k, "4s": k}, cs_t.ncore)
        _, mf_sub = subspace_rohf(mf, fam, cs_t.ncore, win,
                                  cs_t.ref.occ_closed, cs_t.ref.occ_open,
                                  conv_accelerator="newton")
        de[wtag] = (mf_sub.e_tot - anchor) * 1e3
    assert de["s2"] <= de["s1"] + 1e-6, (key, de)
    assert abs(de["s1"] - d_s1) < 0.2 and abs(de["s2"] - d_s2) < 0.2, \
        (key, de)


@pytest.mark.parametrize("key,de_ref", [("Sc-unc", 1.401), ("Ti-unc", 1.923)])
def test_twin_head_r1_certificate_value(key, de_ref):
    """The head's @r1 confinement cost (the measured certificate,
    2026-11-16; the capture exact at 1.0)."""
    cs_t, mol, mf, _ = _twin_substrate(key)
    anchor = cs_t.ref.occ_e_ref
    fam = _ao_family(mol, CORE_SHELLS, cs_t.shell_bands)
    from atomterm.basis import series_window
    s = mf.get_ovlp()
    win, _ = series_window(mol, [mf], s, ("3d", "4s"), cs_t.prior_bands,
                           {"3d": 1, "4s": 1}, cs_t.ncore)
    _, mf_sub = subspace_rohf(mf, fam, cs_t.ncore, win,
                              cs_t.ref.occ_closed, cs_t.ref.occ_open,
                              conv_accelerator="newton")
    d = (mf_sub.e_tot - anchor) * 1e3
    assert abs(d - de_ref) < 0.2, (key, d)
    assert min(mf_sub.subspace["capture"]) > 0.999


@pytest.mark.parametrize(
    "key,de_ref", [("Sc-unc", 148.625), ("Ti-unc", 199.040)])
def test_twin_pure_d_r1_certificate_value(key, de_ref):
    """The pure-d member's @r1 3d-span confinement cost: the
    head-derived single radial cannot host the d-richer configuration's
    mean field (the certificate's poverty class, 2026-11-16)."""
    cs_t, mol, mf_head, _ = _twin_substrate(key)
    spec = cs_t.series[2]
    _, mol_m, mf_m, tag = _twin_substrate(key, spec)
    fam = _ao_family(mol, CORE_SHELLS, cs_t.shell_bands)
    from atomterm.basis import series_window
    s = mf_head.get_ovlp()
    win, _ = series_window(mol, [mf_head], s, ("3d",), cs_t.prior_bands,
                           {"3d": 1}, cs_t.ncore)
    _, mf_sub = subspace_rohf(mf_m, fam, cs_t.ncore, win,
                              spec.closed, spec.open,
                              conv_accelerator="newton")
    d = (mf_sub.e_tot - spec.e_ref) * 1e3
    # the fresh draws on the cached verify-rungs reproduce the
    # certificate values to all printed digits (2026-11-17); the band
    # absorbs a fresh regeneration's draw-to-draw variation
    assert abs(d - de_ref) < 1.0, (key, tag, d)
    assert int(round(np.sum(mf_sub.mo_occ))) == mol_m.nelectron
