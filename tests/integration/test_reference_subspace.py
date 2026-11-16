"""The subspace-ROHF pins (ticket 37 Stage 3): the confined-valence
runner's exactness at the degenerate limit, its closed-shell
cross-check against confined_core_orbitals, and the twins'
per-contraction certificate values.  The substrates build through the
certified seeded product path (seconds); the solves themselves are
sub-second."""

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
    """The degenerate limit: subspace == the whole AO space, the solve
    reproduces the free reference EXACTLY (exact-by-construction: any
    deviation is a bug, not physics)."""
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
    """The closed-shell cross-check: the valence-empty solve on the
    family span vs confined_core_orbitals's own Roothaan loop -- two
    independent constructions of one confined problem."""
    from dataclasses import replace as dc_replace
    from pyscf import scf
    from atomterm.system import build_mol
    cfg = CAMPAIGN["C"].atom_config()
    # the C4+ ion (2e = the one confined pair): the nuclear field is
    # atom-blind to the charge; the counting follows mol.nelectron
    mol = build_mol(dc_replace(cfg, charge=4, spin=0))
    n = mol.nao_nr()
    s = mol.intor_symmetric("int1e_ovlp")
    fam_rows = select_valence_ao(mol, ["1s"]) + select_valence_ao(mol, ["2s"])
    eye = np.eye(n)
    fam_raw = eye[:, fam_rows]
    fam = lowdin(fam_raw, s)
    mf_ion = scf.ROHF(mol)   # the ion-like substrate (reads no orbitals)
    c_conf, rep = confined_core_orbitals(mf_ion, fam, 1)
    e_loop = closed_shell_core_energy(mf_ion, c_conf)
    _, mf_sub = subspace_rohf(mf_ion, fam_raw, 1, np.zeros((n, 0)),
                              ("l0", "l1"), ())
    assert abs(mf_sub.e_tot - e_loop) < 1e-6
    assert rep["comm"] < 1e-9


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
    """Both twins, head + pure-d member: the whole-space legs land the
    declared anchors exactly (the degenerate limit on the diseased
    bases)."""
    cs = CAMPAIGN[key]
    pure = cs.series[2]
    for spec in (None, pure):
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


@pytest.mark.parametrize("key,de_ref", [("Sc-unc", 1.401), ("Ti-unc", 1.923)])
def test_twin_head_r1_certificate_value(key, de_ref):
    """The head's @r1 3d4s confinement cost -- the measured certificate
    (fresh-draw values: Sc +1.401 / Ti +1.923 mEh, 2026-11-16); the
    window column capture is exact (1.0)."""
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
