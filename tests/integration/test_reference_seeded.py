"""The seeded-MOM mode's campaign pins: the twins' head transfers land
their declared basins EXACTLY through the product path
(:func:`atomterm.reference.run_rohf`'s seeded mode -- the projected
6-31G source carries the configuration, the MOM hook follows it, the
character post-check verifies, the basin gate anchors), and the
series members land their declared basins through the same path (the
diis+newton recipe: Newton alone crosses into the deeper sibling
basin, plain DIIS never converges).  The verify chks under ``rungs/``
are these runs' own caches (gitignored): a fresh checkout re-runs the
flow -- seconds -- into the same basins (the seed is the
reproducibility unit here, standing in for the bit-stability the
contracted bases get from their frozen references)."""

import pytest

from dataclasses import replace

from atomterm.cases import CAMPAIGN
from atomterm.reference import (ReferenceSpec, occupation_postcheck,
                                run_rohf, run_uhf, uhf_contamination,
                                uhf_seed_density)

pytestmark = pytest.mark.fast


def test_seeded_head_transfer_lands_the_ti_unc_anchor_exactly():
    cs = CAMPAIGN["Ti-unc"]
    mol, mf = run_rohf(cs.atom_config(),
                       "elements/Ti-unc/rungs/rohf_seeded_verify.chk",
                       ref=cs.ref)
    assert mf.converged
    assert abs(mf.e_tot - cs.ref.occ_e_ref) < 1e-6
    assert occupation_postcheck(mol, mf, cs.ref.occ_closed,
                                cs.ref.occ_open)


def test_seeded_head_transfer_lands_the_sc_unc_anchor_exactly():
    cs = CAMPAIGN["Sc-unc"]
    mol, mf = run_rohf(cs.atom_config(),
                       "elements/Sc-unc/rungs/rohf_seeded_verify.chk",
                       ref=cs.ref)
    assert mf.converged
    assert abs(mf.e_tot - cs.ref.occ_e_ref) < 1e-6
    assert occupation_postcheck(mol, mf, cs.ref.occ_closed,
                                cs.ref.occ_open)


def test_seeded_series_members_land_their_declared_basins():
    """Every declared non-head member of both twins' series: the fresh
    seeded run (gate + post-check armed) lands the member's declared
    basin -- the reproducibility certificate behind the @s windows."""
    for key in ("Sc-unc", "Ti-unc"):
        cs = CAMPAIGN[key]
        for spec in cs.series[1:]:
            cfg = replace(cs.atom_config(), spin=spec.spin)
            ref = ReferenceSpec(occ_closed=spec.closed, occ_open=spec.open,
                                occ_e_ref=spec.e_ref,
                                seed_basis=cs.series_seed,
                                source_closed=spec.source_closed,
                                source_open=spec.source_open,
                                conv_accelerator=spec.accel)
            mol, mf = run_rohf(
                cfg, f"elements/{key}/rungs/series/{spec.tag}.verify.chk",
                ref=ref)
            assert mf.converged, (key, spec.tag)
            assert abs(mf.e_tot - spec.e_ref) < 1e-4, (key, spec.tag)
            assert occupation_postcheck(mol, mf, spec.closed, spec.open), \
                (key, spec.tag)


def test_seeded_uhf_re_lands_the_ti_unc_anchor_basin():
    """The UHF anchor row through the seed (the fresh-SCF disease's
    UHF class: the weight hook's ranking is the disease vector, the
    pre-repair cached record sat 11.4 mEh shallow of the censused
    basin) -- the spin-split projected density seeds the fresh run,
    the gate anchors the declared basin, the spin-pure character
    reads in the contamination."""
    cs = CAMPAIGN["Ti-unc"]
    cfg = cs.atom_config()
    dm0 = uhf_seed_density(
        cfg, cs.ref, "elements/Ti-unc/rungs/uhf_seeded_verify.chk")
    mol, mf = run_uhf(cfg,
                      "elements/Ti-unc/rungs/uhf_seeded_verify.chk",
                      ref=cs.ref, dm0=dm0)
    assert mf.converged
    assert abs(mf.e_tot - cs.ref.occ_e_ref_uhf) < 1e-4
    s2, _ = uhf_contamination(mf)
    s_pure = cs.spin / 2 * (cs.spin / 2 + 1)
    assert s2 < s_pure + 0.05
