"""The seeded-MOM mode's campaign pins: the twins' head transfers land
their declared basins EXACTLY through the product path
(:func:`atomterm.reference.run_rohf`'s seeded mode -- the projected
6-31G source carries the configuration, the MOM hook follows it, the
character post-check verifies, the basin gate anchors).  The verify
chks under ``rungs/`` are these runs' own caches (gitignored): a fresh
checkout re-runs the flow -- seconds -- into the same basins (the seed
is the reproducibility unit here, standing in for the bit-stability
the contracted bases get from their frozen references)."""

import pytest

from atomterm.cases import CAMPAIGN
from atomterm.reference import occupation_postcheck, run_rohf

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
