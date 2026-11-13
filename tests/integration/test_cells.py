"""The cell layer: key addressing, the store contract, and the CI
ordering bridge.

The bridge is the load-bearing piece: PySCF FCI vectors and
``detbasis`` enumerations order determinants differently for norb >= 4
with >= 2 electrons in a block (they coincide on carbon's norb=3
window -- why the difference went unnoticed). The certificate test
verifies the mapping through the eigen-equation against an
independently built Slater-Condon H: generic integrals admit no
permutation automorphism, so a wrong mapping cannot satisfy it.
"""

import json

import numpy as np
import pytest
import re
from pathlib import Path

from atomterm.cells import (
    CellData,
    CellKey,
    cell_chk,
    grid_keys,
    lane_keys,
    lattice_keys,
    load_cell,
    load_lane_cell,
    planned,
    run_lane_cell,
    route_keys,
    spine_keys,
    radial_refused,
    void_reason,
    window_grid,
)
from atomterm.cases import CAMPAIGN, campaign_spec, scenarios_for
from dataclasses import replace
from atomterm.detbasis import determinant_basis, slater_condon_matrix
from atomterm.operators import ci_in_detbasis_order, pyscf_ci_permutation
from atomterm.windows import span_overlap

# mixed tiers in one module (1:1 with src): the key/bridge/store tests
# are fast; the carbon run pins are slow -- marked per test, so the
# module-level mark would leak SCF fixtures into `pytest -m fast`

fast = pytest.mark.fast
slow = pytest.mark.slow


# --------------------------------------------------------------- the key

@fast
def test_address_round_trip():
    keys = [
        CellKey("C", "2p", (1, 1)),
        CellKey("C", "2p", (1, 1), core="raw"),
        CellKey("c", "2p", (2, 0), method="sa-casscf", mask=None, substrate="pur"),
        CellKey("Ti", "3d4s", (3, 1), method="sa-casscf", start="ao",
                mask="radial", substrate="pur", core="avg"),
        CellKey("Ti", "3d4s", (2, 2), method="casscf", start="ao",
                mask="virtmix"),
        CellKey("Sc", "3d4s", (1, 1), method="casscf", mask="frozen",
                core="frac"),
        CellKey("U", "3f2d", (4, 2), method="sa-casscf", start="ao"),
        CellKey("U14e", "9f10f11f12d13dno2", (4, 0), method="sa-casscf",
                start="no"),
    ]
    for k in keys:
        assert CellKey.from_address(k.address) == k
        assert k.address.count("_") == 7
    assert keys[0].address == "c_2p_1-1_casci_ao_-_pol_avg"
    assert keys[3].address == "ti_3d4s_3-1_sa-casscf_ao_radial_pur_avg"


@fast
def test_docstring_example_address_parses():
    """The module docstring's decoded example is a REAL address: it
    round-trips through ``from_address`` (the illustration it replaced
    was invalid -- a plain window under a no start, and a 'free' mask
    token the grammar does not define)."""
    import atomterm.cells as cells_mod

    m = re.search(r"^    (ti_\S+)$", cells_mod.__doc__, re.M)
    assert m, "the decoded example address is missing from the docstring"
    assert CellKey.from_address(m.group(1)).address == m.group(1)


@fast
def test_key_validation_and_canonicalization():
    # the pur-core degeneracy rule: raw on pur canonicalizes to avg
    k = CellKey("C", "2p", (1, 1), substrate="pur", core="raw")
    assert k.core == "avg"
    # frac is a different REFERENCE -- meaningless on the purified span
    with pytest.raises(ValueError, match="frac"):
        CellKey("C", "2p", (1, 1), substrate="pur", core="frac")
    # a mask needs an optimizer method
    with pytest.raises(ValueError, match="mask"):
        CellKey("C", "2p", (1, 1), method="casci", mask="radial")
    # the radial mask rotates member-labeled AO slots: ao-start only
    with pytest.raises(ValueError, match="radial mask"):
        CellKey("C", "2p", (1, 1), method="casscf", start="mo",
                mask="radial")
    with pytest.raises(ValueError, match="radial mask"):
        CellKey("C", "2pno1", (1, 1), method="casscf", start="no",
                mask="radial")
    # mo x pur: VALID since the window-projection builder (T1a) -- the
    # canonical columns are re-orthogonalized on the pur core span; the
    # pur-core canonicalization applies to it like every pur key
    k = CellKey("Ti", "3d4s", (3, 1), method="casscf", start="mo",
                substrate="pur", core="raw")
    assert k.core == "avg"
    assert k.address == "ti_3d4s_3-1_casscf_mo_-_pur_avg"
    for bad in [dict(method="fci"), dict(start="aos"), dict(mask="relax"),
                dict(substrate="polarized"), dict(core="averaged")]:
        with pytest.raises(ValueError, match="not in"):
            CellKey("C", "2p", (1, 1), **bad)
    with pytest.raises(ValueError, match="sector"):
        CellKey("C", "2p", (0, 1))
    with pytest.raises(ValueError, match="shell label"):
        CellKey("C", "3d_4s", (1, 1))
    # the window label parses to shell selectors
    assert CellKey("Ti", "3d4s", (3, 1)).shell_patterns == ["3d", "4s"]
    assert CellKey("C", "2p", (1, 1)).shell_patterns == ["2p"]


@fast
def test_no_window_grammar():
    """NO windows: '<source>no<k>' -- the group prefix is part of the
    window identity (differently-sized NO cells get different
    addresses); start and label must agree."""
    k = CellKey("Ti", "3d4sno2", (3, 1), start="no")
    assert k.no_source == ("3d4s", 2)
    assert k.shell_patterns == ["3d", "4s"]     # the SOURCE shells
    assert CellKey.from_address(k.address) == k
    assert k.address == "ti_3d4sno2_3-1_casci_no_-_pol_avg"
    with pytest.raises(ValueError, match="agree"):
        CellKey("Ti", "3d4s", (3, 1), start="no")
    with pytest.raises(ValueError, match="agree"):
        CellKey("Ti", "3d4sno2", (3, 1))        # ao start, no-tail label


@fast
def test_source_tag_grammar():
    """A prior-axis window carries its distillation source as the tag:
    '@r<k>' (the head configuration's block) / '@s<k>' (the series
    average), k uniform over the span; the tag is part of the address
    (differently-sourced windows never collide), stays out of the span
    selectors, and is NOT a core -- the one @-slot holds one or the
    other.  The NO kind of a tagged source names it through the
    '<source>no<k>' grammar unchanged."""
    k = CellKey("Sc", "3d4s@s2", (2, 1), method="casci")
    assert CellKey.from_address(k.address) == k
    assert k.address == "sc_3d4s@s2_2-1_casci_ao_-_pol_avg"
    assert k.source_tag == ("s", 2)
    assert k.core_tag is None
    assert k.shell_patterns == ["3d", "4s"]
    assert CellKey("Sc", "3d4s@r1", (2, 1)).source_tag == ("r", 1)
    # the tag is not a core: core resolution refuses it (the tag
    # routes as a source, so the key carries no core at all)
    with pytest.raises(ValueError, match="no core tag"):
        CellKey("Sc", "3d4s@s2", (2, 1)).core_shells()
    # a lane core tag is not a source
    assert CellKey("Ti", "3d3p4d@semi", (5, 3), substrate="pur"
                   ).source_tag is None
    # the NO kind of the series window: the source label carries the
    # tag, the group prefix rides after it
    kn = CellKey("Sc", "3d4s@s2no2", (2, 1), start="no")
    assert kn.no_source == ("3d4s@s2", 2)
    assert kn.source_tag == ("s", 2)
    assert kn.core_tag is None
    assert CellKey.from_address(kn.address) == kn


@fast
def test_prior_windows_declared_on_the_campaign():
    """The prior axis derives from the ONE declaration: Sc/Ti declare
    their series + prior realizations, and the planned prior windows
    parse into source tags the cells can resolve (the @s windows'
    series requirement validated at the case)."""
    for name in ("Sc", "Ti"):
        spec = CAMPAIGN[name]
        assert spec.series is not None
        tagged = [w for w in spec.planned_windows if "@" in w]
        assert tagged, f"{name} declares no prior windows"
        for w in tagged:
            k = CellKey(name, w, (1, 0), method="casci")
            assert k.source_tag is not None
            assert k.source_tag[1] <= 2
        # the untagged planned windows (the real spans) stay untagged
        for w in ("3d4d",):
            if w in spec.planned_windows:
                assert CellKey(name, w, (1, 0)).source_tag is None


@fast
def test_lane_window_grammar():
    """A lane window carries its core as a tag: the tag is part of the
    address (never collides with a free-route label), stays out of the
    span selectors, and resolves through the core vocabulary; NO
    windows take no tag (their core is the source's)."""
    k = CellKey("ti", "3d3p4d@semi", (5, 3), substrate="pur")
    assert CellKey.from_address(k.address) == k
    assert k.address == "ti_3d3p4d@semi_5-3_casci_ao_-_pur_avg"
    assert k.shell_patterns == ["3d", "3p", "4d"]
    assert k.core_tag == "semi"
    assert k.core_shells() == ("1s", "2s", "2p", "3s", "4s")
    # a NO window with a core tag is structurally unparseable (the NO
    # tail is terminal -- the '<source>no<k>' grammar admits no tag
    # after it), so the label never agrees with start='no'
    with pytest.raises(ValueError, match="agree"):
        CellKey("Ti", "3d4sno2@semi", (3, 1), start="no")
    with pytest.raises(ValueError, match="core tag"):
        CellKey("Ti", "3d4s", (3, 1)).core_shells()
    with pytest.raises(ValueError, match="neither"):
        CellKey("Ti", "3d3p@nope", (5, 3)).core_shells()


@fast
def test_atom_grammar_admits_one_variant_hyphen():
    """The campaign key is a plain token or ONE hyphenated variant
    suffix (the uncontracted twins): the hyphen rides inside the
    address's first ``_``-slot unharmed (round-trip exact), and the
    malformed forms -- leading, trailing, double hyphen, pure digits --
    are unrepresentable."""
    k = CellKey.from_address("ti-unc_3d4s@ar_3-1_casci_ao_-_pur_avg")
    assert k.atom == "Ti-unc"
    assert k.address == "ti-unc_3d4s@ar_3-1_casci_ao_-_pur_avg"
    assert CellKey("ti-unc", "3d4s", (3, 1)).atom == "Ti-unc"  # normalized
    for bad in ("Ti--unc", "-tiunc", "tiunc-", "1ti", "ti_unc", ""):
        with pytest.raises(ValueError, match="campaign key"):
            CellKey(bad, "3d4s", (3, 1))


@fast
def test_field_axis_default_keeps_every_store_name():
    """The field axis is an OPTIONAL ninth slot, appended only when
    folded: an eight-slot address parses as neutral and re-serializes
    byte-identical (the stores minted before the axis need no
    migration), and a folded cell can never collide with its neutral
    twin (the identity the folded probe ran under is now unambiguous)."""
    pre_axis = "ti-unc_3d4s@ar_3-1_casci_ao_-_pur_avg"
    k = CellKey.from_address(pre_axis)
    assert k.field == "neutral"
    assert k.address == pre_axis
    assert CellKey.from_address(k.address) == k
    folded = CellKey.from_address(pre_axis + "_folded")
    assert folded.field == "folded"
    assert CellKey.from_address(folded.address) == folded
    assert folded.address != k.address


@fast
def test_field_axis_folded_requires_tagged_window():
    """A folded key names a folded reference, and that field is
    defined BY its confined span -- an untagged window has no span to
    define it, so the combination is unrepresentable."""
    CellKey("Ti-unc", "3d4s@ar", (3, 1), substrate="pur", field="folded")
    with pytest.raises(ValueError, match="core-tagged"):
        CellKey("Ti-unc", "3d4s", (3, 1), field="folded")
    with pytest.raises(ValueError, match="not in"):
        CellKey("Ti-unc", "3d4s@ar", (3, 1), field="frozen")


@fast
def test_lane_keys_derive_from_the_declaration():
    """The lane legs and their SECTORS derive from the declaration:
    activating a shell moves its electrons into the window, so the
    sector follows the core (Ti semi grows (3,1) -> (5,3), the
    M-shell cores to (6,4); the Ar legs keep the case sector); Sc's
    five-leg family (d1: one electron fewer per M-shell leg)."""
    ti = {k.window: k.sector for k in lane_keys(CAMPAIGN["Ti"])}
    assert ti["3d4s@ar"] == (3, 1)
    assert ti["3d4d@ar"] == (3, 1)
    assert ti["3d3p@semi"] == (5, 3)
    assert ti["3d3p4d@semi"] == (5, 3)
    assert ti["3s3p3d@1s2s2p4s"] == (6, 4)
    assert ti["3d3p4s@1s2s2p3s"] == (6, 4)
    assert [(k.window, k.sector) for k in lane_keys(CAMPAIGN["Sc"])] == [
        ("3d4s@ar", (2, 1)), ("3d@ar", (2, 1)),
        ("3d3p@semi", (4, 3)),
        ("3d3p4d@semi", (4, 3)), ("3d3p4s@1s2s2p3s", (5, 4)),
        ("3s3p3d@1s2s2p4s", (5, 4))]
    # U14e's lanes: BACK at the CHEMICAL width (the case distills --
    # the lane window rides _window_columns, CAS(4,12)/(6,13), not
    # the retired 36-column band span); the (4,15) section carries
    # none (its window shares the 6p rows with the 5-pair core)
    assert [(k.window, k.sector) for k in lane_keys(CAMPAIGN["U14e"])] == [
        ("5f6d@6s6p7s", (4, 0)), ("5f6d7s@6s6p", (5, 1))]


@fast
def test_void_reason_gates_scenariod_cells():
    """The run gate: the anchor-first scenario voids EXACTLY the
    planned windows' optimizer cells (both label forms); the CASCI
    anchors and the primary grid run.  A void cell stays addressable
    (the key constructs; only the run is gated)."""
    spec = CAMPAIGN["Ti"]
    sc = scenarios_for(spec)
    # the DERIVED same-l class (scenarios_for): every window whose
    # shell composition repeats an l -- the old hand-list missed
    # 3d4p4d (caught when the derivation landed, 2026-11-11)
    same_l = ("3d4d", "3d4dno2", "3d4p4d", "3d4p4dno3",
              "3d4s4p4d", "3d4s4p4dno4")
    for k in planned(spec):
        if k.window in same_l:
            assert void_reason(k, sc) is not None, k.address
        else:
            assert (void_reason(k, sc) is not None) == \
                (k.method != "casci"), k.address
    assert all(void_reason(k, sc) is None for k in grid_keys(spec))
    reason = void_reason("ti_3d4d_3-1_casscf_ao_-_pol_avg", sc)
    assert reason.startswith("anchor-first") and "2026-11-02" in reason
    reason = void_reason("ti_3d4d_3-1_casci_ao_-_pol_avg", sc)
    assert reason.startswith("same-l-window") and "2026-11-10" in reason


NEXT_P = {"2p": "3p", "3p": "4p"}


@fast
def test_pshell_prior_axis_is_declared_on_c_and_si():
    """The p-shell prior axis mirrors the Sc grammar: the head = the
    case's own anchor, one promotion member (censused, one basin,
    occupation-safe), the widened one-higher-p chain, and the three
    prior windows guarded by prior-anchor-first (censused
    2026-11-12)."""
    for atom, shell in (("C", "2p"), ("Si", "3p")):
        spec = CAMPAIGN[atom]
        assert [m.tag for m in spec.series] == ["p2", "s1p3"]
        assert spec.series[0].e_ref == pytest.approx(
            json.loads((Path("elements") / atom
                        / "cells.json").read_text())["rohf_energy"])
        assert spec.prior_bands == {shell: (shell, NEXT_P[shell])}
        # C's @s2 is RETIRED (2026-11-13): the p2+s1p3 series average
        # is rank-1 -- the second radial eigenvalue 1.9e-4 against the
        # leading 2.5, the construction refuses.  Si's 3p average
        # carries a real second radial (5.1e-2, measured) and keeps it.
        expected = ((f"{shell}@r1", f"{shell}@s1") if atom == "C"
                    else (f"{shell}@r1", f"{shell}@s1", f"{shell}@s2"))
        assert spec.planned_windows == expected
        assert any(sc.name == "prior-anchor-first"
                   for sc in scenarios_for(spec))


@fast
def test_same_l_guard_derives_from_the_declared_windows():
    """The same-l guard is DERIVED (scenarios_for), never hand-listed
    per atom: Sc's five same-l lattice spans and Ti's three planned
    ones (incl. 3d4p4d, which the old hand-list missed) void at every
    method; the distinct-l lattice spans stay runnable."""
    for atom, windows in (("Sc", ("3d4d", "3d4s4d", "3d4s4p4d",
                                  "3d4s4p4d5d", "3d4s4p4d4f")),
                          ("Ti", ("3d4d", "3d4p4d", "3d4s4p4d"))):
        spec = CAMPAIGN[atom]
        sc = scenarios_for(spec)
        for k in lattice_keys(spec) + planned(spec):
            if k.window in windows or k.window.startswith(
                    tuple(w + "no" for w in windows)):
                assert void_reason(k, sc) is not None, k.address
        for k in lattice_keys(spec):
            if k.window == "3d" or k.window == "3d4s4p":
                assert void_reason(k, sc) is None, k.address


@fast
def test_radial_refused_covers_the_distilled_cases():
    """The construction refusal is KEY x SPEC: every radial cell of a
    distilled case (declared shell_radials -- the TM campaign) or a
    prior window (@r/@s tag) refuses; a raw-row case runs radial.
    regen's queue skips on it and the skeleton holds no place -- one
    predicate, one home."""
    for atom in ("Sc", "Sc-unc", "Ti", "Ti-unc"):
        spec = CAMPAIGN[atom]
        for k in grid_keys(spec) + planned(spec) + lattice_keys(spec):
            if k.mask == "radial":
                assert radial_refused(k, spec), k.address
    spec = CAMPAIGN["C"]
    k = CellKey("C", "2p", (2, 0), method="casscf", mask="radial",
                start="ao")
    assert not radial_refused(k, spec)


@fast
def test_void_reason_voids_every_twins_s_window_cell():
    """The twins' @s guard is LABEL-level: every cell of an @s window
    voids (whatever its sector -- the planned_sectors mirror moved the
    3d-span priors once and the sector-bearing prefix leaked 18 cells,
    caught 2026-11-11), and no @r1 cell is caught by it."""
    for atom in ("Sc-unc", "Ti-unc"):
        spec = CAMPAIGN[atom]
        sc = scenarios_for(spec)
        for k in planned(spec):
            tag = k.window.partition("@")[2]
            voided = void_reason(k, sc) is not None
            if tag.startswith("s"):
                assert voided, k.address
            elif tag.startswith("r"):
                assert not voided or "series-source" not in (
                    void_reason(k, sc) or ""), k.address


@fast
def test_u14e_chemistry_declaration():
    """The U14e render chemistry lives in the declaration: the window
    display names (chemical shells -- the uncontracted basis counts
    radial primitives), the core shell sets (the (5,1) promoted
    window cores 6s6p -- four pairs -- the (4,0) ones 6s6p7s), the
    chemical filling note, and the primary-family suppression."""
    spec = CAMPAIGN["U14e"]
    assert spec.window_names["10f13d"] == "5f6d"
    assert spec.window_names["10f13d7s"] == "5f6d7s"
    assert spec.window_names["10f13d14p"] == "5f6d6p"
    assert spec.core_chemistry == {5: "Rn·7s", 4: "Rn"}
    assert spec.render_primary is False
    assert spec.anchor_note.startswith("forced filling: 6s²6p⁶7s²")


@fast
def test_void_reason_voids_the_u14e_optimizer_class():
    """U14e's widened gate: the multi-solution class is void at the key
    on EVERY window (the primary included -- ``run_casscf=False``'s
    spirit made explicit) and the frac cores are void with it (the
    frac reference has no converged basin on the ECP atom); exactly
    the 36 remaining CASCI cells run (the 20 anchors + the kind
    block's mo/no cells at both ruler cores, 4 per window -- the tier
    families retired from the declaration 2026-11-06; the kind
    amendment 2026-11-09)."""
    spec = CAMPAIGN["U14e"]
    sc = scenarios_for(spec)
    keys = grid_keys(spec) + planned(spec)
    assert len(keys) == 184   # the virtual-mix keys retired (-4)
    for k in keys:
        void = void_reason(k, sc) is not None
        expect = k.method != "casci" or k.core == "frac"
        assert void == expect, k.address
    assert sum(not void_reason(k, sc) for k in keys) == 36


@fast
def test_spine_keys_contract():
    """The spine = README · Methods verbatim: the three cores; frozen x
    purified under both objectives (the certified class); one free row
    per substrate (the basin detectors)."""
    keys = spine_keys("Ti", "3d4s", (3, 1))
    assert [k.address for k in keys] == [
        "ti_3d4s_3-1_casci_ao_-_pol_raw",
        "ti_3d4s_3-1_casci_ao_-_pol_avg",
        "ti_3d4s_3-1_casci_ao_-_pol_frac",
        "ti_3d4s_3-1_casscf_ao_frozen_pur_avg",
        "ti_3d4s_3-1_sa-casscf_ao_frozen_pur_avg",
        "ti_3d4s_3-1_casscf_ao_-_pol_avg",
        "ti_3d4s_3-1_sa-casscf_ao_-_pur_avg",
    ]


@fast
def test_window_grid_is_the_one_declaration():
    """The declared key-set derived from CaseSpec (T1b): the full grid in
    render order, the spine a subset of it, the named groups filters.
    mo x pur IS in the grid (the projection builder, T1a); radial stays
    ao-only; the NO start addresses the full-span group prefix with
    k = the window's shell count."""
    from atomterm.cases import CAMPAIGN

    spec = CAMPAIGN["Ti"]
    grid = grid_keys(spec)
    addrs = [k.address for k in grid]
    # the manifest count: 29 established + the 4 mo x pur cells
    # + the moc block (2 CASCI + 8 optimizer keys)
    # + the kind-completion block (mo/no CASCI at both ruler cores)
    # -- the virtual-mix diagnostic retired (2026-11-12): 47 -> 46
    assert len(grid) == 46
    assert len(set(addrs)) == 46
    assert addrs[:4] == [           # the CASCI core block, pur third
        "ti_3d4s_3-1_casci_ao_-_pol_raw",
        "ti_3d4s_3-1_casci_ao_-_pol_avg",
        "ti_3d4s_3-1_casci_ao_-_pur_avg",
        "ti_3d4s_3-1_casci_ao_-_pol_frac",
    ]
    assert addrs[4:10] == [         # the kind block at each ruler core
        "ti_3d4s_3-1_casci_mo_-_pol_avg",
        "ti_3d4s_3-1_casci_moc_-_pol_avg",
        "ti_3d4sno2_3-1_casci_no_-_pol_avg",
        "ti_3d4s_3-1_casci_mo_-_pur_avg",
        "ti_3d4s_3-1_casci_moc_-_pur_avg",
        "ti_3d4sno2_3-1_casci_no_-_pur_avg",
    ]
    mopur = [a for a in addrs if a.split("_")[4] == "mo"
             and a.split("_")[6] == "pur"]
    assert mopur == [                   # method-major, like every block
        "ti_3d4s_3-1_casci_mo_-_pur_avg",
        "ti_3d4s_3-1_casscf_mo_-_pur_avg",
        "ti_3d4s_3-1_casscf_mo_frozen_pur_avg",
        "ti_3d4s_3-1_sa-casscf_mo_-_pur_avg",
        "ti_3d4s_3-1_sa-casscf_mo_frozen_pur_avg",
    ]
    # no radial x mo/no keys ever (void at the key)
    assert not [a for a in addrs if "_radial_" in a and "_ao_" not in a]
    # the NO windows: full-span group prefixes, k = shell count
    assert "ti_3d4sno2_3-1_casscf_no_-_pol_avg" in addrs
    # the spine is a subset view; the virtual-mix diagnostic is
    # UNDECLARED (retired 2026-11-12 -- the rows omitted, the mask
    # architecturally possible: the key constructs)
    spine = {k.address for k in spine_keys("Ti", "3d4s", (3, 1))}
    assert spine <= set(addrs)
    assert not [k for k in grid if k.mask == "virtmix"]
    # planned: the declared windows' full grids -- the real spans at
    # the primary's shape, the prior-axis tagged windows one SHORT
    # (they carry no frac core), each with its NO sibling label
    pl = planned(spec)
    assert len(pl) == 2 * len(grid) + 6 * (len(grid) - 1)
    assert {k.window for k in pl} == {
        "3d4d", "3d4dno2", "3d4s4p4d", "3d4s4p4dno4",
        "3d@r1", "3d@r1no1", "3d@s1", "3d@s1no1", "3d@s2", "3d@s2no1",
        "3d4s@r1", "3d4s@r1no2", "3d4s@s1", "3d4s@s1no2",
        "3d4s@s2", "3d4s@s2no2"}
    assert all(k.address not in addrs for k in pl)   # disjoint windows


@fast
def test_planned_sector_override_shrinks_the_fold():
    """A planned window's declared sector overrides the case default in
    its whole grid -- the promoted-pair windows' core count follows the
    sector (ncore_of(14, 6) = 4: the full core minus the 7s pair);
    windows without the override stay at the default, and the NO-form
    keys inherit their window's sector."""
    keys = planned(CAMPAIGN["U14e"])
    by_window = {}
    for k in keys:
        by_window.setdefault(k.window, []).append(k)
    assert all(k.sector == (4, 0) for k in by_window["10f13d"])
    assert all(k.sector == (4, 0) for k in by_window["10f13d14p"])
    assert all(k.sector == (5, 1) for k in by_window["10f13d7s"])
    # the tier windows (no override) and the NO-forms keep their
    # window's sector
    assert all(k.sector == (4, 0)
               for w, ks in by_window.items() if w.startswith("9f")
               for k in ks)
    no7 = [k for k in by_window if k.startswith("10f13d7sno")]
    assert no7 and all(k.sector == (5, 1)
                       for k in by_window[no7[0]])


@fast
def test_lattice_keys_derive_from_the_declaration():
    """The lattice derives from CaseSpec.lattice_windows at the
    declared lattice_sector, carrying the SAME factor grid as the
    primary window (one grammar home: core groups, kinds,
    treatments -- the pending placeholders hold the unrun places);
    undeclared atoms derive nothing; the sector is a declaration,
    not a constant (Sc's valence series differs from Ti's
    free-ion series)."""
    from atomterm.cases import CAMPAIGN, campaign_spec
    from atomterm.cells import window_grid

    assert lattice_keys(CAMPAIGN["C"]) == []
    ti = lattice_keys(CAMPAIGN["Ti"])
    # the full grammar per span: the core-field controls, both
    # substrates, every kind, treatments
    assert {k.core for k in ti} == {"raw", "avg"}
    assert {k.substrate for k in ti} == {"pol", "pur"}
    assert {k.start for k in ti} == {"ao", "mo", "moc", "no"}
    assert {"casci", "casscf", "sa-casscf"} <= {k.method for k in ti}
    # the valence series at its own sector: CAS(3,·), Ms = 1/2
    sc = lattice_keys(CAMPAIGN["Sc"])
    assert sc[0].address == "sc_3d_2-1_casci_ao_-_pol_raw"
    assert {k.sector for k in sc} == {(2, 1)}
    # the grammar IS the primary's (window_grid over a sector view,
    # frac excluded at the source)
    from dataclasses import replace as _re
    view = _re(campaign_spec("Sc"), planned_sectors={"3d": (2, 1)})
    per_span = {k.address for k in sc if k.window in ("3d", "3dno1")}
    grid = {k.address for k in window_grid(view, "3d")
            if k.core != "frac"}
    assert per_span == grid


# ------------------------------------------------------------- the bridge

SECTORS = [(3, 1, 1), (4, 2, 1), (5, 3, 1), (6, 3, 2), (6, 4, 2)]


@fast
def test_ci_permutation_is_a_permutation():
    from math import comb

    for norb, na, nb in SECTORS:
        perm = pyscf_ci_permutation(norb, na, nb)
        assert len(perm) == comb(norb, na) * comb(norb, nb)
        assert np.array_equal(np.sort(perm), np.arange(len(perm)))


@fast
def test_orders_coincide_only_in_trivial_sectors():
    """The empirical fact that kept the bug invisible: carbon's norb=3
    window forces pyscf's integer-ascending string order to coincide
    with detbasis' lexicographic order; norb >= 4 breaks it."""
    assert np.array_equal(pyscf_ci_permutation(3, 1, 1), np.arange(9))
    assert not np.array_equal(pyscf_ci_permutation(4, 2, 1), np.arange(24))


def _random_eri(norb, rng):
    """A generic chemist-notation eri with the full 8-core symmetry
    (average over the orbit under the pair transpositions) -- PySCF's
    contraction assumes it."""
    a = rng.standard_normal((norb,) * 4)
    seen, fr = [], [a]
    while fr:
        x = fr.pop()
        if any(np.array_equal(x, s) for s in seen):
            continue
        seen.append(x)
        for p in ((1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)):
            fr.append(np.ascontiguousarray(x.transpose(p)))
    eri = sum(seen) / len(seen)
    for p in ((1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)):
        assert np.allclose(eri, eri.transpose(p))
    return eri


@fast
@pytest.mark.parametrize("norb,na,nb", SECTORS)
def test_bridge_certified_by_eigen_equation(norb, na, nb):
    """THE bridge certificate: a PySCF-solved CI vector, mapped into
    detbasis order, satisfies the eigen-equation of the independently
    built Slater-Condon Hamiltonian. No ordering assumption anywhere --
    generic integrals admit no permutation automorphism."""
    from pyscf import fci
    from atomterm.operators import spin_block_one_body, spin_block_two_body

    rng = np.random.default_rng(20261027 + 10 * norb + na)
    h1 = rng.standard_normal((norb, norb))
    h1 = (h1 + h1.T) / 2
    eri = _random_eri(norb, rng)

    e_fci, ci = fci.direct_spin1.kernel(h1, eri, norb, (na, nb),
                                        conv_tol=1e-13)
    dets = determinant_basis(norb, na, nb)
    # the owners' spin-orbital expansion; detbasis eats the PHYSICIST
    # arrangement of the same eri; pyscf takes the chemist one
    ham = slater_condon_matrix(
        dets, one_body=spin_block_one_body(h1),
        two_body=spin_block_two_body(eri.transpose(0, 2, 1, 3)))

    w = ci_in_detbasis_order(ci, norb, na, nb)
    e_mine = float(np.real(w.conj() @ ham @ w))
    assert abs(e_mine - e_fci) < 1e-9
    resid = ham @ w - e_fci * w
    assert np.linalg.norm(resid) / np.linalg.norm(ham) < 1e-10


@fast
def test_configuration_summary_uses_the_bridge(tmp_path):
    """The consumer fix: a pyscf-ordered unit amplitude must label the
    determinant at ITS string address, not dets[k]."""
    from pyscf import gto
    from pyscf.fci import cistring
    from atomterm.analysis import configuration_summary

    mol = gto.M(atom="C 0 0 0", basis="cc-pvdz", spin=2, verbose=0)
    s = mol.intor_symmetric("int1e_ovlp")
    # a mixed-l window (2s + 2p) so column identity matters
    cols = [1, 2, 3, 4]
    c_win = np.eye(mol.nao_nr())[:, cols]
    na, nb = 2, 1
    # pyscf address 5 of the (2,1)/(4) sector -> its determinant
    a_strs = cistring.make_strings(range(4), na)
    b_strs = cistring.make_strings(range(4), nb)
    ia, ib = divmod(5, len(b_strs))
    sa, sb = a_strs[ia], b_strs[ib]
    occ_a = [i for i in range(4) if sa >> i & 1]
    occ_b = [i for i in range(4) if sb >> i & 1]
    # expected label: the l-mass of exactly those columns
    from atomterm.windows import shell_weights
    names, wl = shell_weights(mol, c_win, s)
    occ_cols = sorted(occ_a + occ_b)
    mass = {n: sum(wl[c, k] for c in occ_cols)
            for k, n in enumerate(names)}
    expected = "".join(f"{n}^{mass[n]:.2f}" for n in names if mass[n] > 0.01
                       ) + " (|c|=1.00)"
    v = np.zeros(len(a_strs) * len(b_strs))
    v[5] = 1.0
    [label] = configuration_summary(mol, s, [v], c_win, na, nb, top=1)
    assert label == expected


# -------------------------------------------------------- the store layer

@fast
def test_load_cell_synthetic_roundtrip(tmp_path):
    from pyscf import lib

    key = CellKey("Xe", "5p", (1, 1))  # any atom label; synthetic store
    chk = cell_chk(key, tmp_path)
    chk.parent.mkdir(parents=True)
    rng = np.random.default_rng(7)
    mo = rng.standard_normal((6, 6)) + 1j * rng.standard_normal((6, 6))
    energies = np.array([-1.0, -0.5, -0.25])
    ci = [rng.standard_normal(9), rng.standard_normal(9)]
    for name, val in (("energies", energies),
                      ("ci_stack", np.asarray(ci)),
                      ("core_energy", 5.25), ("mo_coeff", mo),
                      ("ncas", 3)):
        lib.chkfile.save(str(chk), f"atomterm/{key.address}/{name}", val)

    data = load_cell(key, tmp_path)
    assert isinstance(data, CellData)
    assert data.key == key
    assert np.array_equal(data.mo_coeff, mo)
    assert np.array_equal(data.energies, energies)
    assert data.core_energy == 5.25
    assert data.ncas == 3
    assert len(data.ci) == 2
    assert np.allclose(data.ci[1], ci[1])


@fast
def test_load_cell_prefers_sa_casscf_roots_and_reports_missing(tmp_path):
    from pyscf import lib

    key = CellKey("Xe", "5p", (1, 1), method="sa-casscf")
    chk = cell_chk(key, tmp_path)
    chk.parent.mkdir(parents=True)
    for name in (key.address, f"{key.address}_roots"):
        lib.chkfile.save(str(chk), f"atomterm/{name}/energies",
                         np.array([-2.0 if "roots" in name else -1.0]))
        lib.chkfile.save(str(chk), f"atomterm/{name}/ci_stack",
                         np.zeros((1, 9)))
        lib.chkfile.save(str(chk), f"atomterm/{name}/core_energy", 1.0)
    lib.chkfile.save(str(chk), f"atomterm/{key.address}/mo_coeff",
                     np.eye(6))
    assert load_cell(key, tmp_path).energies[0] == -2.0

    with pytest.raises(FileNotFoundError, match="run it first"):
        load_cell(CellKey("Xe", "5p", (2, 0)), tmp_path)


# ------------------------------------------------- the run pins (carbon)

def _same_facts(a, b, tol=1e-12):
    """Two facts objects describe the same run: energies, meters, labels."""
    assert a.tag == b.tag
    assert np.allclose(a.energies, b.energies, atol=tol), (a.energies, b.energies)
    assert np.allclose(a.energies_pyscf, b.energies_pyscf, atol=tol)
    for f in ("comm_norm_h_lz", "max_de_vs_pyscf",
              "ortho_err", "n_dev", "win_lz", "win_lperp",
              "ground_l2", "ground_lz2"):
        assert abs(getattr(a, f) - getattr(b, f)) <= tol, (f, getattr(a, f),
                                                           getattr(b, f))
    assert a.ground_term == b.ground_term
    assert a.terms == b.terms
    # configuration labels modulo the |c| amplitudes: inside a degenerate
    # cluster the Davidson gauge is chaotic even on pinned inputs (the
    # run_casci straggle note) -- the l-mass content is the invariant
    strip = lambda cs: [re.sub(r" \(\|c\|=[0-9.]+\)", "", c) for c in cs]
    assert strip(a.configurations) == strip(b.configurations)
    if a.core_anisotropy is not None or b.core_anisotropy is not None:
        assert a.core_anisotropy == pytest.approx(b.core_anisotropy, abs=1e-12)
    assert a.analysis.comm_norm_h_l2 == pytest.approx(
        b.analysis.comm_norm_h_l2, abs=tol)


@slow
def test_run_cell_matches_direct_path(carbon, tmp_path):
    """run_cell(casci cell) == the direct run_casci + stage_facts path --
    the same construction the rungs registry uses, byte-comparable on a
    pinned reference (CASCI on pinned MOs is deterministic)."""
    from atomterm.casscf import run_casci
    from atomterm.stages import stage_facts
    from atomterm.system import ncore_of
    from atomterm.windows import build_full_mo, occupied_core_idx
    from atomterm.cells import run_cell, cell_payload

    _, mol, mf, act_idx = carbon
    na = nb = 1
    ncore = ncore_of(mol.nelectron, 2)
    mo = build_full_mo(mf, occupied_core_idx(mf, ncore), act_idx)
    direct = run_casci(mf, len(act_idx), (na, nb), mo,
                       str(tmp_path / "direct.chk"), nroots=9, name="casci_avg")
    key = CellKey("C", "2p", (na, nb))       # casci/ao/-/pol/avg
    facts_direct = stage_facts(
        mol, mf, direct, direct.mo_coeff[:, ncore:ncore + len(act_idx)],
        ncore, na, nb, tag=key.address, spherical_core=True,
        anisotropy_ref=(mf, None))

    facts = run_cell(key, mf, workdir=tmp_path)
    assert facts.tag == key.address == "c_2p_1-1_casci_ao_-_pol_avg"
    _same_facts(facts, facts_direct)
    payload = cell_payload(key, facts)
    assert payload["cell"] == key.address

    # the cache hit re-analyzes the stored run -- same facts again
    _same_facts(run_cell(key, mf, workdir=tmp_path), facts)


@slow
def test_analyze_cell_matches_insession(carbon, tmp_path):
    """THE store certificate: analyze_cell (load_cell + facts) reproduces
    the in-session run exactly -- energies, meters, AND the stored-CI
    configurations (the bridge applied to real stored data)."""
    from atomterm.cells import analyze_cell, run_cell

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1), core="raw")   # the detector cell
    insession = run_cell(key, mf, workdir=tmp_path, fresh=True)
    reloaded = analyze_cell(key, mol, mf, workdir=tmp_path)
    _same_facts(reloaded, insession)
    assert reloaded.core_anisotropy == insession.core_anisotropy


@slow
def test_posthoc_properties_from_stored_run(carbon, tmp_path):
    """Post-hoc physics from the STORED run alone (load_cell + the
    bridge): the ground root's <L^2>/<S^2> from the stored PySCF CI and
    the energy from the stored window reproduce the payload. Gauge-safe
    observables only -- <Lz^2> is member-gauge-dependent inside a
    degenerate cluster and is NOT compared."""
    from atomterm.cells import cell_payload, load_cell, run_cell
    from atomterm.operators import build_operators
    from atomterm.system import ncore_of

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1))
    facts = run_cell(key, mf, workdir=tmp_path)
    data = load_cell(key, tmp_path)
    ncore = ncore_of(mol.nelectron, 2)
    c_act = data.mo_coeff[:, ncore:ncore + 3]
    ops = build_operators(mol, mf, c_act, ncore, 1, 1)
    w = ci_in_detbasis_order(data.ci[0], 3, 1, 1)
    l2 = float(np.real(w @ ops["L2"].real @ w))
    s2 = float(np.real(w @ ops["S2"].real @ w))
    e = float(np.real(w @ ops["H"] @ w)) + data.core_energy
    p = cell_payload(key, facts)
    assert abs(l2 - p["ground_l2"]) < 1e-9      # 2.0: gauge-safe in the 3P cluster
    assert abs(s2 - 2.0) < 1e-9                 # S=1 across the whole sector
    assert abs(e - p["e_ground"]) < 1e-9


@slow
def test_run_cell_no_start_same_subspace(carbon, tmp_path):
    """The NO start, end to end: on carbon's 2p source the term-averaged
    density's occupation groups form ONE 3-column m-set, so the no1
    window IS the source span in the NO gauge -- the CASCI must be
    basis-invariant (exact diagonalization inside one subspace; fact:
    CASCI is start-irrelevant, window-dependent)."""
    from atomterm.cells import load_cell, run_cell

    _, mol, mf, _ = carbon
    ao = run_cell(CellKey("C", "2p", (1, 1)), mf, workdir=tmp_path)
    no = run_cell(CellKey("C", "2pno1", (1, 1), start="no"), mf,
                  workdir=tmp_path)
    assert load_cell(CellKey("C", "2pno1", (1, 1), start="no"),
                     tmp_path).ncas == 3
    assert np.allclose(ao.energies, no.energies, atol=1e-9)
    assert no.win_lz == pytest.approx(0.0, abs=1e-12)   # m-set: closed span
    assert no.win_lperp == pytest.approx(0.0, abs=1e-12)


@fast
def test_sa_budget_longer_than_sector_falls_back(carbon, tmp_path):
    """A declared sa_weights block LONGER than the sector cannot host
    its roots (the occupation-group windows: CAS(3,2) cannot host a
    term's 2L+1) -- the run falls back to the MEASURED ground block
    (the default rule, window-local) instead of crashing the solver."""
    from atomterm.cells import run_cell

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1), method="sa-casscf")
    facts = run_cell(key, mf, workdir=tmp_path,
                     budgets={"sa_weights": [1 / 11] * 11})
    n = len(facts.energies_pyscf)
    assert 1 < n < 11        # the measured block anchors, not the declared


@slow
def test_sa_casscf_run_completes_the_roots_store(carbon, tmp_path):
    """store_name= pins the WRITER side of the SA store (the aa21071
    bug): the SA run returns the per-ROOTS pass, and the cell-layer
    completion must land under ``<address>_roots`` -- completing
    ``<address>`` instead pairs the ensemble's energies with the
    roots' CI, and the loader's preferred store then fails loud
    (missing ci_stack)."""
    from pyscf import lib

    from atomterm.cells import run_cell

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1), method="sa-casscf")
    facts = run_cell(key, mf, workdir=tmp_path)
    data = load_cell(key, tmp_path)          # prefers <address>_roots
    assert np.allclose(data.energies, facts.energies_pyscf, atol=1e-10)
    chk = cell_chk(key, tmp_path)
    ci = lib.chkfile.load(str(chk), f"atomterm/{key.address}_roots/ci_stack")
    assert ci is not None and len(ci) == len(facts.energies_pyscf)
    # the optimizer's macro-cycle count persists (the callback is the
    # instrument -- PySCF stores no count); CI runs read 0
    assert data.cycles > 0
    ci_key = CellKey("C", "2p", (1, 1), method="casci")
    run_cell(ci_key, mf, workdir=tmp_path)
    assert load_cell(ci_key, tmp_path).cycles == 0


@slow
def test_kernel_cell_payload_meters_match_the_dense_path(carbon, tmp_path):
    """The apply-based meters of a kernel-only payload read the STORED
    window through the run's own core resolution: the ground ⟨L²⟩
    agrees with the dense analysis payload's, and the state residual
    is present at the certificate floor (the free route's big-cell
    columns; here pinned at a small cell where both paths live)."""
    from atomterm.casscf import run_casci
    from atomterm.cells import (cell_payload, kernel_cell_payload,
                                run_cell)
    from atomterm.windows import select_valence_ao

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1), method="casci",
                  substrate="pur", core="avg")
    facts = run_cell(key, mf, workdir=tmp_path)
    dense = cell_payload(key, facts)
    kern = kernel_cell_payload(key, mol, mf, workdir=tmp_path)
    assert kern["energies_pyscf"] == dense["energies_pyscf"]
    assert abs(kern["ground_l2"] - dense["ground_l2"]) < 1e-9
    # observed on the frozen carbon: 1.2e-16 / 5.0e-16 (the apply floor)
    assert 0.0 <= kern["res_h_lz"] < 1e-8
    assert 0.0 <= kern["res_h_l2"] < 1e-8
    # no substrate given: the kernel-only shape (NaN meters, no
    # residual) -- the free-standing energies class
    bare = kernel_cell_payload(key, workdir=tmp_path)
    assert "res_h_lz" not in bare and "res_h_l2" not in bare


@slow
def test_frozen_mask_applies_to_the_mo_start(carbon, tmp_path):
    """Masks are START-INDEPENDENT: a frozen-window
    MO-start cell must converge WITH its start columns -- if the mask
    resolution slipped back into the ao branch, this cell would run
    free and the window would rotate away (carbon's radial relaxation,
    the theta* gain)."""
    from atomterm.cells import run_cell
    from atomterm.system import ncore_of

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1), method="casscf", start="mo",
                  mask="frozen")
    run_cell(key, mf, workdir=tmp_path)
    ncore = ncore_of(mol.nelectron, 2)
    start = mf.mo_coeff[:, ncore:ncore + 3]
    conv = load_cell(key, tmp_path).mo_coeff[:, ncore:ncore + 3]
    s = mf.get_ovlp()
    for i in range(3):
        # |<start_i|S|conv_i>| == 1: the same column up to phase
        assert abs(start[:, i] @ (s @ conv[:, i])) == pytest.approx(
            1.0, abs=1e-8), f"window column {i} moved under a frozen mask"


@slow
def test_run_cell_mo_start_on_pur_substrate(carbon, tmp_path):
    """The mo x pur cell end to end (T1a): canonical window+virtual
    columns re-orthogonalized on the pur core through build_full_mo's
    one treatment. On carbon the s-core x p-window are exactly
    S-orthogonal (different l) -- the projection is the no-op case; the
    defect-scale case is Ti's, measured when the gated batch runs.
    Machinery pins: the S-orthonormal substrate, the pur core block,
    CASCI internal consistency."""
    from atomterm.cells import load_cell, run_cell
    from atomterm.core import purified_core_orbitals
    from atomterm.system import ncore_of

    _, mol, mf, _ = carbon
    key = CellKey("C", "2p", (1, 1), start="mo", substrate="pur")
    facts = run_cell(key, mf, workdir=tmp_path)
    assert facts.ortho_err < 1e-10
    assert facts.max_de_vs_pyscf < 1e-9
    ncore = ncore_of(mol.nelectron, 2)
    c_pur, _ = purified_core_orbitals(mol, mf, ncore)
    data = load_cell(key, tmp_path)
    s = mf.get_ovlp()
    err = np.abs(data.mo_coeff.conj().T @ s @ data.mo_coeff
                 - np.eye(mol.nao_nr())).max()
    assert err < 1e-10
    # substrate = the core SPAN (PySCF canonicalizes the inactive block
    # WITHIN its span after the kernel (span-preserving; it hits
    # the ao-pur cells identically): the stored core spans the pur core
    assert span_overlap(data.mo_coeff[:, :ncore], c_pur, s) == \
        pytest.approx(1.0, abs=1e-10)


@fast
def test_pick_no_window_selects_by_group_number():
    """k counts GROUPS (window identity), but no_group_windows returns
    a COMPACT list -- prefixes too small to open the sector are
    skipped, so labels can be non-contiguous. Selection is by label,
    never by index (the Ti {3d,4s}/(3,1) lattice failure class: the
    1-column s prefix skipped, no2 = the full 6-column span)."""
    from atomterm.cells import _pick_no_window
    from atomterm.windows import Window

    wins = [Window("no2", None, np.eye(6, 6), True),
            Window("no4", None, np.eye(6, 6), True)]   # no1/no3 skipped
    assert _pick_no_window(wins, 2).label == "no2"
    assert _pick_no_window(wins, 4).c_win.shape[1] == 6
    with pytest.raises(ValueError, match="usable group prefixes"):
        _pick_no_window(wins, 1)
    with pytest.raises(ValueError, match="no3 is not among them"):
        _pick_no_window(wins, 3)


@slow
def test_run_lane_cell_carbon_roundtrip(carbon_d2h, tmp_path):
    """The lane runner on the frozen D2h carbon twin: the declared leg
    hosts its term exactly (all 2L+1 components read L(L+1), S(S+1),
    degenerate), the payload round-trips the store, and the runner
    refuses non-lane keys, wrong sectors, and absent core shells."""
    cfg, mol, mf, _ = carbon_d2h
    key = lane_keys(CAMPAIGN["C"])[0]
    assert (key.window, key.sector) == ("2p@1s2s", (2, 0))
    p = run_lane_cell(key, mol, mf, CAMPAIGN["C"].term, workdir=tmp_path)
    assert p["cell"] == key.address and p["term"] == [1, "g"]
    assert len(p["energies"]) == 3 and p["ptp"] < 1e-9
    for v in p["l2"]:
        assert abs(v - 2.0) < 1e-8          # 3P: L(L+1) = 2
    for v in p["s2"]:
        assert abs(v - 2.0) < 1e-8
    # the state residual: the closed window's term components are
    # symmetry eigenstates -- the residual sits at the apply floor
    # (observed on the frozen D2h twin: 1.1e-16 / 2.2e-16)
    assert 0.0 <= p["res_h_lz"] < 1e-8
    assert 0.0 <= p["res_h_l2"] < 1e-8
    # the store round-trips (and a fresh run is a cache hit)
    q = load_lane_cell(key, workdir=tmp_path)
    assert q["total"] == p["total"] and q["energies"] == p["energies"]
    assert "mo_coeff" not in q and "cycles" not in q  # CI: no optimizer
    assert run_lane_cell(key, mol, mf, CAMPAIGN["C"].term,
                         workdir=tmp_path) == p
    # refusals: a free-route key, a sector disagreeing with the core,
    # a core shell the basis lacks
    with pytest.raises(ValueError, match="lane cell"):
        run_lane_cell(CellKey("C", "2p", (2, 0)), mol, mf, (1, "g"),
                      workdir=tmp_path)
    with pytest.raises(ValueError, match="active electrons"):
        run_lane_cell(CellKey("C", "2p@1s2s", (3, 1), substrate="pur"),
                      mol, mf, (1, "g"), workdir=tmp_path)
    with pytest.raises(ValueError, match="no 5s band"):
        run_lane_cell(CellKey("C", "2p@1s5s", (2, 0), substrate="pur"),
                      mol, mf, (1, "g"), workdir=tmp_path)


@slow
def test_run_lane_cell_multicolumn_fold(carbon_d2h, carbon, tmp_path):
    """The multi-column core end to end: the SAME lane leg with the
    shell layer declaring the 1s2s chemistry over three s bands (two
    pairs over a three-column span) still hosts its term EXACTLY (the
    frozen core is the reference's OWN core restricted to the span --
    the term certificates are core-blind), while the core axes read
    the CHEMISTRY: the payload's pair count is 2 (never the span
    width 3) with the in-span restriction certificate.  A declared
    band shared by window and core refuses at the ROW level (the
    layer's own ambiguity, not the labels')."""
    cfg, mol, mf, _ = carbon_d2h
    key = lane_keys(CAMPAIGN["C"])[0]
    layer = {"2s": ("2s", "3s")}
    p = run_lane_cell(key, mol, mf, CAMPAIGN["C"].term, shell_bands=layer,
                      workdir=tmp_path)
    assert p["core_pairs"] == 2
    assert p["core_in_span_min"] > 1 - 1e-4   # the span hosts the core
    # (observed 0.999994: the core's orthogonality tails reach the
    # valence bands' rows past even the full-family span)
    assert p["ptp"] < 1e-9
    for v in p["l2"]:
        assert abs(v - 2.0) < 1e-8              # 3P: L(L+1) = 2 exact
    for v in p["s2"]:
        assert abs(v - 2.0) < 1e-8
    assert 0.0 <= p["res_h_lz"] < 1e-8
    # the 1:1 leg (no layer): same term, its own span (the energies
    # differ -- DIFFERENT core spans)
    q = run_lane_cell(key, mol, mf, CAMPAIGN["C"].term, workdir=tmp_path,
                      fresh=True)
    assert q["core_pairs"] == 2
    assert q["total"] != p["total"]
    # the partition is CHEMICAL: a shell in both sides refuses...
    with pytest.raises(ValueError, match="sit in both"):
        run_lane_cell(CellKey("C", "2p2s@1s2s", (2, 0), substrate="pur"),
                      mol, mf, (1, "g"), workdir=tmp_path)
    # ...while shared BANDS are physical (the core's 2s realized over
    # the window's 3s band): the projection owns the shared rows
    # -- the leg RUNS, its own model core
    p_shared = run_lane_cell(
        CellKey("C", "2p3s@1s2s", (2, 0), substrate="pur"), mol, mf,
        (1, "g"), shell_bands={"2s": ("2s", "3s")}, workdir=tmp_path)
    assert p_shared["core_pairs"] == 2
    assert p_shared["ptp"] < 1e-9
    for v in p_shared["l2"]:
        assert abs(v - 2.0) < 1e-8
    # a band shared by TWO FOLD SHELLS duplicates its rows -- the span
    # is the UNION (a duplicated row would make the Gram exactly
    # singular; pinned after the Sc-unc lane hit it live)
    p_union = run_lane_cell(
        CellKey("C", "2p@1s2s", (2, 0), substrate="pur"), mol, mf,
        (1, "g"), shell_bands={"1s": ("1s", "2s"), "2s": ("2s",)},
        workdir=tmp_path)
    assert p_union["core_pairs"] == 2
    assert p_union["ptp"] < 1e-9


@slow
def test_no_arm_refuses_the_m_unequal_source_density():
    """The NO arm's construction certificate: a source CI whose ground
    block averages a SPLIT multiplet (the canonical-core field on the
    anisotropic substrate -- the eigenstates carry inter-term
    admixture, their average is m-scalar only to that order) is
    REFUSED at the closure meter, never silently stored -- the exact
    case that fired live on the uncontracted twin (measured: block 5
    roots, groups [1,5], closures 2.7e-2/4.0e-2 against the m-complete
    class at <= 4e-9)."""
    import numpy as np
    from pyscf import scf

    from atomterm.cases import campaign_spec
    from atomterm.natural import no_group_windows
    from atomterm.windows import (build_full_mo, occupied_core_idx,
                                  window_distilled)

    spec = campaign_spec("Sc-unc")
    mol, rec = scf.chkfile.load_scf("elements/Sc-unc/rungs/rohf.chk")
    mf = scf.RHF(mol)
    mf.__dict__.update(rec)
    c_occ = mf.mo_coeff[:, np.asarray(mf.mo_occ) > 0]
    win, _ = window_distilled(mol, c_occ, mf.get_ovlp(), ["3d", "4s"],
                              spec.shell_bands, {"3d": 1, "4s": 1})
    cols = win.columns if hasattr(win, "columns") else win
    # the source on the CANONICAL core (c_core=None) -- the split
    # substrate; the run path cores the purified span instead
    mo = build_full_mo(mf, occupied_core_idx(mf, 9), cols)
    with pytest.raises(ValueError, match="not rotation-closed"):
        no_group_windows(mol, mf, 9, mo[:, 9:9 + cols.shape[1]], 2, 1)


def test_unc_lanes_load_exact_from_the_store():
    """The validation lanes re-verified from the COMMITTED stores (the
    heavy-pin pattern: fresh legs run in the campaign; the verdicts
    load and re-check their meters): Sc-unc hosts 2D and Ti-unc 3F
    EXACTLY through the multi-column cores -- every component at its
    L(L+1) and S(S+1), degenerate to the ptp floor, the core pair
    count chemical (9), the Roothaan converged, and the dense
    adjudication's corrections recorded (the pinned solver's
    near-degenerate skips -- the branch the live runs exercised)."""
    from atomterm.cells import load_lane_cell

    # the distilled twins' sectors are tiny -- the dense adjudication
    # eigen-solves every block and the pinned solver agrees to the
    # machine floor (observed: 8.9e-16 / 2.2e-15, no corrections);
    # U14e's 17-component lane is the live-skip class (the recorded
    # correction's deviation, observed 3.0e-14 ptp era: > 1e-3)
    cases = (
        ("Sc-unc", "sc-unc_3d4s@ar_2-1_casci_ao_-_pur_avg", 6.0, 0.75, 9,
         "exact"),
        ("Ti-unc", "ti-unc_3d4s@ar_3-1_casci_ao_-_pur_avg", 12.0, 2.0, 9,
         "exact"),
        ("U14e", "u14e_5f6d@6s6p7s_4-0_casci_ao_-_pur_avg", 72.0, 6.0, 5,
         "skip"),
    )
    for key_token, addr, l2, s2, pairs, cls in cases:
        p = load_lane_cell(CellKey.from_address(addr))
        assert p["core_pairs"] == pairs               # the chemistry
        # the frozen-core certificate of its era: the repaired lane
        # (the SCF core) carries the in-span weight (observed
        # 0.99991: the core's tails reach past the full-family span);
        # the committed confined-era stores carry the Roothaan's
        # commutator
        assert (p.get("core_in_span_min", 0) > 1 - 1e-3
                or p["core_comm"] < 1e-9)
        assert p["ptp"] < 1e-9                         # degenerate
        for v in p["l2"]:
            assert abs(v - l2) < 1e-8
        for v in p["s2"]:
            assert abs(v - s2) < 1e-8
        assert 0.0 <= p["res_h_lz"] < 1e-8
        if cls == "exact":
            assert p["anchor_dev"] < 1e-12
            assert not p["anchor_corrected"]
        else:
            assert p["anchor_corrected"]                # the live skip
            assert p["anchor_dev"] > 1e-3               # honestly recorded


@slow
def test_confined_core_reproduces_the_ion_rhf(carbon):
    """The EXTERNAL anchor of the confined-core construction (the one
    certificate the internal identities cannot give): on the WHOLE AO
    space the mini-Roothaan over n_pairs pairs is the closed-shell
    ion of the same charge in the same basis -- its RHF energy, to
    the Roothaan convergence floor.  The frozen reference supplies
    the field pieces; the ion SCF is the independent construction."""
    from pyscf import gto, scf

    from atomterm.core import closed_shell_core_energy, confined_core_orbitals
    from atomterm.windows import lowdin

    _, mol, mf, _ = carbon
    ion_mol = gto.M(atom="C 0 0 0", basis="ccpvdz", charge=2, verbose=0)
    ion = scf.RHF(ion_mol)
    ion.conv_tol = 1e-12
    ion.kernel()
    assert ion.converged
    c_span = lowdin(np.eye(mol.nao_nr()), mf.get_ovlp())
    c_core, report = confined_core_orbitals(mf, c_span, 2)
    assert report["comm"] < 1e-10
    assert abs(closed_shell_core_energy(mf, c_core)
               - ion.e_tot) < 1e-8          # C2+ closed shell, 2 pairs


@slow
def test_route_axis_twins(carbon, carbon_d2h, tmp_path):
    """The route axis: route_keys spawns every leg's four twins (the
    substrate axis naming the route); the SPAN-IDENTITY certificate --
    the free CI twin (run_cell, the C1 driver) reproduces the pinned
    leg's ground component to solver tolerance (same window, same
    core, different construction); run_cell's core branch refuses the
    pur class (run_lane_cell's), raw/frac, non-ao starts, and wrong
    sectors; the pinned SA runner delivers the two-pass payload with
    its term exact."""
    from dataclasses import replace

    from atomterm.cells import cell_payload, run_cell
    _, mol_c1, mf_c1, _ = carbon
    _, mol_d2h, mf_d2h, _ = carbon_d2h
    cs = CAMPAIGN["C"]
    leg = lane_keys(cs)[0]

    # the derivation: 1 pinned CI + 4 twins per leg, routes by substrate
    rk = route_keys(cs)
    assert len(rk) == 5 and rk[0].address == leg.address
    tails = sorted(k.address[len(leg.address) - len("casci_ao_-_pur_avg"):]
                   for k in rk[1:])
    assert tails == ["casci_ao_-_pol_avg", "casscf_ao_-_pol_avg",
                     "sa-casscf_ao_-_pol_avg", "sa-casscf_ao_-_pur_avg"]

    # the span-identity certificate at CI level: the free driver's
    # ground TOTAL equals the pinned leg's TERM total (the components
    # are degenerate -- same span, same core, different construction;
    # observed on the frozen twins: 7.1e-15)
    pinned = run_lane_cell(leg, mol_d2h, mf_d2h, cs.term, workdir=tmp_path)
    free = cell_payload(
        replace(leg, substrate="pol"),
        run_cell(replace(leg, substrate="pol"), mf_c1, mol_c1,
                 workdir=str(tmp_path)))
    assert abs(free["energies_pyscf"][0] - pinned["total"]) < 1e-8

    # the core branch's refusals
    for kw, msg in (
            (dict(substrate="pur"), "run_lane_cell"),
            (dict(substrate="pol", core="raw"), "raw/frac"),
            (dict(substrate="pol", start="mo"), "ao/pol/avg"),
            (dict(substrate="pol"), "active electrons")):
        k = (replace(leg, window="2p@1s2s", sector=(3, 1), **kw)
             if msg == "active electrons" else replace(leg, **kw))
        with pytest.raises(ValueError, match=msg):
            run_cell(k, mf_c1, mol_c1, workdir=str(tmp_path))

    # the pinned SA runner: the two-pass payload, term exact
    sa = run_lane_cell(replace(leg, method="sa-casscf"), mol_d2h, mf_d2h,
                       cs.term, workdir=tmp_path)
    assert "ensemble" in sa and sa["ptp"] < 1e-8  # observed: 1.3e-14
    # the SA TOTAL is the kernel ensemble (the optimizer rotates
    # same-irrep core/window columns -- the readout's absolutes are of
    # the run's own partition; only the kernel bookkeeping compares)
    assert sa["total"] == sa["ensemble"]
    assert all(abs(v - 2.0) < 1e-7 for v in sa["l2"])   # 3P exact
    assert sa["total"] < min(pinned["energies"]) + 1.0  # relaxation helps
    # the converged window PERSISTS: the stored payload round-trips
    # the full MO matrix (restart + orbital-level route comparisons)
    # and the optimizer's cycle count
    sa_q = load_lane_cell(replace(leg, method="sa-casscf"),
                          workdir=tmp_path)
    assert sa_q["mo_coeff"].shape == (mol_d2h.nao_nr(),) * 2
    assert np.array_equal(sa_q["mo_coeff"], sa["mo_coeff"])
    assert sa_q["cycles"] == sa["cycles"] > 0
    # the convergence precision persists beside the count (the conv
    # column's data; 0.0 is the legitimately dead-still final cycle)
    assert sa_q["conv_dE"] == sa["conv_dE"] >= 0.0


def methods_module():
    """The methods-table builder, loaded once per call site (a script,
    not a package import)."""
    import importlib.util
    from pathlib import Path as _Path

    loader = importlib.util.spec_from_file_location(
        "build_methods_tables",
        _Path(__file__).resolve().parents[2] / "scripts"
        / "build_methods_tables.py")
    mod = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(mod)
    return mod


def _methods_payload(e, roots=1, kind="casci", n_det=90):
    """A minimal run payload for methods-table rendering."""
    return {"energies_pyscf": [e] * roots, "n_det": n_det, "kind": kind,
            "ground_l2": 12.0, "comm_norm_h_lz": 1e-9,
            "comm_norm_h_l2": 1e-8, "win_lz": 1e-15,
            "win_lperp": 2e-15}


@fast
def test_methods_table_is_byte_deterministic():
    """THE LADDER (one table per atom, between the markers): the
    section builds byte-identically from the same store, once."""
    mod = methods_module()
    cs = CAMPAIGN["Ti"]
    keys = grid_keys(cs)
    agree = {k.address: _methods_payload(-848.4) for k in keys}
    text = "\n".join(mod.methods_section(cs, agree, -848.4065))
    assert text == "\n".join(mod.methods_section(cs, agree, -848.4065))
    assert text.count(mod.BEGIN) == 1 and text.count(mod.END) == 1


@fast
def test_methods_anchors_lead_and_families_follow():
    """The ROHF anchor leads; the UHF anchor row follows when the
    store carries the record (DeltaE against the same zero, the meter
    in the note); the first family follows directly."""
    mod = methods_module()
    cs = CAMPAIGN["Ti"]
    keys = grid_keys(cs)
    agree = {k.address: _methods_payload(-848.4) for k in keys}

    def body(cells, uhf=None):
        return "\n".join(mod.methods_section(cs, cells, -848.4065, uhf=uhf))

    rows = [r for r in body(agree).splitlines() if r.startswith("|")]
    assert rows[2].startswith("| ROHF (anchor) |")
    with_uhf = body(agree, uhf={"e_tot": -848.41, "s2": 4.02,
                                "sigma_min": 0.98})
    rows_uhf = [r for r in with_uhf.splitlines() if r.startswith("|")]
    assert rows_uhf[2].startswith("| ROHF (anchor) |")
    cells_uhf = [c.strip() for c in rows_uhf[3].split("|")[1:-1]]
    assert cells_uhf[0] == "UHF (anchor)"
    assert cells_uhf[7] == "−3.5"
    assert "⟨S²⟩=4.0200" in cells_uhf[17]
    # the anchor block set apart by the SECTION rule, the first
    # family follows
    assert all(c.strip() == "===" for c in rows_uhf[4].split("|")[1:-1])
    assert rows_uhf[5].startswith("| CASCI |")


@fast
def test_methods_primary_renders_core_groups_beside_optimizers():
    """EVERY comparison row a row -- the three rendered core variants'
    CI rows beside their optimizer variants, the method bare +
    restriction only, the window cell in slash form; the frac
    reference class and the NO start (its span IS its source's,
    re-gauged) stay store-only; free optimizer rows carry no lane
    marker."""
    mod = methods_module()
    cs = CAMPAIGN["Ti"]
    keys = grid_keys(cs)
    agree = {k.address: _methods_payload(-848.4) for k in keys}
    text = "\n".join(mod.methods_section(cs, agree, -848.4065))
    rows = [r for r in text.splitlines() if r.startswith("|")]
    ci = [r for r in rows
          if r.startswith("| CASCI | free | 3d4s / D |")
          and r.split("|")[6].strip() != "—"]   # data rows only
    assert len(ci) == 3            # raw/avg/pur (the frac class: store-only)
    assert not any("/ NO |" in r for r in rows)  # the NO class: store-only
    assert not any("frac" in r for r in rows)
    # the IDENTICAL per-start order: the D start's CASSCF rows (free,
    # then the masks in declaration order) before its SA rows
    d_opt = [r.split("|")[1].strip() for r in rows
             if "3d4s / D |" in r and (r.startswith("| CASSCF")
                                       or r.startswith("| SA"))]
    # the virtual-mix diagnostic retired from the declaration
    # (2026-11-12): six optimizer rows at the D start
    assert d_opt[:6] == ["CASSCF", "CASSCF, radial", "CASSCF, frozen",
                         "SA-CASSCF",
                         "SA-CASSCF, radial", "SA-CASSCF, frozen"]
    mo_opt = [r.split("|")[1].strip() for r in rows
              if "3d4s / MO |" in r and (r.startswith("| CASSCF")
                                         or r.startswith("| SA"))]
    seq = ["CASSCF", "CASSCF, frozen", "SA-CASSCF", "SA-CASSCF, frozen"]
    assert mo_opt[:4] == seq                    # the pol group's order
    assert mo_opt[4:] == seq                    # pur repeats it identically
    opt = [r for r in rows
           if r.startswith("| CASSCF | free | 3d4s / D |")
           and r.split("|")[6].strip() != "—"]      # data rows only
    assert len(opt) == 2 and not any("·lane" in r for r in opt)  # Ar + Ar·pur substrates


@fast
def test_methods_local_minimum_note_derived_from_the_sa_floor():
    """A free single-root CASSCF above its span's SA ground carries
    the derived note; below-floor and constrained rows never do."""
    mod = methods_module()
    cs = replace(CAMPAIGN["Sc"], render_spans=None)

    keys = [*grid_keys(cs), *route_keys(cs)]
    cells = {}
    for k in keys:
        if k.method == "sa-casscf" and k.mask is None and k.core == "avg":
            cells[k.address] = _methods_payload(-759.7495, kind="sa-casscf")
        elif k.method == "casscf" and k.core == "avg":
            above = k.window == "3d4s@ar" and k.substrate == "pol"
            cells[k.address] = _methods_payload(
                -759.7436 if above else -759.7680,        # above: flagged
                kind="casscf")
    text = "\n".join(mod.methods_section(cs, cells, -759.7363))
    rows = [l for l in text.splitlines() if l.startswith("| CASSCF")]
    flagged = [l for l in rows if "local minimum" in l]
    assert flagged and all("Ar·lane" in l for l in flagged)
    assert all("frozen" not in l and "radial" not in l for l in flagged)
    assert not any("local minimum" in l for l in rows
                   if "frozen" in l or "radial" in l)


@fast
def test_methods_sc_and_sc_unc_line_by_line_parity():
    """Sc's and Sc-unc's tables carry the SAME row structure (the
    method+window sequence, rules included) -- the line-by-line
    review identity, a maintained invariant; Sc's legacy radial rows
    (a construction class since refused) are the ONE exception, the
    twins' frozencore variety (their declared optimizer_masks --
    the swap class lives there) the second: both sides filter to
    the shared skeleton before comparing."""
    mod = methods_module()
    for a, b in (("Sc", "Sc-unc"), ("Ti", "Ti-unc")):
        _assert_pair_parity(mod, a, b)


def _assert_pair_parity(mod, a, b):
    skel = []
    for atom in (a, b):
        spec = campaign_spec(atom)
        data = json.loads(Path(f"elements/{atom}/cells.json").read_text())
        text = "\n".join(mod.methods_section(
            spec, data["cells"], data["rohf_energy"],
            uhf=data.get("uhf")))
        rows = []
        for l in text.splitlines():
            m = re.match(r"^\| (\S[^|]*?) \| \S[^|]*? \| ([^|]+?) \|", l)
            if m and m.group(1) not in ("method", "ROHF (anchor)",
                                        "UHF (anchor)"):
                if set(l.replace("|", "").replace(" ", "")) <= {"-", "="}:
                    rows.append("RULE")
                else:
                    rows.append(f"{m.group(1).strip()}|{m.group(2).strip()}")
        # the twins' declared frozencore variety is theirs alone --
        # the parity reads the shared skeleton
        rows = [r for r in rows if "frozen core" not in r]
        skel.append(rows)
    assert [r for r in skel[0] if "radial" not in r] == skel[1], (
        f"{a} vs {b} skeletons differ")


@fast
def test_methods_meter_and_ptp_columns_render():
    """The closure columns carry the span certificate; ptp computes
    where a root set exists and stays dash on single roots; the conv
    column carries precision/count."""
    mod = methods_module()
    cs = CAMPAIGN["Ti"]
    keys = grid_keys(cs)
    agree = {k.address: _methods_payload(-848.4) for k in keys}

    def body(cells):
        return "\n".join(mod.methods_section(cs, cells, -848.4065))

    text = body(agree)
    rows = [r for r in text.splitlines() if r.startswith("|")]
    opt = [r for r in rows
           if r.startswith("| CASSCF | free | 3d4s / D |")]
    assert opt[0].split("|")[14].strip() == "1e-15"
    assert opt[0].split("|")[15].strip() == "2e-15"
    assert opt[0].split("|")[13].strip() == "—"  # the single root: no set
    sa = next(k for k in keys if k.method == "sa-casscf"
              and k.mask is None and k.start == "ao"
              and k.substrate == "pol")
    traced = dict(agree)
    traced[sa.address] = {**_methods_payload(-848.4, roots=5),
                          "conv_dE": 2.3e-9, "cycles": 5}
    row = [r for r in body(traced).splitlines()
           if r.startswith("| SA-CASSCF | free | 3d4s / D |")][0]
    assert "| 1e-9/5 |" in row and "| 0e+00 |" in row


@fast
def test_methods_cycle_cap_note():
    """The cap note fires by IDENTITY on payloads carrying the cap
    (the store's ``max_cycle`` field -- the run's actual mc cap) and
    keeps the legacy ``== 200`` key for payloads predating the field;
    a mid-cap count never flags."""
    mod = methods_module()
    cs = CAMPAIGN["Ti"]
    keys = grid_keys(cs)
    agree = {k.address: _methods_payload(-848.4) for k in keys}
    sa = next(k for k in keys if k.method == "sa-casscf"
              and k.mask is None and k.start == "ao"
              and k.substrate == "pol")

    def row_of(extra):
        cells = dict(agree)
        cells[sa.address] = {**_methods_payload(-848.4, roots=5),
                             "conv_dE": 2.3e-9, **extra}
        row = [r for r in "\n".join(
                   mod.methods_section(cs, cells, -848.4065)).splitlines()
               if r.startswith("| SA-CASSCF | free | 3d4s / D |")][0]
        return "cycle cap" in row

    assert row_of({"cycles": 5, "max_cycle": 5})       # identity fires
    assert not row_of({"cycles": 5, "max_cycle": 200})  # below the cap
    assert row_of({"cycles": 200})                     # legacy key fires
    assert not row_of({"cycles": 5})                   # legacy, mid-cap


@fast
def test_methods_ci_only_store_renders_pending_placeholders():
    """No optimizer run cells anywhere: the CI rows render and every
    unrun declared cell holds its PLACE as a pending row (the table
    doubles as the run plan)."""
    mod = methods_module()
    cs = CAMPAIGN["Ti"]
    keys = grid_keys(cs)
    text = "\n".join(mod.methods_section(
        cs, {k.address: _methods_payload(-848.4) for k in keys
             if k.method == "casci"}, -848.4065))
    assert "| CASCI | free | 3d4s / D |" in text
    opt = [l for l in text.splitlines()
           if l.startswith(("| CASSCF", "| SA"))]
    assert opt and all(l.endswith("| pending |") for l in opt)


@fast
def test_methods_prior_rows_join_the_span_family_in_comparison_order():
    """The same-CAS comparison: a prior window's CI rows render BESIDE
    the incumbent's at its span, ordered incumbent -> r -> s1 -> s2.
    The equal-CAS family (D, @r1, @s1) is contiguous -- no rule row
    inside; the DOUBLED window (@s2, a different CAS) opens its own
    === section (the CAS sectioning, 2026-11-12) -- in the lattice
    family of a lattice span and in the primary's core groups of the
    primary span (render_spans unset -- the pin proves the FULL
    ordering)."""
    mod = methods_module()
    cs = replace(CAMPAIGN["Sc"], render_spans=None)
    keys = [k for k in [*lattice_keys(cs), *grid_keys(cs),
                        *(k for w in cs.planned_windows
                          for k in window_grid(cs, w))]
            if k.method == "casci"]
    text = "\n".join(mod.methods_section(
        cs, {k.address: _methods_payload(-759.5) for k in keys}, -759.7363))
    lines = text.splitlines()

    def idx(window, core="Ar"):
        # DATA rows only (the dets column carries a number) -- the
        # lane families' free twins share the window/core display
        return next(i for i, l in enumerate(lines)
                    if l.startswith(f"| CASCI | free | {window} / D |")
                    and l.split("|")[5].strip() == core
                    and l.split("|")[6].strip() != "—")

    for span in ("3d", "3d4s"):
        a, b, c, d = (idx(span), idx(f"{span}@r1"),
                      idx(f"{span}@s1"), idx(f"{span}@s2"))
        assert a < b < c < d
        # the equal-CAS family is contiguous
        assert not any(l.startswith("| --- ") for l in lines[a:c + 1])
        # the doubled window opens its own CAS section (===)
        assert any(l.startswith("| === ") for l in lines[c + 1:d + 1])


@fast
def test_methods_table_diet_repeats_do_not_render():
    """The rendering diet: the prior windows' control rows, the
    lattice's per-span control, the frac reference class, and the
    whole NO start (its span IS its source's, re-gauged -- a
    certificate, not a row) never render; the store keeps
    everything.  Pinned with the spec's render_spans review filter
    unset."""
    mod = methods_module()
    cs = replace(CAMPAIGN["Sc"], render_spans=None)
    keys = [k for k in [*lattice_keys(cs), *grid_keys(cs),
                        *(k for w in cs.planned_windows
                          for k in window_grid(cs, w))]
            if k.method == "casci"]
    text = "\n".join(mod.methods_section(
        cs, {k.address: _methods_payload(-759.5) for k in keys}, -759.7363))
    rows = [l for l in text.splitlines() if l.startswith("| CASCI")]
    assert not any("/ NO" in l for l in rows)   # the NO class: store-only
    assert not any("frac" in l for l in rows)
    # prior windows: their D rows in EVERY core group where they ran
    # (here: the raw and avg groups), never their control clones
    prior_rows = [l for l in rows if "@" in l.split("|")[3]]
    assert prior_rows
    assert {l.split("|")[5].strip() for l in prior_rows} >= {"Ar", "Ar·raw"}
    assert all("/ D |" in l for l in prior_rows)
    # the controls render at EVERY span now (the lattice sections
    # carry the same structure as the primary -- the lattice diet
    # clause retired with the structural round)
    mo_spans = {l.split("|")[3].strip().split(" /")[0]
                for l in rows if "/ MO |" in l}
    assert {"3d4s", "3d"} <= mo_spans


@fast
def test_methods_prior_span_without_incumbent_renders_its_own_family():
    """A prior span with no incumbent family (no lattice, not the
    primary) renders as its own comparison family -- beside the
    primary's, separated by the rule row."""
    mod = methods_module()
    cs = CAMPAIGN["Ti"]
    keys = [k for k in (k for w in cs.planned_windows
                        for k in window_grid(cs, w))
            if k.method == "casci" and k.source_tag is not None]
    text = "\n".join(mod.methods_section(
        cs, {k.address: _methods_payload(-848.2) for k in keys}, -848.4065))
    lines = text.splitlines()

    def idx(window):
        return next(i for i, l in enumerate(lines)
                    if l.startswith(f"| CASCI | free | {window} / D |")
                    and "·" not in l.split("|")[5])

    a, b, d = idx("3d@r1"), idx("3d@s1"), idx("3d@s2")
    assert a < b < d
    # the equal-CAS family contiguous; the doubled window opens its
    # own CAS section (the CAS sectioning, 2026-11-12)
    assert not any(l.startswith("| --- ") for l in lines[a:b + 1])
    assert any(l.startswith("| === ") for l in lines[b + 1:d + 1])


@fast
def test_methods_lane_families():
    """The one table carries the lane families: pinned CI, free CI,
    the free CASSCF/SA, and the pinned SA in one family block; the
    tagged cores carry the lane marker; un-run lane families are
    absent; families are separated by the rule row (the spec's
    render_spans review filter unset -- the lane spans are not in
    Sc's review tuple)."""
    mod = methods_module()
    cs = replace(CAMPAIGN["Sc"], render_spans=None)

    def p(e, pinned=False, ci=False):
        d = {"n_det": 204490, "comm_norm_h_lz": float("nan"),
             "comm_norm_h_l2": float("nan"), "res_h_lz": 1e-14,
             "res_h_l2": 2e-14, "win_lz": 1e-15, "win_lperp": 1e-15}
        if pinned or ci:
            d.update({"total": e, "l2": [6.0] * 5, "ptp": 5e-8,
                      "energies": [e] * 5})
        else:
            d["energies_pyscf"] = [e]
            d["ground_l2"] = 6.0
        return d

    routes = route_keys(cs)
    fam = [k for k in routes if k.window == "3d3p4d@semi"]
    cells = {k.address: p(-759.85,
                          pinned=(k.substrate == "pur"
                                  and k.method == "sa-casscf"),
                          ci=(k.method == "casci"))
             for k in fam}
    text = "\n".join(mod.methods_section(cs, cells, -759.7363))
    lines = text.splitlines()
    assert sum(l.startswith("| method |") for l in lines) == 1
    assert lines[3].startswith("| ROHF (anchor) |")
    rules = [l for l in lines
             if set(l.replace("|", "").replace(" ", "")) == {"-"}]
    assert len(rules) >= 2    # the column SEP + the anchor rule
    # the FIRST family leads without a rule of its own
    after = lines[lines.index(rules[1]) + 1]
    assert after.startswith("| CASCI")
    assert sum(l.startswith("| CASCI | pinned | 3d3p4d / D |")
               for l in lines) == 1
    # the review-stage span filter matches the SPAN, never the raw
    # label -- a tagged lane window renders when its span is listed
    cs_f = replace(cs, render_spans=("3d3p4d",))
    text_f = "\n".join(mod.methods_section(cs_f, cells, -759.7363))
    assert sum(l.startswith("| CASCI | pinned | 3d3p4d / D |")
               for l in text_f.splitlines()) == 1
    assert sum(l.startswith("| CASCI | free | 3d3p4d / D |")
               for l in lines) == 1
    assert sum(l.startswith("| CASSCF | free | 3d3p4d / D |")
               for l in lines) == 1
    assert sum(l.startswith("| SA-CASSCF | pinned | 3d3p4d")
               for l in lines) == 1
    assert "·lane" in [l for l in lines
                       if l.startswith("| SA-CASSCF | pinned | 3d3p4d")
                       ][0]
    assert "3d3p / AO" not in text       # un-run families absent



@fast
def test_swap_note_wording_and_band_counter_demotion():
    """The renderer's two display contracts (the review round, this
    wave): a swap_basin payload notes 'swap' -- the user's wording,
    never the longer form; a chemically named window (window_names)
    renders the chemical name and demotes its basis-band counters to
    the note column."""
    mod = methods_module()

    # the swap wording: 'swap', and the longer form nowhere
    cs = replace(CAMPAIGN["Sc"], render_spans=None)
    keys = [*grid_keys(cs), *route_keys(cs)]
    cells = {}
    for k in keys:
        if k.method == "casci" and k.start == "ao" and k.core == "avg":
            cells[k.address] = _methods_payload(-759.7342,
                                                kind="casci")
    swap_key = next(k for k in keys if k.method == "casscf"
                    and k.mask is None and k.start == "ao")
    p = _methods_payload(-759.7593, kind="casscf")
    p["swap_basin"] = True
    p["converged"] = False
    cells[swap_key.address] = p
    text = "\n".join(mod.methods_section(cs, cells, -759.7363))
    assert "swap basin" not in text
    flagged = [l for l in text.splitlines() if "; swap" in l]
    assert flagged and all(l.rstrip("| ").endswith("swap") for l in flagged)

    # the band-counter demotion: the chemical name renders, the
    # counters land in the note column
    cs_u = CAMPAIGN["U14e"]
    ukey = CellKey("U14e", "10f13d", (4, 0))
    ucells = {ukey.address: _methods_payload(-50.35, kind="casci")}
    utext = "\n".join(mod.methods_section(cs_u, ucells, -50.35))
    row = [l for l in utext.splitlines()
           if "5f6d" in l and "pending" not in l][0]
    assert "bands 10f+13d" in row
