"""Cases: the declared per-atom campaign and the launch seams
(CaseSpec / run_case / run_atom). Unit tier: every test here raises
before any SCF runs; the executed pipeline is pinned by
tests/integration/test_atoms.py."""

import dataclasses

import pytest

from atomterm.cases import (CAMPAIGN, CaseSpec, Scenario, run_atom,
                            scenarios_for)
from atomterm.system import AtomConfig

pytestmark = pytest.mark.fast


def _spec(**kw):
    base = dict(symbol="C", basis={}, nelec_act=2, ncore=1,
                shell_patterns=["2p"], spin=0)
    base.update(kw)
    return CaseSpec(**base)


def test_casespec_is_the_intent_manifest():
    spec = _spec()
    # the declaration carries everything the campaign needs
    names = {f.name for f in dataclasses.fields(CaseSpec)}
    assert {"symbol", "basis", "ecp", "spin", "nelec_act", "ncore",
            "shell_patterns", "ref", "sa_weights", "run_casscf",
            "nroots_primary", "planned_windows", "lattice_windows",
            "lattice_sector", "scan", "term", "lane_windows",
            "scenarios"} <= names
    # immutable: a declaration is not mutated by running it
    with pytest.raises(Exception):
        spec.symbol = "N"


def test_casespec_atom_config_roundtrip():
    spec = _spec(tag="c_trial")
    cfg = spec.atom_config()
    assert isinstance(cfg, AtomConfig)
    assert cfg.spin == 0


def test_lattice_declaration_needs_its_sector():
    # the lattice sector says WHAT the series studies (free-ion vs
    # valence); windows without it refuse at the declaration
    with pytest.raises(ValueError, match="lattice sector"):
        _spec(lattice_windows=("2p",))
    _spec(lattice_windows=("2p",), lattice_sector=(2, 0))


def test_run_atom_rejects_parity_violations(tmp_path):
    # spin 2S and the active-electron count must share parity -- the
    # Ms-sector construction depends on it; raises before any SCF
    with pytest.raises(ValueError, match="parity"):
        run_atom(AtomConfig("C", "ccpvdz", spin=1), ["2p"], 2,
                 workdir=str(tmp_path))


def test_scenario_grammar_refusals():
    # the declaration's own grammar: empty name/cells/reason all
    # refuse at construction
    ok = dict(name="s", cells=("ti_3d4d_",), reason="dated claim")
    Scenario(**ok)
    for kw in (dict(name=""), dict(cells=()), dict(reason="")):
        with pytest.raises(ValueError):
            Scenario(**{**ok, **kw})


def test_scenario_gates():
    # the run gate's vocabulary: prefix gating, enabled = the whole
    # scenario's gate
    on = Scenario(name="s", cells=("ti_3d4d_",), reason="r")
    off = Scenario(name="s", cells=("ti_3d4d_",), reason="r",
                   enabled=False)
    assert on.gates("ti_3d4d_3-1_casscf_ao_-_pol_avg")
    assert not on.gates("ti_3d4s_3-1_casscf_ao_-_pol_avg")
    assert not off.gates("ti_3d4d_3-1_casscf_ao_-_pol_avg")


def test_scenarios_for_merges_by_name(monkeypatch):
    # campaign defaults + the row's own, merged by name -- the row
    # wins; new names append in declaration order
    import atomterm.cases as cases

    monkeypatch.setattr(
        cases, "DEFAULT_SCENARIOS",
        {"d": Scenario(name="d", cells=("c_2p_",), reason="default")})
    row = _spec(scenarios=(
        Scenario(name="d", cells=("c_2p_",), reason="row wins"),
        Scenario(name="own", cells=("c_",), reason="own")))
    got = scenarios_for(row)
    assert [s.name for s in got] == ["d", "own"]
    assert got[0].reason == "row wins"


def test_casespec_validates_lane_windows():
    # lane windows validate AT CONSTRUCTION: the label carries a
    # resolvable core tag that stays out of its window span, and the
    # legs need the term target declared
    ti = dict(symbol="Ti", basis="ccpvdz", nelec_act=4, ncore=9,
              shell_patterns=["3d", "4s"], spin=2, term=(3, "g"))
    CaseSpec(**ti, lane_windows=("3d3p@semi",))          # the semicore leg
    CaseSpec(**ti, lane_windows=("3s3p3d@1s2s2p4s",))    # the M-shell core
    for lw in (("3d3p",),        # no core tag
               ("3d3p@nope",),   # unresolvable tag
               ("3p@ar",)):      # 3p folded AND active
        with pytest.raises(ValueError):
            CaseSpec(**ti, lane_windows=lw)
    with pytest.raises(ValueError):   # legs without the term target
        CaseSpec(**{k: v for k, v in ti.items() if k != "term"},
                 lane_windows=("3d3p@semi",))


def test_shell_radials_must_cover_every_declared_window():
    """The compression grammar: every window shell the case names
    (primary, lattice, lane spans, planned) carries its k -- a mixed
    distilled/raw window is refused at the declaration, not at the
    first run."""
    ti = dict(symbol="Ti", basis="ccpvdz", nelec_act=4, ncore=9,
              shell_patterns=["3d", "4s"], spin=2, term=(3, "g"))
    CaseSpec(**ti, shell_radials={"3d": 1, "4s": 1},
             lattice_windows=("3d",), lattice_sector=(2, 0),
             lane_windows=("3d4s@ar",))
    for rad, kw in (
            ({"3d": 1}, {}),                                   # primary
            ({"3d": 1, "4s": 1}, {"lattice_windows": ("3d4d",),
                                  "lattice_sector": (2, 0)}),   # lattice
            ({"3d": 1, "4s": 1}, {"planned_windows": ("3d4p",)})):
        with pytest.raises(ValueError, match="shell_radials"):
            CaseSpec(**ti, shell_radials=rad,
                     lane_windows=("3d4s@ar",), **kw)


def test_twins_declare_the_seeded_honest_series():
    """The Stage-2 declarations: the twins' FULL series speaks the
    honest l-vocabulary (the shell-counter labels would select 1-5 AO
    slivers -- the disease's own declaration bug); the heads declare
    the SOURCE basis's own vocabulary -- the whole-l closed rank on a
    contracted source ties occupied valence (4s ~1 on l0) against
    virtual polarization (4p ~1 on l1) and steals virtuals into the
    closed block (measured 2026-11-14: the Ti/6-31G l-hook source
    diverged at E = -838.91 with eps +0.09 virtuals paired); the
    members declare their own configurations' source vocabulary (the
    4s closed only where the configuration closes it), the mixed
    members pin l0 in their open spans (the 4s open column's
    character lives in l0), and every member carries the
    diis+newton recipe (Newton alone crosses into the deeper sibling
    basin -- the d3 fall to d2s1 at -128 mEh, reproduced 2026-11-16;
    plain DIIS oscillates unconverged)."""
    for key in ("Sc-unc", "Ti-unc"):
        cs = CAMPAIGN[key]
        assert cs.series_seed == "6-31g"
        # the series covers the valence configurations: the head,
        # the mixed d(n-1)s1, and the pure dn member
        assert len(cs.series) == 3
        head, mixed, pure = cs.series
        assert mixed.tag.endswith("s1")
        assert pure.tag == f"d{cs.nelec_act}"
        for m in cs.series:
            assert m.closed == ("l0", "l1")
            assert m.e_ref is not None
        # the mixed member: the 4s open column's character lives in l0
        assert mixed.open == ("l0", "l2")
        assert pure.open == ("l2",)
        # the members' recipe (head never re-runs; its accel is dead)
        for m in cs.series[1:]:
            assert m.accel == "diis+newton"
        assert cs.ref.seed_basis == "6-31g"
        assert cs.ref.source_closed == ("1s", "2s", "2p", "3s", "3p",
                                        "4s")
        assert cs.ref.source_open == ("3d",)
        assert cs.ref.occ_e_ref is not None
        # the members' sources close the core only (their
        # configurations leave the 4s open or empty)
        for m in cs.series[1:]:
            assert m.source_closed == ("1s", "2s", "2p", "3s", "3p")
            assert m.source_open in (("3d", "4s"), ("3d",))


def test_u14e_declares_the_series_revival():
    """The ECP atom's prior axis (tickets 39 + 40, 2026-11-17): the
    THREE prior windows -- the @s half REVIVED by the self-seeded MOM
    (the census landed both members; the diagnosis re-attributed the
    weight hook's refusal to its character-only premium, not the
    landscape), the series declares the three measured basins
    (head f3d1 + f4 + f2d2, the members seed="self" with the per-span
    counts the flat ranking cannot express), the prior realization is
    the certified band families, and the prior windows' optimizer
    cells are gated (the comparison is CASCI)."""
    cs = CAMPAIGN["U14e"]
    assert cs.prior_bands == {"5f": ("9f", "10f", "11f"),
                              "6d": ("12d", "13d", "14d")}
    head, f4, f2d2 = cs.series
    assert (head.tag, head.spin, head.e_ref) == ("f3d1", 4, -50.38185233)
    assert head.seed is None and head.open_counts is None  # never re-run
    assert (f4.tag, f4.open, f4.e_ref) == ("f4", ("l3",), -50.31571501)
    assert (f4.seed, f4.open_counts) == ("self", (4,))
    assert (f2d2.tag, f2d2.open, f2d2.e_ref) == ("f2d2", ("l3", "l2"),
                                                 -50.34788201)
    assert (f2d2.seed, f2d2.open_counts) == ("self", (2, 2))
    # all members spin 4 (Hund; the head is the case's own reference)
    assert all(s.spin == cs.spin for s in cs.series)
    tagged = [w for w in cs.planned_windows if "@" in w]
    assert tagged == ["5f6d@r1", "5f6d@s1", "5f6d@s2"]
    assert "prior-anchor-first" in {s.name for s in cs.scenarios}
