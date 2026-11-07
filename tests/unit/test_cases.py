"""Cases: the declared per-atom campaign and the launch seams
(CaseSpec / run_case / run_atom). Unit tier: every test here raises
before any SCF runs; the executed pipeline is pinned by
tests/integration/test_atoms.py."""

import dataclasses

import pytest

from atomterm.cases import CaseSpec, Scenario, run_atom, scenarios_for
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
    # resolvable fold tag that stays out of its window span, and the
    # legs need the term target declared
    ti = dict(symbol="Ti", basis="ccpvdz", nelec_act=4, ncore=9,
              shell_patterns=["3d", "4s"], spin=2, term=(3, "g"))
    CaseSpec(**ti, lane_windows=("3d3p@semi",))          # the semicore leg
    CaseSpec(**ti, lane_windows=("3s3p3d@1s2s2p4s",))    # the M-shell fold
    for lw in (("3d3p",),        # no fold tag
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
    compressed/raw window is refused at the declaration, not at the
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
