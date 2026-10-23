"""Cases: the declared per-atom campaign and the new-atom launch
(CaseSpec / run_case). Fast tier -- no SCF anywhere."""

import dataclasses

import pytest

from atomterm.cases import CaseSpec, run_case
from atomterm.system import AtomConfig


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
            "want"} <= names
    # immutable: a declaration is not mutated by running it
    with pytest.raises(Exception):
        spec.symbol = "N"


def test_run_case_demands_term_tables():
    # T1/T2 are the contract; a case without them is not a case
    with pytest.raises(AssertionError):
        run_case(_spec(want=frozenset({"casscf", "sa"})))


def test_casespec_atom_config_roundtrip():
    spec = _spec(tag="c_trial")
    cfg = spec.atom_config()
    assert isinstance(cfg, AtomConfig)
    assert cfg.spin == 0
