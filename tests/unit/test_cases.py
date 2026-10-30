"""Cases: the declared per-atom campaign and the launch seams
(CaseSpec / run_case / run_atom). Unit tier: every test here raises
before any SCF runs; the executed pipeline is pinned by
tests/integration/test_atoms.py."""

import dataclasses

import pytest

from atomterm.cases import CaseSpec, run_atom
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
            "nroots_primary", "planned_windows", "lattice_windows"} <= names
    # immutable: a declaration is not mutated by running it
    with pytest.raises(Exception):
        spec.symbol = "N"


def test_casespec_atom_config_roundtrip():
    spec = _spec(tag="c_trial")
    cfg = spec.atom_config()
    assert isinstance(cfg, AtomConfig)
    assert cfg.spin == 0


def test_run_atom_rejects_parity_violations(tmp_path):
    # spin 2S and the active-electron count must share parity -- the
    # Ms-sector construction depends on it; raises before any SCF
    with pytest.raises(ValueError, match="parity"):
        run_atom(AtomConfig("C", "ccpvdz", spin=1), ["2p"], 2,
                 workdir=str(tmp_path))
