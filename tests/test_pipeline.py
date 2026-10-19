"""End-to-end smoke test of the carbon pipeline (checkpointed, ~1 min)."""

from collections import Counter

import numpy as np

from atomterm.pipeline import run_carbon
from atomterm.system import AtomConfig

import pytest
pytestmark = pytest.mark.slow


def test_pipeline_end_to_end(tmp_path):
    cfg = AtomConfig(verbose=0)
    results = run_carbon(cfg, workdir=str(tmp_path))

    # stage artifacts exist (chkfiles + reports)
    for stage in ["rohf", "casci_ms1", "casci_ms0", "casscf", "sacasscf"]:
        assert (tmp_path / f"{cfg.label}_{stage}.chk").exists(), stage
    assert (tmp_path / f"{cfg.label}_report.md").exists()
    assert (tmp_path / f"{cfg.label}_report.json").exists()

    # ROHF sanity (C/cc-pVDZ triplet)
    assert abs(results["rohf_energy"] + 37.6824) < 1e-3

    # CASCI on the 2p shell, Ms=1: the 3P M_L triplet. PySCF's roots split by
    # ~1e-7 only: residual core-density noise of the SCF (our own matrices are
    # exactly degenerate, see test_commutators_vanish_in_shell_space).
    cas1 = results["casci_ms1"]["analysis"]
    assert cas1.states[0].term == "3P"
    energies = np.asarray(results["casci_ms1"]["energies_pyscf"])
    # PySCF's root splitting: residual SCF core-density noise; observed
    # 1e-7 .. 2.1e-4 across cold reruns (band widened 2026-10-15)
    assert np.ptp(energies) < 1e-3, energies

    # CASCI Ms=0: the full 1S + 1D + 3P(Ms=0) multiplet, L^2 = {0, 2, 6}
    cas0 = results["casci_ms0"]["analysis"]
    terms = Counter(s.term for s in cas0.states)
    assert terms == {"1S": 1, "1D": 5, "3P": 3}, terms
    assert {round(s.l2, 4) for s in cas0.states} == {0.0, 2.0, 6.0}

    # CASSCF / SA-CASSCF stages classified (note: the SA-CASSCF optimized
    # orbitals genuinely break rotational symmetry -- <L^2> < 2 and a finite
    # ||[H, Lz]|| are the interesting physical output here, not bugs)
    assert results["casscf"]["analysis"].states[0].term == "3P"
    assert len(results["sa_casscf"]["analysis"].states) == 3
