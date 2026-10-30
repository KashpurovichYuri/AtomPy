"""Tests of atomterm.stages: the one-pass fold-and-analyze facts module.

The physics content is pinned end-to-end by test_treatments/test_atoms;
these pin the MODULE's contract on the frozen carbon reference (no SCF,
seconds-fast): what stage_facts computes, the payload shape every
persisted rung shares, and the divergence-guard property.
"""

import json

import numpy as np
import pytest

from atomterm.casscf import run_casci
from atomterm.stages import ORTHO_TOL, stage_facts
from atomterm.windows import build_full_mo

pytestmark = pytest.mark.fast


def test_stage_facts_payload_and_invariants(carbon, tmp_path):
    """stage_facts owns the fold: the Ms=0 p^2 sector reproduces the
    textbook 3P/1D/1S table with machine-precision closure, vanishing
    N-deviation, and a payload carrying every persisted rung key."""
    _, mol, mf, act_idx = carbon
    ncore = (mol.nelectron - 2) // 2
    mo = build_full_mo(mf, list(range(ncore)), act_idx)
    cas = run_casci(mf, 3, (1, 1), mo, str(tmp_path / "st.chk"), nroots=9)
    c_act = np.eye(mol.nao_nr())[:, act_idx]
    facts = stage_facts(mol, mf, cas, c_act, ncore, 1, 1, tag="st")

    assert facts.n_dev < 1e-9
    assert facts.ortho_err < 1e-10 and not facts.diverged
    assert facts.terms == {"3P": 3, "1D": 5, "1S": 1}
    assert facts.analysis.comm_norm_h_lz < 1e-12
    assert facts.ground_term == "3P"
    assert facts.window["p"] == 1.0

    payload = facts.payload("casci")
    for key in ("tag", "kind", "n_det", "diverged", "ortho_err", "e_ground",
                "energies_pyscf", "max_de_vs_pyscf", "ground_term",
                "ground_l2", "terms", "comm_norm_h_lz", "comm_norm_h_l2",
                "window", "n_dev", "configurations"):
        assert key in payload, key
    assert payload["e_ground"] == facts.e_ground
    assert payload["diverged"] == (facts.ortho_err > ORTHO_TOL)
    # the payload is json-serializable (the rungs persist it verbatim)
    json.dumps(payload)


def test_divergence_guard_fires_on_broken_orthonormality(carbon, tmp_path):
    """Negative control for the ortho_err guard (ticket 14): a deliberately
    non-S-orthonormal MO matrix (pre-fix-style) must be flagged DIVERGED --
    the energies of such a rung are meaningless and the flag says so."""
    _, mol, mf, act_idx = carbon
    ncore = (mol.nelectron - 2) // 2
    mo = build_full_mo(mf, list(range(ncore)), act_idx)
    mo[:, ncore] = 1.05 * mo[:, ncore]  # break S-orthonormality (~0.1)
    cas = run_casci(mf, 3, (2, 0), mo, str(tmp_path / "bad.chk"), nroots=3,
                    name="bad")
    facts = stage_facts(mol, mf, cas, mo[:, ncore:ncore + 3], ncore, 2, 0,
                        tag="bad")
    assert facts.ortho_err > 1e-6
    assert facts.diverged
    assert facts.payload("casci")["diverged"] is True


def test_stage_facts_anisotropy_ref_roundtrip(carbon, tmp_path):
    """The anisotropy_ref knob: raw vs averaged reference cores of the SAME
    stage differ exactly by the discarded core anisotropy (the CASCI-rung
    detector); without the knob the field is absent."""
    _, mol, mf, act_idx = carbon
    ncore = (mol.nelectron - 2) // 2
    mo = build_full_mo(mf, list(range(ncore)), act_idx)
    cas = run_casci(mf, 3, (2, 0), mo, str(tmp_path / "st2.chk"), nroots=3)
    c_act = np.eye(mol.nao_nr())[:, act_idx]
    plain = stage_facts(mol, mf, cas, c_act, ncore, 2, 0, tag="plain")
    aniso = stage_facts(mol, mf, cas, c_act, ncore, 2, 0, tag="aniso",
                        anisotropy_ref=(mf, None))
    assert plain.core_anisotropy is None
    assert "core_anisotropy" not in plain.payload("casci")
    assert 1e-5 < aniso.core_anisotropy < 1e-2  # carbon's small raw field
    assert aniso.core_anisotropy == aniso.payload("casci")["core_anisotropy"]
