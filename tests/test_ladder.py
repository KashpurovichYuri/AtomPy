"""Ladder tests: per-rung symmetry fingerprints, parameterized over atoms.

Each case runs the full orbital-treatment ladder (ticket 08 on carbon,
ticket 09 on Si) and asserts the MEASURED pattern, with per-atom chaotic
bands -- never assumed exactness:

- raw core: the internal field is real (split + [H, Lz]) but bounded;
- averaged / fractional-occupation cores: exact invariance;
- free CASSCF: symmetry breaking is real but the window stays l-pure
  (fact 8: m-incompleteness, not s/d leakage) and the optimizer stays
  orthonormal (the ticket-09 build_full_mo regression guard);
- SA: restoration to machine precision.

CASSCF gauges are chaotic across identical reruns (fact 10): <L^2> was
observed 0.99..1.92 on carbon -- bands assert "real but not catastrophic",
with the precisely measured values living in the tickets.
"""

from collections import Counter

import numpy as np
import pytest

from atomterm.ladder import run_carbon_ladder, run_ladder
from atomterm.system import AtomConfig

LADDER_CASES = {
    "C": dict(
        cfg=AtomConfig(),
        patterns=["2p"],
        nelec_act=2,
        bands=dict(
            raw_comm=(1e-7, 1e-3),      # observed 3.8e-6 .. 2.0e-5
            split=(1e-6, 1e-3),         # observed 4.4e-5 (avg-removal 2.9e-5)
            avg_dE=(1e-6, 1e-3),        # |E_avg - E_raw|
            frac_dE=(1e-6, 1e-1),       # |E_frac - E_avg|, observed 6.6e-3
            corr_gain=(-0.3, -0.1),     # free CASSCF vs casci_avg: -0.251
            relax_offset=1e-9,          # same window definition: monotone
        ),
    ),
    "Si": dict(
        cfg=AtomConfig("Si", "ccpvdz", spin=2),
        patterns=["3p"],
        nelec_act=2,
        bands=dict(
            raw_comm=(1e-5, 1e-2),      # observed ~6e-4 (bigger core)
            split=(1e-5, 1e-3),         # observed 3.1e-4
            avg_dE=(1e-6, 1e-3),        # observed 3.1e-4
            frac_dE=(1e-6, 1e-1),       # observed 2.4e-3
            corr_gain=(-0.2, -0.05),    # observed -0.127
            # Si's raw 3p window and PySCF's core-orthogonalized window
            # differ physically (~0.3 overlap): the relax rung optimizes
            # the latter, so its total sits ABOVE the raw-shell casci_avg
            # by the window-definition offset (fact 14)
            relax_offset=0.2,
        ),
    ),
}


@pytest.fixture(scope="session")
def ladders(tmp_path_factory):
    out = {}
    for name, case in LADDER_CASES.items():
        workdir = tmp_path_factory.mktemp(f"ladder_{name}")
        if name == "C":
            out[name] = run_carbon_ladder(case["cfg"], workdir=str(workdir))
        else:
            out[name] = run_ladder(
                case["cfg"], case["patterns"], case["nelec_act"],
                workdir=str(workdir), sa_weights=[1 / 3] * 3,
                frac_nelec_open=2, frac_n_open=3)
    return out


@pytest.mark.parametrize("name", list(LADDER_CASES))
def test_casci_fingerprints(name, ladders):
    """Raw core: split + [H, Lz] present and bounded; avg/frac: exact."""
    case, r = LADDER_CASES[name], ladders[name]
    lo, hi = case["bands"]["raw_comm"]
    for tag, comm_bound in [("casci_raw", (lo, hi)),
                            ("casci_avg", (0.0, 1e-8)),
                            ("casci_frac", (0.0, 1e-8))]:
        rung = r[tag]
        assert rung["terms"] == {"3P": 3}, (name, tag)
        assert rung["ortho_err"] < 1e-10, (name, tag)
        assert comm_bound[0] <= rung["comm_norm_h_lz"] < comm_bound[1], (name, tag)
        assert rung["configurations"][0].startswith("p^2.00"), (name, tag)
    # [H, L^2] is trivially zero in a complete single-irrep sector
    assert r["casci_raw"]["comm_norm_h_l2"] < 1e-12
    # 3P split: {M_L=0} off the {M_L=+-1} pair
    roots = np.sort(r["casci_raw"]["energies_pyscf"])
    lo, hi = case["bands"]["split"]
    assert lo < roots[2] - roots[0] < hi, name
    assert abs(roots[2] - roots[1]) < (roots[2] - roots[0]) / 10, name
    # averaging shifts the spectrum by exactly the splitting scale
    lo, hi = case["bands"]["avg_dE"]
    assert lo < abs(r["casci_avg"]["e_ground"] - r["casci_raw"]["e_ground"]) < hi
    # frac-occ core: invariant by construction, physically distinct reference
    lo, hi = case["bands"]["frac_dE"]
    assert lo < abs(r["casci_frac"]["e_ground"] - r["casci_avg"]["e_ground"]) < hi
    # discarded core anisotropy: measured small (density norm)
    assert 1e-4 < r["casci_raw"]["core_anisotropy"] < 1e-1, name


@pytest.mark.parametrize("name", list(LADDER_CASES))
def test_core_relaxation(name, ladders):
    """Frozen window + relaxed core: real physics, zero escape, exact symmetry."""
    case, r = LADDER_CASES[name], ladders[name]
    relax = r["casscf_relax"]
    assert relax["comm_norm_h_lz"] < 1e-3, name
    assert relax["ortho_err"] < 1e-10, name
    assert relax["ground_l2"] == pytest.approx(2.0, abs=1e-4), name
    assert relax["terms"] == {"3P": 3}, name
    assert relax["window"]["p"] > 0.95, name
    assert relax["configurations"][0].startswith("p^2.00"), name
    # monotone within one window definition; on Si the window-definition
    # offset (fact 14) dominates, bounded per case
    assert relax["e_ground"] <= \
        r["casci_avg"]["e_ground"] + case["bands"]["relax_offset"], name


@pytest.mark.parametrize("name", list(LADDER_CASES))
def test_free_casscf_and_virtmix(name, ladders):
    """Free CASSCF breaks symmetry WITHOUT leaving the shell (fact 8:
    m-incompleteness, l-pure window) and stays orthonormal; virtmix lands
    in the same basin as free CASSCF (same-basin assertion is the
    regression guard for the ticket-09 bug class)."""
    case, r = LADDER_CASES[name], ladders[name]
    lo, hi = case["bands"]["corr_gain"]
    e_ref = None
    for tag in ("casscf_shell", "casscf_virtmix", "casscf_pyscf"):
        rung = r[tag]
        assert rung["ground_term"] == "3P", (name, tag)
        assert rung["ortho_err"] < 1e-10, (name, tag)
        assert 1e-4 < rung["comm_norm_h_lz"] < 1.0, (name, tag)
        assert 0.9 < rung["ground_l2"] < 2.0, (name, tag)
        assert rung["window"]["p"] > 0.9, (name, tag)
        d_e = rung["e_ground"] - r["casci_avg"]["e_ground"]
        assert lo < d_e < hi, (name, tag)
        if e_ref is None:
            e_ref = rung["e_ground"]
        assert abs(rung["e_ground"] - e_ref) < 1e-2, (name, tag)


@pytest.mark.parametrize("name", list(LADDER_CASES))
def test_sa_restoration(name, ladders):
    """SA over the term's M_L components restores L to machine precision."""
    _, r = LADDER_CASES[name], ladders[name]
    for tag in ("sa", "sa_pyscf"):
        if tag not in r:
            continue
        sa = r[tag]
        assert sa["ortho_err"] < 1e-10, (name, tag)
        assert sa["comm_norm_h_lz"] < 1e-7, (name, tag)
        assert sa["ground_l2"] > 1.9999, (name, tag)
        assert sa["terms"] == {"3P": 3}, (name, tag)
