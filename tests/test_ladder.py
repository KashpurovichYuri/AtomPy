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

The Sc tests (ticket 11) pin the TRANSITION-METAL departures from that
story: a TM core polarizes at the percent level, the free-CASSCF window
escapes into p (an 2S-like basin, still orthonormal/variational), and the
shell-start SA restores <L^2> only -- the window's s member swaps with a
core p, so [H, Lz] stays large (m-incomplete window), while the
canonical-start SA collapses outright (start dependence).
"""

from collections import Counter

import numpy as np
import pytest

from conftest import SC_EXACT_TERMS

from atomterm.ladder import run_ladder
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
            relax_window_offset=1e-9,   # same window definition: monotone
            # tightened bands restored post-review (ticket 15 step 1): the
            # ticket-10 merge had unified them to the loosest atom
            exact_comm=1e-9,            # avg/frac commutators, carbon-exact
            relax_comm=1e-9,
            relax_l2_abs=1e-6,
            relax_purity=0.999,
            sa_comm=1e-9,
            sa_l2=1.99999,
            sa_nroots=True,             # assert len(energies_pyscf) == 3
            basin_tol=1e-4,
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
            relax_window_offset=0.2,
            exact_comm=1e-8,            # avg/frac commutators
            relax_comm=1e-3,
            relax_l2_abs=1e-4,
            relax_purity=0.95,
            sa_comm=1e-7,
            sa_l2=1.9999,
            sa_nroots=False,
            basin_tol=1e-3,
        ),
    ),
}


@pytest.fixture(scope="session")
def ladders(tmp_path_factory):
    out = {}
    for name, case in LADDER_CASES.items():
        workdir = tmp_path_factory.mktemp(f"ladder_{name}")
        out[name] = run_ladder(
            case["cfg"], case["patterns"], case["nelec_act"],
            workdir=str(workdir), sa_weights=[1 / 3] * 3,
            frac_nelec_open=2, frac_n_open=3)
    return out


@pytest.mark.parametrize("name", list(LADDER_CASES))
def test_casci_fingerprints(name, ladders):
    """Raw core: split + [H, Lz] present and bounded; avg/frac: exact."""
    case, r = LADDER_CASES[name], ladders[name]
    b = case["bands"]
    lo, hi = b["raw_comm"]
    for tag, comm_bound in [("casci_raw", (lo, hi)),
                            ("casci_avg", (0.0, b["exact_comm"])),
                            ("casci_frac", (0.0, b["exact_comm"]))]:
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
    b = case["bands"]
    relax = r["casscf_relax"]
    assert relax["comm_norm_h_lz"] < b["relax_comm"], name
    assert relax["ortho_err"] < 1e-10, name
    assert relax["ground_l2"] == pytest.approx(2.0, abs=b["relax_l2_abs"]), name
    assert relax["terms"] == {"3P": 3}, name
    assert relax["window"]["p"] > b["relax_purity"], name
    assert relax["configurations"][0].startswith("p^2.00"), name
    # monotone within one window definition; on Si the window-definition
    # offset (fact 14) dominates, bounded per case
    assert relax["e_ground"] <= \
        r["casci_avg"]["e_ground"] + b["relax_window_offset"], name


def _assert_same_basin(rungs: dict, tags: tuple, r, ref: str, tol: float):
    """All ``tags`` land in one optimizer basin (|E - E_ref| < tol each)."""
    e_ref = rungs[ref]["e_ground"]
    for tag in tags:
        assert abs(rungs[tag]["e_ground"] - e_ref) < tol, tag


@pytest.mark.parametrize("name", list(LADDER_CASES))
def test_free_casscf_and_virtmix(name, ladders):
    """Free CASSCF breaks symmetry WITHOUT leaving the shell (fact 8:
    m-incompleteness, l-pure window) and stays orthonormal; virtmix lands
    in the same basin as free CASSCF (same-basin assertion is the
    regression guard for the ticket-09 bug class)."""
    case, r = LADDER_CASES[name], ladders[name]
    b = case["bands"]
    lo, hi = b["corr_gain"]
    for tag in ("casscf_shell", "casscf_virtmix", "casscf_pyscf"):
        rung = r[tag]
        assert rung["ground_term"] == "3P", (name, tag)
        assert rung["ortho_err"] < 1e-10, (name, tag)
        assert 1e-4 < rung["comm_norm_h_lz"] < 1.0, (name, tag)
        assert 0.9 < rung["ground_l2"] < 2.0, (name, tag)
        assert rung["window"]["p"] > 0.9, (name, tag)
        d_e = rung["e_ground"] - r["casci_avg"]["e_ground"]
        assert lo < d_e < hi, (name, tag)
    _assert_same_basin(r, ("casscf_shell", "casscf_virtmix", "casscf_pyscf"),
                       r, "casscf_shell", b["basin_tol"])


@pytest.mark.parametrize("name", list(LADDER_CASES))
def test_sa_restoration(name, ladders):
    """SA over the term's M_L components restores L to machine precision."""
    case, r = LADDER_CASES[name], ladders[name]
    b = case["bands"]
    for tag in ("sa", "sa_pyscf"):
        if tag not in r:
            continue
        sa = r[tag]
        assert sa["ortho_err"] < 1e-10, (name, tag)
        assert sa["comm_norm_h_lz"] < b["sa_comm"], (name, tag)
        assert sa["ground_l2"] > b["sa_l2"], (name, tag)
        assert sa["terms"] == {"3P": 3}, (name, tag)
        if b["sa_nroots"]:
            assert len(sa["energies_pyscf"]) == 3, (name, tag)


# --------------------------------------------------------------------------
# Sc (ticket 11): the transition-metal ladder. Same machinery, different
# physics -- every assertion below is the measured Sc pattern across three
# independent cold runs (see AGENTS.md fact 15).
# --------------------------------------------------------------------------

@pytest.fixture(scope="session")
def sc_ladder(tmp_path_factory):
    return run_ladder(
        AtomConfig("Sc", "ccpvdz", spin=1), ["3d", "4s"], 3,
        workdir=str(tmp_path_factory.mktemp("ladder_Sc")),
        sa_weights=[1 / 5] * 5, frac_nelec_open=3, frac_n_open=6)


def test_sc_casci_exact_closure_and_term_table(sc_ladder):
    """Raw-shell CAS(3,6) over {3d, 4s} is exactly rotation-closed (a
    complete d m-set + one s): averaged/frac cores commute at 1e-13 and
    the full 90-state sector reproduces the free-ion term tables of
    3d1 4s2, 3d2 4s1 and 3d3. The RAW core field is percent-level on a
    TM (vs 1e-5 on C/Si) but still a bounded perturbation."""
    r = sc_ladder
    for tag in ("casci_avg", "casci_frac"):
        assert r[tag]["comm_norm_h_lz"] < 1e-10, tag
        assert r[tag]["comm_norm_h_l2"] < 1e-10, tag
        assert r[tag]["terms"] == SC_EXACT_TERMS, tag
        assert r[tag]["ground_term"] == "2D", tag
        assert r[tag]["ground_l2"] == pytest.approx(6.0, abs=1e-8), tag
        w = r[tag]["window"]
        assert abs(w.get("s", 0.0) - 1 / 6) < 1e-4, tag
        assert abs(w.get("d", 0.0) - 5 / 6) < 1e-4, tag

    raw = r["casci_raw"]
    assert 1e-3 < raw["comm_norm_h_lz"] < 1e-1
    assert 1e-2 < raw["comm_norm_h_l2"] < 0.5
    assert 6.0001 < raw["ground_l2"] < 6.001
    assert raw["ground_term"] == "2D"
    assert 1e-3 < raw["core_anisotropy"] < 1e-1
    # PySCF's own (unaveraged) fold splits the 2D quintet at the mEh scale
    e = np.sort(np.asarray(raw["energies_pyscf"]))
    assert 1e-4 < e[4] - e[0] < 1e-2

    # frac core: exactly invariant (PySCF's quintet degenerate to 1e-10),
    # a DIFFERENT reference (+30 mEh) whose softer field reshuffles the low
    # spectrum: the 2D terms of 3d3 / 3d1 4s2 / 3d2 4s1 become
    # quasi-degenerate at the bottom and their ORDERING flips across SCF
    # chaos (measured: d^3 vs s d^2 ground between identical reruns)
    frac = r["casci_frac"]
    ef = np.sort(np.asarray(frac["energies_pyscf"]))
    assert ef[4] - ef[0] < 1e-10
    assert 1e-3 < abs(frac["e_ground"] - r["casci_avg"]["e_ground"]) < 1e-1


def test_sc_relax_keeps_window_exact(sc_ladder):
    """Frozen window + relaxed core: the rung that stays honest on a TM --
    s+d window purity preserved exactly, sharp term table, only the core's
    own polarization remains (bounded commutator)."""
    relax = sc_ladder["casscf_relax"]
    assert relax["ortho_err"] < 1e-10
    assert relax["ground_term"] == "2D"
    assert 5.9990 < relax["ground_l2"] < 6.0002
    assert relax["terms"] == SC_EXACT_TERMS
    w = relax["window"]
    assert w.get("s", 0.0) + w.get("d", 0.0) > 0.999
    assert 1e-4 < relax["comm_norm_h_lz"] < 1e-2
    assert relax["max_de_vs_pyscf"] < 5e-3


def test_sc_free_casscf_escapes_the_shell(sc_ladder):
    """All three free-CASSCF starts land in one dominantly-L=0 basin
    (~52 mEh above casci_avg, 2S-CHARACTER state of the escaped window --
    an optimizer artifact, not a physical Sc term): the window borrows p
    character massively (s/p/d mix), L breaks at the 0.3 level -- but the
    optimizer stays orthonormal and variational (the ticket-09 guard
    holds; this is a real multi-solution, not a construction bug)."""
    r = sc_ladder
    for tag in ("casscf_virtmix", "casscf_shell", "casscf_pyscf"):
        rung = r[tag]
        assert rung["ortho_err"] < 1e-10, tag
        assert rung["ground_l2"] < 0.5, tag
        assert 0.1 < rung["comm_norm_h_lz"] < 1.0, tag
        assert rung["window"].get("p", 0.0) > 0.3, tag
        assert rung["max_de_vs_pyscf"] < 5e-3, tag
        assert not rung["diverged"], tag
        d_e = rung["e_ground"] - r["casci_avg"]["e_ground"]
        assert 0.03 < d_e < 0.08, tag
    _assert_same_basin(r, ("casscf_virtmix", "casscf_shell", "casscf_pyscf"),
                       r, "casscf_virtmix", 5e-3)


def test_sc_sa_restores_l2_not_lz(sc_ladder):
    """Shell-start SA over the five 2D M_L components restores <L^2> = 6 --
    but NOT [H, Lz]: the optimizer swaps the window's 4s member for a core
    3p (window p-weight exactly 1/6, core reorganized 4s+5p, anisotropy
    0.30): an m-INCOMPLETE window. Our averaged fold of that strongly
    polarized core vs PySCF's raw one shows as max|dE| ~ 0.8 (the
    polarized-core detector at Eh scale). The canonical-start SA collapses
    outright -- the start dependence is real on a TM."""
    r = sc_ladder
    sa = r["sa"]
    assert sa["ortho_err"] < 1e-10
    assert 5.99 < sa["ground_l2"] < 6.02
    assert sa["comm_norm_h_lz"] > 0.1
    assert 0.1 < sa["window"].get("p", 0.0) < 0.25
    assert 0.7 < sa["window"].get("d", 0.0) < 0.9
    assert 0.1 < sa["max_de_vs_pyscf"] < 5.0

    sa_p = r["sa_pyscf"]
    assert sa_p["ground_l2"] < 0.5
    assert sa_p["ortho_err"] < 1e-10
