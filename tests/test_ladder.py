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

from conftest import SC_EXACT_TERMS, TI_EXACT_TERMS, TI_OCC

from atomterm.drivers import build_full_mo, run_casci
from atomterm.ladder import run_ladder
from atomterm.stages import stage_facts
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
            radial_purity=0.999,        # the mask's structural guarantee
            sa_radial_comm=1e-8,        # measured 8e-11
            sa_radial_l2=1.99999,
            radial_free_tol=1e-5,       # radial == free minimum (measured 2e-8)
            sa_radial_dE=1e-4,          # vs the free SA (measured 2.7e-5)
            # purified-substrate rungs (fact 20): everything closes at
            # machine precision -- free SA included (measured 0.0)
            pur_comm=1e-10,
            pur_sa_comm=1e-8,
            pur_core_dE=2e-5,           # gate (d): measured 3.0e-6
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
            radial_purity=0.999,
            sa_radial_comm=1e-2,        # measured 1.6e-5 (gauge chaos headroom)
            sa_radial_l2=1.99,
            radial_free_tol=1e-3,       # radial == free (measured 3e-7; the
                                        # free rung's own basin wobbles)
            sa_radial_dE=5e-4,          # vs the free SA (measured 3.4e-5)
            # purified-substrate rungs (fact 20): machine precision here
            # too -- the free SA heals BELOW its canonical 4e-6 (fact 19b)
            # all the way to the certificate (measured 3.3e-16)
            pur_comm=1e-10,
            pur_sa_comm=1e-8,
            pur_core_dE=5e-4,           # gate (d): measured 1.28e-4
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
            frac_nelec_open=2, frac_n_open=3, purified=True)
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


def _assert_window(rung: dict, s: float, d: float, tol: float, tag: str = ""):
    """The rung's window carries exactly the given s/d l-weights."""
    w = rung["window"]
    assert abs(w.get("s", 0.0) - s) < tol, tag
    assert abs(w.get("d", 0.0) - d) < tol, tag


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


@pytest.mark.parametrize("name", list(LADDER_CASES))
def test_radial_rung_full_gain_without_l_leakage(name, ladders):
    """The same-(l, member) radial constraint (fact 16): the ENTIRE
    free-CASSCF correlation gain is radially reachable (radial rung ==
    free rung to µEh) with ZERO l-leakage (exact shell purity, the mask's
    structural guarantee). Which MEMBER gauge the single-root optimizer
    lands in is chaotic (C measured 5e-11, Si 1.2e-1 -- a flat direction),
    so the commutators are documented facts, not pinned here; sa_radial
    anchors the gauge (L exactly restored, commutators per-case bounded).
    """
    case, r = LADDER_CASES[name], ladders[name]
    b = case["bands"]
    rad = r["casscf_radial"]
    assert rad["ground_term"] == "3P", name
    assert rad["ortho_err"] < 1e-10 and not rad["diverged"], name
    assert rad["window"]["p"] > b["radial_purity"], name  # no l-leakage at all
    assert abs(rad["e_ground"] - r["casscf_shell"]["e_ground"]) < b["radial_free_tol"], name
    lo, hi = b["corr_gain"]
    assert lo < rad["e_ground"] - r["casci_avg"]["e_ground"] < hi, name
    # imposed-symmetry SA: same window family, L restored
    sar = r["sa_radial"]
    assert sar["ground_term"] == "3P", name
    assert sar["ortho_err"] < 1e-10, name
    assert sar["window"]["p"] > b["radial_purity"], name
    assert sar["ground_l2"] > b["sa_radial_l2"], name
    assert sar["comm_norm_h_lz"] < b["sa_radial_comm"], name
    assert abs(sar["e_ground"] - r["sa"]["e_ground"]) < b["sa_radial_dE"], name


@pytest.mark.parametrize("name", list(LADDER_CASES))
def test_purified_rungs_exact_closure(name, ladders):
    """Fact 20 (C/Si): on the purified core substrate EVERY optimized rung
    that keeps the window closed does so at machine precision -- the
    fact-17 input baseline (Si 3.4e-3) is gone. The FREE SA closes too:
    on p2 atoms the radial direction has a restoring force (fact 18a), so
    SA's stationary point is the symmetric window. casci_avg_pur
    additionally shows max|dE| = 0 (the purified core needs no averaging:
    PySCF's own fold and ours coincide) and gate (d): the purification
    shifts the reference by the discarded core polarization only."""
    case, r = LADDER_CASES[name], ladders[name]
    b = case["bands"]
    ci = r["casci_avg_pur"]
    assert ci["comm_norm_h_lz"] < b["pur_comm"], name
    assert ci["max_de_vs_pyscf"] < 1e-8, name          # measured 0.0
    assert ci["core_anisotropy"] < 1e-12, name          # m-scalar by build
    assert ci["terms"] == {"3P": 3}, name
    assert ci["ground_l2"] == pytest.approx(2.0, abs=1e-9), name
    for tag, comm_b in [("casscf_relax_pur", b["pur_comm"]),
                        ("sa_relax_pur", b["pur_comm"]),
                        ("sa_pur", b["pur_sa_comm"]),
                        ("sa_radial_pur", b["pur_sa_comm"])]:
        rung = r[tag]
        assert rung["comm_norm_h_lz"] < comm_b, (name, tag)
        assert rung["ground_term"] == "3P", (name, tag)
        assert rung["ground_l2"] == pytest.approx(2.0, abs=1e-6), (name, tag)
        assert rung["ortho_err"] < 1e-10 and not rung["diverged"], (name, tag)
        assert rung["window"]["p"] > 0.999, (name, tag)
    # the SA family reaches the full correlation gain (fact 18a optimum)
    lo, hi = b["corr_gain"]
    for tag in ("sa_pur", "sa_radial_pur"):
        assert lo < r[tag]["e_ground"] - r["casci_avg"]["e_ground"] < hi, \
            (name, tag)
    # single-root radial pur twin: same flat member gauge as the canonical
    # rung (fact 16a) -- purity and gain are structural, commutators not
    rad = r["casscf_radial_pur"]
    assert rad["ground_term"] == "3P", (name, rad["ground_term"])
    assert rad["window"]["p"] > 0.999, name
    assert rad["ortho_err"] < 1e-10 and not rad["diverged"], name
    assert lo < rad["e_ground"] - r["casci_avg"]["e_ground"] < hi, name
    # gate (d): same window, canonical vs purified core -- the reference
    # moved by the core polarization energy only (not a model change)
    d_e = r["casci_avg_pur"]["energies_pyscf"][0] \
        - r["casci_avg"]["energies_pyscf"][0]
    assert 0.0 < d_e < b["pur_core_dE"], (name, d_e)


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
        sa_weights=[1 / 5] * 5, frac_nelec_open=3, frac_n_open=6,
        purified=True)


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


def test_sc_purified_rungs_unbroken_lz_and_l2(sc_ladder):
    """Fact 20, the ticket-13 deliverable: Sc CASSCF and SA-CASSCF with
    UNBROKEN Lz and L^2. On the purified substrate the input-window
    baseline (5.7e-3 / 2.1e-2, fact 17) vanishes, so:

    - relax_pur / sa_relax_pur (window frozen): commutators at MACHINE
      precision (measured 1.3e-13 / 6.7e-14, 3.6e-13) with the full
      90-state free-ion term table and max|dE| ~ 0 -- the rung pair that
      answers "SA-CASSCF with exact rotational symmetry on a TM";
    - sa_radial_pur: the flat member gauge of fact 18d has NO restoring
      force (m-symmetric radial rotation gains nothing), so the optimizer
      stops at a ~1e-6 window tilt (measured comm 5.7e-7 / 1.8e-5,
      convergence-tolerance-independent -- an optimizer floor, not a
      baseline);
    - sa_pur (free): the 4s<->3p swap ATTRACTOR persists on the purified
      substrate -- it is an optimizer basin, not a core-span artifact
      (bimodal across SCF chaos: branch on the window p-weight, as
      always)."""
    r = sc_ladder
    ci = r["casci_avg_pur"]
    assert ci["comm_norm_h_lz"] < 1e-10 and ci["comm_norm_h_l2"] < 1e-10
    assert ci["terms"] == SC_EXACT_TERMS
    assert ci["max_de_vs_pyscf"] < 1e-8          # measured 0.0
    assert ci["core_anisotropy"] < 1e-12         # m-scalar by construction
    # gate (d) on Sc: +0.67 mEh = the [Ar] percent-polarization energy
    d_e = ci["energies_pyscf"][0] - r["casci_avg"]["energies_pyscf"][0]
    assert 0.0 < d_e < 2e-3, d_e

    for tag, comm_b in [("casscf_relax_pur", 1e-10),
                        ("sa_relax_pur", 1e-10)]:
        rung = r[tag]
        assert rung["comm_norm_h_lz"] < comm_b, tag
        assert rung["comm_norm_h_l2"] < comm_b, tag
        assert rung["ground_l2"] == pytest.approx(6.0, abs=1e-6), tag
        assert rung["terms"] == SC_EXACT_TERMS, tag
        assert rung["ground_term"] == "2D", tag
        assert rung["ortho_err"] < 1e-10 and not rung["diverged"], tag
        w = rung["window"]
        assert abs(w.get("s", 0.0) - 1 / 6) < 1e-6, tag
        assert abs(w.get("d", 0.0) - 5 / 6) < 1e-6, tag
    # relaxation moves the energy by the core-polarization scale, but its
    # SIGN in our bookkeeping is not physics: the optimizer buys
    # polarization energy in PySCF's RAW fold (measured relax +0.47 mEh,
    # sa_relax -0.09 mEh vs the CASCI in our averaged fold) -- the
    # max|dE| ~ 1e-3 column is precisely that detector. Bound, don't pin.
    for tag in ("casscf_relax_pur", "sa_relax_pur"):
        assert abs(r[tag]["e_ground"] - ci["e_ground"]) < 5e-3, tag
    assert r["sa_relax_pur"]["max_de_vs_pyscf"] < 1e-5  # measured 4.8e-8

    sar = r["sa_radial_pur"]
    assert sar["comm_norm_h_lz"] < 1e-4, sar["comm_norm_h_lz"]  # meas 5.7e-7
    assert sar["comm_norm_h_l2"] < 1e-3, sar["comm_norm_h_l2"]  # meas 1.8e-5
    assert sar["ground_l2"] == pytest.approx(6.0, abs=1e-4), sar
    assert sar["terms"] == SC_EXACT_TERMS, sar
    assert sar["ortho_err"] < 1e-10 and not sar["diverged"], sar

    sa = r["sa_pur"]
    assert sa["ortho_err"] < 1e-10
    assert 5.99 < sa["ground_l2"] < 6.02
    if sa["window"].get("p", 0.0) > 0.05:  # the swap basin (measured p=1/6)
        assert sa["comm_norm_h_lz"] > 0.1, sa
        assert abs(sa["window"].get("p", 0.0) - 1 / 6) < 1e-2, sa["window"]
    else:  # the healthy basin
        assert sa["comm_norm_h_lz"] < 1e-3, sa


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
    """Shell-start SA over the five 2D M_L components restores <L^2> = 6
    always, and is BIMODAL across SCF chaos (fact 15 amended): either it
    swaps the window's 4s member for a core 3p (window p-weight exactly
    1/6, an m-INCOMPLETE window, [H, Lz] ~ 0.3; our averaged fold of that
    strongly polarized core vs PySCF's raw one shows as max|dE| ~ 0.8 --
    the polarized-core detector at Eh scale), or it converges healthy
    (window {3d, 4s}, commutators at the ~1e-3 input-window baseline of
    fact 17, like sa_radial). The canonical-start SA collapses to the
    dominantly-L=0 basin -- the start dependence is real on a TM."""
    r = sc_ladder
    sa = r["sa"]
    assert sa["ortho_err"] < 1e-10
    assert 5.99 < sa["ground_l2"] < 6.02
    if sa["window"].get("p", 0.0) > 0.05:  # the swapped-window basin
        assert sa["comm_norm_h_lz"] > 0.1
        assert 0.1 < sa["window"].get("p", 0.0) < 0.25
        assert 0.7 < sa["window"].get("d", 0.0) < 0.9
        assert 0.1 < sa["max_de_vs_pyscf"] < 5.0
    else:  # the healthy basin (observed in cold runs after ticket 14)
        assert sa["comm_norm_h_lz"] < 5e-3
        assert abs(sa["window"].get("s", 0.0) - 1 / 6) < 1e-2
        assert abs(sa["window"].get("d", 0.0) - 5 / 6) < 1e-2

    sa_p = r["sa_pyscf"]
    assert sa_p["ground_l2"] < 0.5
    assert sa_p["ortho_err"] < 1e-10


def test_sc_radial_rung_blocks_the_escape(sc_ladder):
    """Fact 16 on the TM: the same-(l, member) mask keeps the window
    exactly {3d 5/6, 4s 1/6} (no p-route to the escape), so the radial
    ceiling sits ~25 mEh ABOVE the escaped free minimum (stable across
    cold runs: +24.4, +25.1) -- the escape is an l/window-identity
    phenomenon, not radially reachable correlation. The single-root radial
    GAUGE is chaotic (cold runs: comm 2e-3..1.6, <L^2> 3.5..5.4, the
    ground term even relabels 2D->2P -- a flat member direction, exactly
    like the free rung's m-breaking), so only its energy is pinned.
    sa_radial is the anchored counterpart: reproducible to µEh, <L^2> =
    5.9997, commutators at the INPUT window's core-orthogonalization
    baseline (~1e-3, fact 17), NOT machine precision. All Sc CASSCF-family
    energies carry the window-definition offset vs casci_avg (fact 14)."""
    r = sc_ladder
    rad, sar, shell = r["casscf_radial"], r["sa_radial"], r["casscf_shell"]
    for rung in (rad, sar):
        assert abs(rung["window"].get("s", 0.0) - 1 / 6) < 1e-3
        assert abs(rung["window"].get("d", 0.0) - 5 / 6) < 1e-3
        assert rung["window"].get("p", 0.0) < 1e-3
        assert rung["ortho_err"] < 1e-10 and not rung["diverged"]
        assert 0.03 < rung["e_ground"] - r["casci_avg"]["e_ground"] < 0.09
    # the escape is BELOW the radial ceiling (both minima are well-defined
    # even though the single-root radial gauge is not; observed gap ~25 mEh,
    # band widened for the escape energy's own wobble, fact 10)
    assert 0.0 < rad["e_ground"] - shell["e_ground"] < 0.05
    # imposed-symmetry SA: L restored, commutators at the 1e-3 baseline
    assert sar["ground_term"] == "2D"
    assert sar["ground_l2"] == pytest.approx(6.0, abs=1e-3)
    assert sar["comm_norm_h_lz"] < 1e-2  # meas 3.3e-3; fresh draws to 5.3e-3


@pytest.fixture(scope="session")
def sc_big_window_ladder(tmp_path_factory):
    """CAS(3,9) over the complete {3d, 4s, 4p} window (324 dets): is the
    escape just 4p correlation hunger?"""
    return run_ladder(
        AtomConfig("Sc", "ccpvdz", spin=1), ["3d", "4s", "4p"], 3,
        workdir=str(tmp_path_factory.mktemp("ladder_Sc_cas311")),
        sa_weights=[1 / 5] * 5, rung_tags=["casci_avg", "casscf_shell"])


def test_sc_complete_window_still_escapes(sc_ladder, sc_big_window_ladder):
    """With 4p IN the window the CASCI gains ~25 mEh of static correlation
    at EXACT closure (complete shells) -- yet the free CASSCF still leaves
    for the dominantly-L=0 basin (window grabs even f character): the
    escape is an optimizer attractor, not a physical correlation target,
    and a complete START window is not a complete CONVERGED window."""
    small, big = sc_ladder, sc_big_window_ladder
    ci = big["casci_avg"]
    assert ci["comm_norm_h_lz"] < 1e-10 and ci["comm_norm_h_l2"] < 1e-10
    assert ci["ground_term"] == "2D"
    assert ci["ground_l2"] == pytest.approx(6.0, abs=1e-8)
    assert abs(ci["window"].get("p", 0.0) - 3 / 9) < 1e-3
    d_static = ci["e_ground"] - small["casci_avg"]["e_ground"]
    assert -0.035 < d_static < -0.015  # ~-25 mEh of 4p static correlation
    shell = big["casscf_shell"]
    assert shell["ground_l2"] < 1.5  # escaped anyway
    assert not shell["diverged"] and shell["ortho_err"] < 1e-10


# ==========================================================================
# Ti (ticket 12, option B): CAS(4,6)/{3d, 4s}, Ms=1 sector (120 dets) --
# does the whole settled story (fact 15 TM phenomenology, fact 18 scans,
# fact 20 purification) transfer to a 4-valence-electron d2 TM?
# ==========================================================================


@pytest.fixture(scope="session")
def ti_ladder(tmp_path_factory):
    return run_ladder(
        AtomConfig("Ti", "ccpvdz", spin=2), ["3d", "4s"], 4,
        workdir=str(tmp_path_factory.mktemp("ladder_Ti")),
        sa_weights=[1 / 7] * 7, frac_nelec_open=4, frac_n_open=6,
        purified=True, **TI_OCC)


def test_ti_casci_exact_closure_and_term_table(ti_ladder):
    """Raw-shell CAS(4,6) over {3d, 4s} is exactly rotation-closed: the
    averaged/frac cores commute at 1e-13 and the 120-state sector
    reproduces the d2 + d3 s + d4 term tables (ground 3F, the first F
    ground term). The RAW core field is percent-level (like Sc) and the
    frac reference is a different (+38 mEh) but exactly invariant one."""
    r = ti_ladder
    for tag in ("casci_avg", "casci_frac"):
        assert r[tag]["comm_norm_h_lz"] < 1e-10, tag
        assert r[tag]["comm_norm_h_l2"] < 1e-10, tag
        assert r[tag]["terms"] == TI_EXACT_TERMS, tag
        assert r[tag]["ground_term"] == "3F", tag
        assert r[tag]["ground_l2"] == pytest.approx(12.0, abs=1e-8), tag
        _assert_window(r[tag], 1 / 6, 5 / 6, 1e-4, tag)

    raw = r["casci_raw"]
    assert 1e-3 < raw["comm_norm_h_lz"] < 1e-1
    assert 1e-2 < raw["comm_norm_h_l2"] < 0.5
    # unlike Sc (whose raw core shifts <L^2> to 6.0001), Ti's weaker raw
    # polarization (anisotropy 7.8e-3 vs Sc's 1.3e-2) leaves the ground
    # <L^2> essentially unshifted (measured 12.0 - 1e-13)
    assert raw["ground_l2"] == pytest.approx(12.0, abs=1e-3)
    assert raw["ground_term"] == "3F"
    assert 1e-3 < raw["core_anisotropy"] < 1e-1
    # PySCF's own (unaveraged) fold splits the 3F septet at the mEh scale
    # (measured 4.1e-3)
    e = np.sort(np.asarray(raw["energies_pyscf"]))
    assert 1e-4 < e[6] - e[0] < 1e-2

    # frac core: exactly invariant (septet degenerate to its numerical
    # floor, measured 6.5e-11), a DIFFERENT reference (+38 mEh measured;
    # band widened for the SCF-chaos wobble)
    frac = r["casci_frac"]
    ef = np.sort(np.asarray(frac["energies_pyscf"]))
    assert ef[6] - ef[0] < 1e-8
    assert 1e-3 < abs(frac["e_ground"] - r["casci_avg"]["e_ground"]) < 1e-1


def test_ti_relax_keeps_window_exact(ti_ladder):
    """Frozen window + relaxed core: on Ti (reference gated to the clean
    ROHF basin by occ_e_ref) the relax rung is nearly closed -- measured
    2.3e-8, versus Sc's fact-17 ~1e-3 baseline: Ti's good-basin input
    window carries only a ~1e-8 core-orthogonalization defect (the dirty
    +11.2 mEh SCF basin would show 1e-2-class defects; the basin gate
    keeps tests deterministic)."""
    relax = ti_ladder["casscf_relax"]
    assert relax["ortho_err"] < 1e-10
    assert relax["ground_term"] == "3F"
    assert relax["ground_l2"] == pytest.approx(12.0, abs=1e-3)  # meas 12.000000
    assert relax["terms"] == TI_EXACT_TERMS
    w = relax["window"]
    assert w.get("s", 0.0) + w.get("d", 0.0) > 0.999
    assert relax["comm_norm_h_lz"] < 1e-6  # measured 2.3e-8
    assert relax["max_de_vs_pyscf"] < 5e-3


def test_ti_purified_rungs_unbroken_lz_and_l2(ti_ladder):
    """Fact 20's deliverable transfers to the d2 atom: on the purified
    [Ar] substrate the frozen-window rungs close at machine precision
    (relax_pur 6.2e-14, sa_relax_pur 1.4e-13/8.1e-13) with the full
    120-state term table, <L^2> = 12 and unbroken Lz -- Ti CASSCF and
    SA-CASSCF with exact rotational symmetry. sa_radial_pur sits at the
    flat-gauge optimizer floor (fact 21: d2 tilts DESCEND, so no
    restoring force tightens the member gauge; measured 1.4e-4, slightly
    above Sc's 5.7e-7). The free-sa swap attractor survives the purified
    substrate (bimodal across SCF chaos: branch on window p-weight)."""
    r = ti_ladder
    ci = r["casci_avg_pur"]
    assert ci["comm_norm_h_lz"] < 1e-10 and ci["comm_norm_h_l2"] < 1e-10
    assert ci["terms"] == TI_EXACT_TERMS
    assert ci["max_de_vs_pyscf"] < 1e-8          # measured 1.1e-13
    assert ci["core_anisotropy"] < 1e-12         # m-scalar by construction
    # the pur reference shift: the discarded core polarization (+1.36 mEh
    # measured in PySCF's own folds -- same scale as Sc's +0.67 mEh)
    d_e = ci["energies_pyscf"][0] - r["casci_avg"]["energies_pyscf"][0]
    assert 0.0 < d_e < 5e-3, d_e

    for tag in ("casscf_relax_pur", "sa_relax_pur"):
        rung = r[tag]
        assert rung["comm_norm_h_lz"] < 1e-10, tag
        assert rung["comm_norm_h_l2"] < 1e-10, tag
        assert rung["ground_l2"] == pytest.approx(12.0, abs=1e-6), tag
        assert rung["terms"] == TI_EXACT_TERMS, tag
        assert rung["ground_term"] == "3F", tag
        assert rung["ortho_err"] < 1e-10 and not rung["diverged"], tag
        _assert_window(rung, 1 / 6, 5 / 6, 1e-6, tag)
        assert abs(rung["e_ground"] - ci["e_ground"]) < 5e-3, tag
    assert r["sa_relax_pur"]["max_de_vs_pyscf"] < 1e-5  # measured 8e-8

    sar = r["sa_radial_pur"]
    assert sar["comm_norm_h_lz"] < 1e-3, sar["comm_norm_h_lz"]  # meas 1.4e-4
    assert sar["comm_norm_h_l2"] < 1e-2, sar["comm_norm_h_l2"]  # meas 5.3e-4
    assert sar["ground_l2"] == pytest.approx(12.0, abs=1e-4), sar
    assert sar["terms"] == TI_EXACT_TERMS, sar
    assert sar["ortho_err"] < 1e-10 and not sar["diverged"], sar

    sa = r["sa_pur"]
    assert sa["ortho_err"] < 1e-10
    # SA restores L in both basins; the swap mixes slightly (observed
    # 12.011..12.025 across fresh ROHF draws)
    assert sa["ground_l2"] == pytest.approx(12.0, abs=0.05)
    if sa["window"].get("p", 0.0) > 0.03:  # the swap basin (measured p=1/18)
        assert sa["comm_norm_h_lz"] > 0.1, sa
    else:  # the healthy basin
        assert sa["comm_norm_h_lz"] < 1e-3, sa


def test_ti_free_casscf_escapes_the_shell(ti_ladder):
    """All three free-CASSCF starts land in one dominantly-L=0 basin
    (+46 mEh above casci_avg, the window borrowing p massively --
    s/p/d-mixed): the Sc escape phenomenology transfers to d2 unchanged."""
    r = ti_ladder
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


def test_ti_sa_basins_and_radial_rung(ti_ladder):
    """Shell-start SA over the seven 3F M_L components: <L^2> = 12 always;
    bimodal across SCF chaos -- either the exact 4s<->3p swap basin of Sc
    (window p-weight EXACTLY 1/6, s ~ 0, [H, Lz] ~ 1) or healthy. The
    radial-masked rungs keep the window exactly {3d, 4s}: sa_radial is
    reproducible with commutators at the canonical-substrate baseline
    (~1e-3 class) and gains ~11 mEh over relax -- on Ti that gain lives on
    the ORTHOGONALIZED substrate's s slot (a diffuse 6s admixture, fact
    21), NOT in raw radial rotations (raw 3d->4d and 4s->5s/6s all RISE).
    The single-root radial gauge is chaotic (like fact 16c on Sc)."""
    r = ti_ladder
    sa = r["sa"]
    assert sa["ortho_err"] < 1e-10
    # basin-mixed <L^2>: same headroom as sa_pur (the ticket's own lesson)
    assert sa["ground_l2"] == pytest.approx(12.0, abs=0.05)
    if sa["window"].get("p", 0.0) > 0.05:  # the swapped-window basin
        assert sa["comm_norm_h_lz"] > 0.1
        assert abs(sa["window"].get("p", 0.0) - 1 / 6) < 1e-2
        assert 0.7 < sa["window"].get("d", 0.0) < 0.9
        assert 0.1 < sa["max_de_vs_pyscf"] < 5.0
    else:  # the healthy basin
        assert sa["comm_norm_h_lz"] < 5e-3
        assert abs(sa["window"].get("s", 0.0) - 1 / 6) < 1e-2
        assert abs(sa["window"].get("d", 0.0) - 5 / 6) < 1e-2

    sa_p = r["sa_pyscf"]
    assert sa_p["ground_l2"] < 0.5
    assert sa_p["ortho_err"] < 1e-10

    for tag in ("casscf_radial", "sa_radial"):
        rung = r[tag]
        assert abs(rung["window"].get("s", 0.0) - 1 / 6) < 1e-3, tag
        assert abs(rung["window"].get("d", 0.0) - 5 / 6) < 1e-3, tag
        assert rung["window"].get("p", 0.0) < 1e-3, tag
        assert rung["ortho_err"] < 1e-10 and not rung["diverged"], tag
    # the radial gain over the frozen relax window (~11 mEh, orthogonalized
    # substrate; the raw-window fold would show every one-parameter radial
    # rotation RISING -- fact 21)
    assert -0.02 < r["sa_radial"]["e_ground"] - r["casscf_relax"]["e_ground"] < 0.0
    sar = r["sa_radial"]
    assert sar["ground_term"] == "3F"
    assert sar["ground_l2"] == pytest.approx(12.0, abs=1e-3)
    assert sar["comm_norm_h_lz"] < 5e-3
