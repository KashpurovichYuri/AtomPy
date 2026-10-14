"""Carbon ladder tests (ticket 08): pattern assertions per rung.

The ladder reruns every orbital treatment on carbon and asserts the
MEASURED symmetry pattern, not assumed exactness:

- raw core: multiplet split and [H, Lz] ~ eps are PRESENT (the internal
  field is real) but bounded;
- averaged / fractional-occupation cores: exact invariance;
- free CASSCF: symmetry breaking is real but window escape is NOT (the
  window stays l-pure -- fact 8's breaking is m-incompleteness, not s/d
  leakage);
- SA(3): restoration to machine precision.

Numbers are linted loosely against chaotic CASSCF/BLAS reproduction
(AGENTS.md fact 10); the ticket carries the precisely measured values.
"""

from collections import Counter

import numpy as np
import pytest

from atomterm.ladder import run_carbon_ladder
from atomterm.system import AtomConfig


@pytest.fixture(scope="session")
def ladder(tmp_path_factory):
    cfg = AtomConfig()
    return cfg, run_carbon_ladder(
        cfg, workdir=str(tmp_path_factory.mktemp("ladder")))


def test_casci_raw_core_fingerprint(ladder):
    """The un-averaged core field breaks symmetry measurably but mildly."""
    _, r = ladder
    raw = r["casci_raw"]
    # 3P term letters survive the ~1e-5 perturbation
    assert raw["terms"] == {"3P": 3}
    # [H, Lz] nonzero (the internal field) yet small
    assert 1e-7 < raw["comm_norm_h_lz"] < 1e-3
    # [H, L^2] is trivially zero here: a complete L=1 irrep, L^2 = 2*Id
    assert raw["comm_norm_h_l2"] < 1e-12
    # the 3P multiplet is split: {M_L=0} off the {M_L=+-1} pair
    roots = np.sort(raw["energies_pyscf"])
    split = roots[2] - roots[0]
    pair = abs(roots[2] - roots[1])
    assert 1e-6 < split < 1e-3
    assert pair < split / 10
    # and the discarded anisotropy is a measured small quantity: the
    # density norm carries the quadrupole at ~2e-3 while the energy feels
    # it at ~4e-5 (second-order attenuation)
    assert 1e-4 < raw["core_anisotropy"] < 1e-1


def test_casscf_virtual_mix(ladder):
    """Virtmix rung: window<->virtual correlation without core<->window
    exchange lands in the same basin as free CASSCF (the Si lesson: that
    exchange, not virtual borrowing, was the identity-changer -- once the
    window construction bug was fixed, even that vanished)."""
    _, r = ladder
    vm = r["casscf_virtmix"]
    assert vm["ortho_err"] < 1e-10
    assert vm["ground_term"] == "3P"
    assert abs(vm["e_ground"] - r["casscf_shell"]["e_ground"]) < 1e-2
    assert 1e-4 < vm["comm_norm_h_lz"] < 1.0


def test_casci_avg_core_exact(ladder):
    """Spherical averaging restores exact invariance."""
    _, r = ladder
    avg = r["casci_avg"]
    assert avg["comm_norm_h_lz"] < 1e-9
    assert avg["comm_norm_h_l2"] < 1e-9
    assert avg["terms"] == {"3P": 3}
    # the averaged fold shifts the spectrum by exactly the splitting scale
    assert 1e-6 < abs(r["casci_avg"]["e_ground"] - r["casci_raw"]["e_ground"]) < 1e-3


def test_casci_frac_occ_core(ladder):
    """Average-of-configuration core: invariant by construction.

    Its energy does NOT match the averaged-ROHF core at O(eps^2) (the
    ticket's original prediction): average-of-configuration is a different
    reference -- open-shell exchange with smeared vs localized valence --
    differing at the mEh scale. Both cores are rotationally invariant; the
    difference is physical, measured, and recorded (AGENTS.md).
    """
    _, r = ladder
    frac = r["casci_frac"]
    assert frac["comm_norm_h_lz"] < 1e-9
    assert frac["comm_norm_h_l2"] < 1e-9
    assert frac["terms"] == {"3P": 3}
    assert 1e-6 < abs(frac["e_ground"] - r["casci_avg"]["e_ground"]) < 1e-1


def test_casscf_core_relaxation(ladder):
    """Frozen window + relaxed core: real physics, zero escape, exact symmetry."""
    _, r = ladder
    relax = r["casscf_relax"]
    assert relax["comm_norm_h_lz"] < 1e-9
    assert relax["ortho_err"] < 1e-10
    assert relax["ground_l2"] == pytest.approx(2.0, abs=1e-6)
    assert relax["terms"] == {"3P": 3}
    # relaxation lowers the energy vs frozen-core CASCI ...
    assert relax["e_ground"] <= r["casci_avg"]["e_ground"] + 1e-9
    # ... while the window stays exactly the raw shell (no span change)
    assert r["casscf_relax"]["window"]["p"] > 0.999


def test_casscf_free_symmetry_breaking(ladder):
    """Free CASSCF breaks symmetry WITHOUT leaving the p shell.

    Refined fact 8: the breaking is m-incompleteness of the optimized
    window (radial p mixtures lacking the full m-set), NOT s/d leakage --
    the window stays l-pure, and the commutators/ <L^2> measure the
    remaining m-openness.
    """
    _, r = ladder
    for tag in ("casscf_shell", "casscf_pyscf"):
        rung = r[tag]
        assert rung["ground_term"] == "3P"
        # chaotic gauge (fact 10): comm 1e-3..1e-1, <L^2> observed
        # 0.99..1.92 across reruns -- bands assert "real but not
        # catastrophic" breaking and a still-mostly-p window
        assert 1e-4 < rung["comm_norm_h_lz"] < 1.0
        assert 0.9 < rung["ground_l2"] < 2.0
        assert rung["window"]["p"] > 0.95
        assert rung["ortho_err"] < 1e-10
    # both starts land in the same basin (same energy within chaos noise)
    assert abs(r["casscf_shell"]["e_ground"] - r["casscf_pyscf"]["e_ground"]) < 1e-4


def test_sa_restoration(ladder):
    """SA over the three 3P M_L components restores L to machine precision."""
    _, r = ladder
    sa = r["sa"]
    assert sa["comm_norm_h_lz"] < 1e-9
    assert sa["ground_l2"] > 1.99999
    assert len(sa["energies_pyscf"]) == 3
    assert sa["terms"] == {"3P": 3}
