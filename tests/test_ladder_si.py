"""Si ladder tests (ticket 09): the "escaped 3S basin" resolved.

Ticket 07 read Si's CASSCF/SA collapse as a competing 3S well (CASSCF
multi-solution, fact 13 then). Ticket 09 traced it to a CONSTRUCTION BUG:
build_full_mo silently returned a non-orthonormal MO matrix on Si (raw 3p
AOs overlap the canonical 2p core MOs by ~0.3, so the complement projector
was inexact and extra junk columns appeared). Fed that matrix, the orbital
optimizer diverged -- unphysical energies below the variational bound.
With the block-safe construction every rung converges and Si behaves
exactly like carbon: CASCI textbook-exact, CASSCF in the 3P basin with a
real correlation gain, SA restoring <L^2> = 2.000000.
"""

import numpy as np
import pytest

from atomterm.ladder import run_ladder
from atomterm.system import AtomConfig


@pytest.fixture(scope="session")
def si_ladder(tmp_path_factory):
    cfg = AtomConfig("Si", "ccpvdz", spin=2)
    return cfg, run_ladder(
        cfg, ["3p"], 2,
        workdir=str(tmp_path_factory.mktemp("ladder_si")),
        sa_weights=[1 / 3] * 3, frac_nelec_open=2, frac_n_open=3)


def test_si_casci_textbook(si_ladder):
    """All CASCI rungs: exact 3p^2 tables, invariance per core treatment."""
    _, r = si_ladder
    for tag in ("casci_raw", "casci_avg", "casci_frac"):
        rung = r[tag]
        assert rung["terms"] == {"3P": 3}
        assert rung["ortho_err"] < 1e-10
        # configuration labels read the p^2 character off the determinants
        assert rung["configurations"][0].startswith("p^2.00")
    assert r["casci_avg"]["comm_norm_h_lz"] < 1e-8
    assert r["casci_frac"]["comm_norm_h_lz"] < 1e-8
    # the raw core field is a larger perturbation on Si than on C
    # (bigger, more polarizable core): measured ~5e-4 vs carbon ~4e-6
    assert 1e-5 < r["casci_raw"]["comm_norm_h_lz"] < 1e-2
    # 3P split, {M_L=0} off the {M_L=+-1} pair
    roots = np.sort(r["casci_raw"]["energies_pyscf"])
    split = roots[2] - roots[0]
    assert 1e-5 < split < 1e-3
    assert abs(roots[2] - roots[1]) < split / 10


def test_si_casscf_rungs_converge_in_3p_basin(si_ladder):
    """Free CASSCF (both starts) and virtmix: same 3P basin, real
    correlation gain, no divergence -- the ticket-07 escape is gone."""
    _, r = si_ladder
    e_pyscf_start = r["casscf_pyscf"]["e_ground"]
    e_ref = r["casci_avg"]["e_ground"]
    for tag in ("casscf_shell", "casscf_virtmix", "casscf_pyscf"):
        rung = r[tag]
        assert rung["ortho_err"] < 1e-10, tag
        assert rung["ground_term"] == "3P", tag
        # correlation gain vs the frozen-orbital CASCI, same all starts
        d_e = rung["e_ground"] - e_ref
        assert -0.2 < d_e < -0.05, tag
        assert abs(rung["e_ground"] - e_pyscf_start) < 1e-3, tag
        assert 1.0 < rung["ground_l2"] < 2.0, tag
        assert rung["window"]["p"] > 0.9, tag


def test_si_sa_restoration(si_ladder):
    """SA over the three 3P M_L components restores L on Si too."""
    _, r = si_ladder
    sa = r["sa_pyscf"]
    assert sa["ortho_err"] < 1e-10
    assert sa["ground_l2"] > 1.9999
    assert sa["comm_norm_h_lz"] < 1e-7
    assert sa["terms"] == {"3P": 3}


def test_si_core_relaxation(si_ladder):
    """Core-relaxation rung: converges, symmetry exact, window l-pure.

    Note the cross-rung dE on Si carries a window-definition offset (~0.09
    Eh): Si's raw 3p AOs overlap the 2p core MOs by ~0.3, so PySCF's
    orthogonalized window and our raw-shell analysis window differ
    physically. Carbon has no such offset (overlaps ~1e-6).
    """
    _, r = si_ladder
    relax = r["casscf_relax"]
    assert relax["ortho_err"] < 1e-10
    assert relax["comm_norm_h_lz"] < 1e-3
    assert relax["ground_l2"] == pytest.approx(2.0, abs=1e-4)
    assert relax["window"]["p"] > 0.95
    assert relax["configurations"][0].startswith("p^2.00")
