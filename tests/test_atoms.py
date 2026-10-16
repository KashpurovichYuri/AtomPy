"""Table-driven end-to-end tests of the generic pipeline (ticket 06).

Each row of ``CASES`` runs one atom through :func:`atomterm.pipeline.run_atom`
and checks the textbook term structure, the exact commutator symmetries of the
raw-shell CASCI stages, and the particle number of the determinant space.
Adding an atom (Si, Sc, Ti, ...) must be a pure data change: a new row here,
no new code.
"""

from collections import Counter

import numpy as np
import pytest

from atomterm.pipeline import run_atom
from atomterm.system import AtomConfig

# ground: expected (term, L, S) of the lowest CASCI state in the primary
# sector; primary_tag / primary: sector chkfile tag and expected term -> count
# over the FULL sector spectrum (ms0_* likewise, None when the sector does not
# exist, e.g. odd active-electron doublets)
CASES = {
    "C": dict(
        cfg=AtomConfig("C", "ccpvdz", spin=2),
        patterns=["2p"],
        nelec_act=2,
        run_casscf=True,
        sa_weights=[1 / 3] * 3,
        ground=("3P", 1, 1),
        primary_tag="casci_ms1",
        primary={"3P": 3},
        ms0_tag="casci_ms0",
        ms0={"1S": 1, "1D": 5, "3P": 3},
    ),
    # ticket 09: with the build_full_mo construction fixed, Si's CASSCF/SA
    # converge in the 3P basin (the old "escape" was a non-orthonormal MO
    # matrix, see AGENTS.md fact 13) -- stages re-enabled and asserted.
    "Si": dict(
        cfg=AtomConfig("Si", "ccpvdz", spin=2),
        patterns=["3p"],
        nelec_act=2,
        run_casscf=True,
        sa_weights=[1 / 3] * 3,
        ground=("3P", 1, 1),
        primary_tag="casci_ms1",
        primary={"3P": 3},
        ms0_tag="casci_ms0",
        ms0={"1S": 1, "1D": 5, "3P": 3},
    ),
    # ticket 11: first doublet, first transition metal. The primary sector
    # (2,1) is the ONLY Ms sector (odd active count) and holds the complete
    # free-ion term table of 3d1 4s2 + 3d2 4s1 + 3d3. Single-root CASSCF is
    # off: on a TM it collapses to an 2S-like basin with the window escaped
    # into p -- the full rung-by-rung story lives in test_ladder.py.
    "Sc": dict(
        cfg=AtomConfig("Sc", "ccpvdz", spin=1),
        patterns=["3d", "4s"],
        nelec_act=3,
        run_casscf=False,
        sa_weights=None,
        nroots_primary=90,
        ground=("2D", 2, 0.5),
        primary_tag="casci_ms0p5",
        primary={"2D": 20, "4F": 14, "2F": 14, "4P": 6, "2G": 18,
                 "2P": 6, "2S": 1, "2H": 11},
        ms0_tag=None,
        ms0=None,
        # the (2,1) sector carries exactly S = 1/2 and S = 3/2
        s_values={0.75, 3.75},
    ),
    # ticket 12 option A: d2 -- the Ti+-ion-like CAS(2,5) over {3d} on the
    # [Ar]+4s2 closed core. First L=3 ground term (3F, Hund), first L=4
    # term (1G, in the Ms=0 sector), first basis with an f shell under a
    # d-analyzed load (4f AOs present; core f-tails ~1e-4). occ_* are the
    # MANDATORY character-forced occupations: default energy-ordered ROHF
    # on Ti is an aufbau runaway (+123.07 mEh, fact 21).
    "Ti": dict(
        cfg=AtomConfig("Ti", "ccpvdz", spin=2),
        patterns=["3d"],
        nelec_act=2,
        run_casscf=False,
        sa_weights=None,
        nroots_primary=10,
        occ_closed=["1s", "2s", "3s", "4s", "2p", "3p"],
        occ_open=["3d"],
        ground=("3F", 3, 1),
        primary_tag="casci_ms1",
        primary={"3F": 7, "3P": 3},
        ms0_tag="casci_ms0",
        ms0={"1S": 1, "1D": 5, "3P": 3, "1G": 9, "3F": 7},
        # the (2,0) sector carries exactly S = 1
        s_values={2.0},
    ),
}


@pytest.fixture(scope="session")
def atom_results(tmp_path_factory):
    """run_atom output per case, computed once per test session."""
    return {
        name: run_atom(
            case["cfg"],
            case["patterns"],
            case["nelec_act"],
            workdir=str(tmp_path_factory.mktemp(f"atom_{name}")),
            sa_weights=case["sa_weights"],
            run_casscf=case["run_casscf"],
            nroots_primary=case.get("nroots_primary", 9),
            occ_closed=case.get("occ_closed"),
            occ_open=case.get("occ_open"),
        )
        for name, case in CASES.items()
    }


@pytest.mark.parametrize("name", list(CASES))
def test_atom_term_tables(name, atom_results):
    """Full-sector term tables, commutator symmetries, N, and ground state."""
    case, results = CASES[name], atom_results[name]

    # every CASCI sector: textbook term counts, exact rotational symmetry of
    # the raw-shell space, correct particle number
    for tag, expected in [
        (case["primary_tag"], case["primary"]),
        (case["ms0_tag"], case.get("ms0")),
    ]:
        if expected is None:
            assert tag is None or tag not in results, (name, tag)
            continue
        res = results[tag]["analysis"]
        terms = Counter(s.term for s in res.states)
        assert terms == expected, (name, tag, terms)
        assert res.comm_norm_h_lz < 1e-8, (name, tag)
        assert res.comm_norm_h_l2 < 1e-8, (name, tag)
        assert results[tag]["n_dev"] < 1e-9, (name, tag)

    # ground state of the primary sector carries the Hund-rule term (L, S)
    ground = results[case["primary_tag"]]["analysis"].states[0]
    assert (ground.term, ground.l, ground.s) == case["ground"], (name, ground)

    # S^2 sharpness: the sector's spin content (ticket 11: "S^2 spectrum
    # contains 0.75" -- pinned as the full S content of the sector)
    if "s_values" in case:
        s2_set = {round(st.s2, 6) for st in results[case["primary_tag"]]["analysis"].states}
        assert s2_set == case["s_values"], (name, sorted(s2_set))

    # optimized-orbital stages keep the same ground term (their residual
    # symmetry breaking is physical output, asserted in test_pipeline)
    if case["run_casscf"]:
        assert results["casscf"]["analysis"].states[0].term == case["ground"][0]
    if case["sa_weights"]:
        assert len(results["sa_casscf"]["analysis"].states) == len(case["sa_weights"])
