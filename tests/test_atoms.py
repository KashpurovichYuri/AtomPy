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
    # ticket 07: same p^2 structure as carbon, one row further down the
    # periodic table -- first p-type core ([Ne] incl. 2p^6) under the
    # spherically averaged folded density. CASSCF/SA disabled: the orbital
    # optimizer escapes the 3p shell (3s mixing pulls it into an ³S-like
    # well -- the CASSCF multi-solution class, see AGENTS.md fact 12).
    # CASCI stages carry the physics assertions and are exact.
    "Si": dict(
        cfg=AtomConfig("Si", "ccpvdz", spin=2),
        patterns=["3p"],
        nelec_act=2,
        run_casscf=False,
        sa_weights=None,
        ground=("3P", 1, 1),
        primary_tag="casci_ms1",
        primary={"3P": 3},
        ms0_tag="casci_ms0",
        ms0={"1S": 1, "1D": 5, "3P": 3},
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

    # optimized-orbital stages keep the same ground term (their residual
    # symmetry breaking is physical output, asserted in test_pipeline)
    if case["run_casscf"]:
        assert results["casscf"]["analysis"].states[0].term == case["ground"][0]
    if case["sa_weights"]:
        assert len(results["sa_casscf"]["analysis"].states) == len(case["sa_weights"])
