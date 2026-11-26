"""constrained_casscf's weight-list contract, at the PySCF boundary:
a ONE-root SA list is the single root -- the measured-block rule can
hand down a single weight (a window whose ground block is one root,
e.g. an occupation-group NO window in a maximal-Ms sector) and the
SA wrapper's ensemble einsum has no axis for the scalar base return.
An EMPTY list is refused at this seam (a caller bug; the plain
kernel is ``weights=None``)."""

import pytest

from atomterm.cases import CAMPAIGN
from atomterm.casscf import constrained_casscf
from atomterm.reference import run_rohf

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def carbon_mf(tmp_path_factory):
    cfg = CAMPAIGN["C"].atom_config()
    root = tmp_path_factory.mktemp("c_casscf")
    _, mf = run_rohf(cfg, str(root / "rohf.chk"))
    return mf


def test_one_root_weights_is_the_single_root(carbon_mf, tmp_path):
    """weights=[w] takes the single-root path: the same energy as
    weights=None, one root reported (the length-1 SA list collapses
    to the plain kernel by contract).  Same driver, same start --
    observed: |dE| = 1.4e-14 (solver noise)."""
    plain = constrained_casscf(carbon_mf, 5, (1, 1), str(tmp_path / "a.chk"))
    one = constrained_casscf(carbon_mf, 5, (1, 1), str(tmp_path / "b.chk"),
                             weights=[1.0])
    assert one.energies[0] == pytest.approx(plain.energies[0], abs=1e-10)
    assert len(one.energies) == 1


def test_empty_weights_refused(carbon_mf, tmp_path):
    """weights=[] raises here, not inside PySCF's SA wrapper (R5:
    the caller bug fails fast at the seam that owns the contract)."""
    with pytest.raises(ValueError, match="empty weights"):
        constrained_casscf(carbon_mf, 5, (1, 1), str(tmp_path / "c.chk"),
                           weights=[])


@pytest.mark.fast
def test_cycle_trace_counts_distinct_imacro():
    """The macro-cycle count = the DISTINCT imacro values.  PySCF
    fires the callback once per micro iteration AND once per macro
    cycle (both carry the 1-based imacro), so an invocation count
    overcounts by the micro overhead; the last callback is the final
    macro-cycle-end call, so dE is the last macro step's true delta."""
    from types import SimpleNamespace

    from atomterm.casscf import _cycle_trace

    mc = SimpleNamespace(callback=None)
    cell = _cycle_trace(mc)
    for env in ({"imacro": 1, "e_tot": -1.0},
                {"imacro": 1, "e_tot": -1.0001},
                {"imacro": 2, "e_tot": -2.0},
                {"imacro": 2, "e_tot": -2.0001},
                {"imacro": 3, "e_tot": -3.0}):
        mc.callback(env)
    assert cell[0] == 3                       # 3 macro cycles, not 5 calls
    assert cell[1] == pytest.approx(0.9999)   # the last macro step's |dE|
