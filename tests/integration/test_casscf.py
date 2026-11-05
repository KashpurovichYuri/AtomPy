"""constrained_casscf's weight-list contract, at the PySCF boundary:
a ONE-root SA list is the single root -- the measured-block rule can
hand down a single weight (a window whose ground block is one root,
e.g. an occupation-group NO window in a maximal-Ms sector) and the
SA wrapper's ensemble einsum has no axis for the scalar base return."""

import pytest

from atomterm.cases import CAMPAIGN
from atomterm.casscf import constrained_casscf
from atomterm.reference import run_rohf


@pytest.fixture(scope="module")
def carbon_mf(tmp_path_factory):
    cfg = CAMPAIGN["C"].atom_config()
    root = tmp_path_factory.mktemp("c_casscf")
    _, mf = run_rohf(cfg, str(root / "rohf.chk"))
    return mf


def test_one_root_weights_is_the_single_root(carbon_mf, tmp_path):
    """weights=[w] takes the single-root path: the same energy as
    weights=None, one root reported (the length-1 SA list collapses
    to the plain kernel by contract)."""
    plain = constrained_casscf(carbon_mf, 5, (1, 1), str(tmp_path / "a.chk"))
    one = constrained_casscf(carbon_mf, 5, (1, 1), str(tmp_path / "b.chk"),
                             weights=[1.0])
    assert one.energies[0] == pytest.approx(plain.energies[0], abs=1e-9)
    assert len(one.energies) == 1
