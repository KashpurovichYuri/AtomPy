"""Integration-tier fixtures: the campaign references (Sc/Ti gated
ROHF). The pinned term tables and ReferenceSpec constants live in
``refs.py`` (a plain module -- imports of the ``conftest`` basename
are order-dependent across tiers and broke the runtime fixture
import); fixtures stay here."""

from atomterm.reference import ReferenceSpec, run_rohf
from atomterm.windows import select_valence_ao
import pytest

from atomterm.system import AtomConfig

from refs import TI_REF


@pytest.fixture(scope="session")
def sc_rohf(tmp_path_factory):
    """(cfg, mol, mf, act_idx) for Sc/cc-pVDZ doublet ROHF: the TM
    reference (3d1 4s2 occupation, canonical [:9] = [Ar])."""
    cfg = AtomConfig("Sc", "ccpvdz", spin=1)
    chk = str(tmp_path_factory.mktemp("sc_rohf") / "rohf.chk")
    mol, mf = run_rohf(cfg, chk)
    return cfg, mol, mf, select_valence_ao(mol, ["3d", "4s"])


@pytest.fixture(scope="session")
def ti_rohf(tmp_path_factory):
    """(cfg, mol, mf, act_idx) for Ti/cc-pVDZ triplet ROHF: the
    character-forced 3d2 4s2 reference (see TI_REF)."""
    cfg = AtomConfig("Ti", "ccpvdz", spin=2)
    chk = str(tmp_path_factory.mktemp("ti_rohf") / "rohf.chk")
    mol, mf = run_rohf(cfg, chk, ref=TI_REF)
    return cfg, mol, mf, select_valence_ao(mol, ["3d", "4s"])
