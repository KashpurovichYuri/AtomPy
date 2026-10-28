"""Tests of atomterm.system: configuration validation and labeling."""

import pytest

from atomterm.system import AtomConfig, build_mol, ncore_of, sector_na_nb

pytestmark = pytest.mark.fast


def test_atom_config_label_encodes_basis_and_spin():
    assert AtomConfig().label == "C_ccpvdz_s2"
    assert AtomConfig("Si", "ccpvdz", spin=0).label == "Si_ccpvdz_s0"


def test_label_includes_charge_when_nonzero():
    assert AtomConfig("C", "ccpvdz", spin=2, charge=1).label == "C_ccpvdz_q1_s2"


def test_build_mol_yields_spherical_atom():
    """The validated product: one atom, spherical basis, requested 2S."""
    mol = build_mol(AtomConfig())
    assert mol.natm == 1
    assert not mol.cart
    assert mol.spin == 2


def test_sector_na_nb_primary_sector_and_parity_guard():
    assert sector_na_nb(2, 2) == (2, 0)
    assert sector_na_nb(3, 1) == (2, 1)
    with pytest.raises(ValueError, match="no Ms sector"):
        sector_na_nb(3, 2)  # 3 active electrons cannot carry 2S = 2


def test_ncore_of_closed_core_and_parity_guard():
    assert ncore_of(6, 2) == 2
    with pytest.raises(ValueError, match="core not closed"):
        ncore_of(9, 4)  # an odd electron cannot sit in a closed core


def test_atom_config_raw_basis_passthrough_and_tag():
    """Custom (dict) basis specs go straight to mol.basis (PySCF parses
    raw input); ``tag`` keeps the label path-safe for dict bases. Raw
    ECP-string passthrough is exercised on the real U block in the
    campaign tests (a fake minimal ECP block does not parse)."""
    from pyscf import gto

    basis = {"H": gto.parse("H S\n 1.0 1.0")}  # minimal raw spec
    cfg = AtomConfig(symbol="H", basis=basis, spin=0, charge=1,
                     tag="custom_pairing")
    mol = build_mol(cfg)
    assert cfg.label == "custom_pairing"
    assert mol.nelectron == 0  # charge=1 on one-electron H
    # a multi-line ECP string is not path-safe in the fallback label --
    # custom bases want tag (documented in the AtomConfig docstring)
    assert "\n" not in AtomConfig(
        symbol="H", basis=basis, ecp="ECP\nH nelec 0\nEND").label
