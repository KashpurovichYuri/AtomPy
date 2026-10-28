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
    mol = build_mol(AtomConfig())
    assert mol.natm == 1
    assert not mol.cart
    assert mol.spin == 2


def test_build_mol_rejects_multi_atom_symbol():
    with pytest.raises(ValueError, match="single atom"):
        build_mol(AtomConfig(symbol="H 0 0 0; H 0 0 1.5"))


def test_sector_na_nb_primary_sector():
    assert sector_na_nb(2, 2) == (2, 0)
    assert sector_na_nb(3, 1) == (2, 1)


def test_sector_na_nb_rejects_impossible_sectors():
    with pytest.raises(ValueError, match="exceeds"):
        sector_na_nb(2, 4)  # spin above the electron count: nb < 0
    with pytest.raises(ValueError, match="no Ms sector"):
        sector_na_nb(3, 2)  # odd difference: no such sector


def test_ncore_of_counts_closed_core_orbitals():
    assert ncore_of(6, 2) == 2
    assert ncore_of(11, 3) == 4  # Sc/ccECP: the explicit 3s2 3p6


def test_ncore_of_rejects_open_core():
    with pytest.raises(ValueError, match="exceed"):
        ncore_of(2, 6)  # more active than total electrons: ncore < 0
    with pytest.raises(ValueError, match="core not closed"):
        ncore_of(9, 4)  # an odd electron cannot sit in a closed core


def test_atom_config_raw_basis_and_tag():
    """Raw basis specs pass through to mol.basis; tag keeps the label
    path-safe (the multiline-ECP case)."""
    from pyscf import gto

    basis = {"H": gto.parse("H S\n 1.0 1.0")}  # minimal raw spec
    cfg = AtomConfig(symbol="H", basis=basis, spin=0, charge=1,
                     tag="custom_pairing")
    mol = build_mol(cfg)
    assert cfg.label == "custom_pairing"
    assert mol.nelectron == 0  # charge=1 on one-electron H
    assert "\n" not in AtomConfig(
        symbol="H", basis=basis, ecp="ECP\nH nelec 0\nEND").label
