"""Tests of atomterm.system: configuration validation and labeling."""

import numpy as np
import pytest

from atomterm.system import AtomConfig, build_mol

pytestmark = pytest.mark.fast


def test_atom_config_label_encodes_basis_and_spin():
    assert AtomConfig().label == "C_ccpvdz_s2"
    assert AtomConfig("Si", "ccpvdz", spin=0).label == "Si_ccpvdz_s0"


def test_build_mol_rejects_off_origin_atom():
    """Lz / L^2 are defined about r = 0 -- the constraint is load-bearing."""
    with pytest.raises(ValueError, match="origin"):
        build_mol(AtomConfig(coord=(0.0, 0.0, 0.1)))


def test_build_mol_yields_origin_centered_spherical_atom():
    """The validated product: one atom, spherical basis, requested 2S."""
    mol = build_mol(AtomConfig())
    assert mol.natm == 1
    assert not mol.cart
    assert mol.spin == 2
    assert np.allclose(mol.atom_coords()[0], 0.0, atol=1e-10)


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
