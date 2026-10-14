"""Tests of atomterm.system: configuration validation and labeling."""

import numpy as np
import pytest

from atomterm.system import AtomConfig, build_mol


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
