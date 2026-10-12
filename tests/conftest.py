"""Shared fixtures: expensive carbon ROHF computed once per test session."""

import numpy as np
import pytest

from atomterm.system import AtomConfig
from atomterm.drivers import run_rohf, select_shell_ao


@pytest.fixture(scope="session")
def carbon(tmp_path_factory):
    """(cfg, mol, mf, act_idx) for C/cc-pVDZ triplet ROHF."""
    cfg = AtomConfig()
    chk = str(tmp_path_factory.mktemp("carbon") / "rohf.chk")
    mol, mf = run_rohf(cfg, chk)
    return cfg, mol, mf, select_shell_ao(mol, "2p")


@pytest.fixture(scope="session")
def p_shell_mol():
    """Minimal single-center atom with an s and a p shell."""
    from pyscf import gto

    return gto.M(
        atom="C 0 0 0",
        basis=[[0, [3.0, 1.0]], [1, [1.0, 1.0]]],
        spin=2,
        verbose=0,
    )


@pytest.fixture(scope="session")
def synthetic_p_operators():
    """Exact one-electron Lx/Ly/Lz in the (px, py, pz) real-harmonic basis.

    Condon-Shortley convention: Lz|px> = i|py>, Lz|py> = -i|px>,
    Lx|py> = i|pz>, Ly|px> = -i|pz>, etc.
    """
    i = 1j
    lx = np.array([[0, 0, 0], [0, 0, -i], [0, i, 0]])
    ly = np.array([[0, 0, i], [0, 0, 0], [-i, 0, 0]])
    lz = np.array([[0, -i, 0], [i, 0, 0], [0, 0, 0]])
    return {"x": lx, "y": ly, "z": lz}
