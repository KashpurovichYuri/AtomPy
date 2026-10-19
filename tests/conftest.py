"""Shared fixtures: expensive reference computations once per test session."""

import numpy as np
import pytest

from atomterm.system import AtomConfig
from atomterm.drivers import (
    ReferenceSpec,
    run_rohf,
    select_shell_ao,
    select_valence_ao,
)

# Shared per-atom data (import from conftest in test modules -- pytest puts
# the tests dir on sys.path; no __init__.py by design).

# Sc CAS(3,6)/{3d,4s}, Ms=1/2 sector: the complete free-ion term table of
# 3d1 4s2 + 3d2 4s1 + 3d3 over the 90 determinants (ticket 11).
SC_EXACT_TERMS = {"2D": 20, "4F": 14, "2F": 14, "4P": 6, "2G": 18,
                  "2P": 6, "2S": 1, "2H": 11}

# Ti CAS(4,6)/{3d,4s}, Ms=1 sector: the complete term table of 3d2 4s2 +
# 3d3 4s1 + 3d4 over the 120 determinants (ticket 12). S >= 1 only (the
# sector has no singlet component); dimension check: d2s2 gives 10 =
# 3F x7 + 3P x3, d3s1 gives 60, d4 gives 50.
TI_EXACT_TERMS = {"3D": 15, "3F": 35, "3G": 18, "3H": 22, "3P": 15,
                  "5D": 5, "5F": 7, "5P": 3}


@pytest.fixture(scope="session")
def carbon(tmp_path_factory):
    """(cfg, mol, mf, act_idx) for C/cc-pVDZ triplet ROHF."""
    cfg = AtomConfig()
    chk = str(tmp_path_factory.mktemp("carbon") / "rohf.chk")
    mol, mf = run_rohf(cfg, chk)
    return cfg, mol, mf, select_shell_ao(mol, "2p")


@pytest.fixture(scope="session")
def sc_rohf(tmp_path_factory):
    """(cfg, mol, mf, act_idx) for Sc/cc-pVDZ doublet ROHF: the TM
    reference (3d1 4s2 occupation, canonical [:9] = [Ar])."""
    cfg = AtomConfig("Sc", "ccpvdz", spin=1)
    chk = str(tmp_path_factory.mktemp("sc_rohf") / "rohf.chk")
    mol, mf = run_rohf(cfg, chk)
    return cfg, mol, mf, select_valence_ao(mol, ["3d", "4s"])


# Ti occupations are forced by character (fact 21 / ticket 12): the default
# energy-ordered ROHF is an aufbau runaway (+123.07 mEh). [Ar]+4s closed,
# 3d open. occ_e_ref gates the SCF basin: with the same forced occupations
# a metastable +11.2 mEh member-gauge basin exists whose polarized core
# smears the input window at the 1e-2 level.
TI_REF = ReferenceSpec(
    occ_closed=("1s", "2s", "3s", "4s", "2p", "3p"),
    occ_open=("3d",),
    occ_e_ref=-848.40648875,
)


@pytest.fixture(scope="session")
def ti_rohf(tmp_path_factory):
    """(cfg, mol, mf, act_idx) for Ti/cc-pVDZ triplet ROHF: the
    character-forced 3d2 4s2 reference (see TI_REF)."""
    cfg = AtomConfig("Ti", "ccpvdz", spin=2)
    chk = str(tmp_path_factory.mktemp("ti_rohf") / "rohf.chk")
    mol, mf = run_rohf(cfg, chk, ref=TI_REF)
    return cfg, mol, mf, select_valence_ao(mol, ["3d", "4s"])


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
