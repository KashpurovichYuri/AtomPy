"""Shared fixtures: the unit tier's small references (one carbon ROHF
by design -- "we can build C even within unit") and the synthetic
no-SCF substrates. Campaign fixtures (Sc/Ti gated ROHF, term tables)
live in integration/conftest.py."""

from atomterm.windows import select_shell_ao
import numpy as np
import pytest

from atomterm.system import AtomConfig, build_mol


@pytest.fixture(scope="session")
def carbon():
    """(cfg, mol, mf, act_idx) for C/cc-pVDZ triplet ROHF -- the
    FROZEN reference (``tests/carbon_reference.py``): bit-stable
    across runs, NO SCF (the state loads from source; floats are
    repr-round-trip exact)."""
    from pyscf import scf

    import carbon_reference as ref

    cfg = AtomConfig()
    mol = build_mol(cfg)
    mf = scf.ROHF(mol)
    mf.e_tot = ref.E_TOT
    mf.mo_energy = ref.MO_ENERGY
    mf.mo_occ = ref.MO_OCC
    mf.mo_coeff = ref.MO_COEFF
    mf.converged = True
    return cfg, mol, mf, select_shell_ao(mol, "2p")


@pytest.fixture(scope="module")
def shell_ops(carbon):
    """Operators in the pure 2p-shell determinant spaces (Ms=1 and Ms=0)."""
    from atomterm.operators import build_operators

    _, mol, mf, act_idx = carbon
    iao = np.eye(mol.nao_nr())
    c_act = iao[:, act_idx]
    ncore = (mol.nelectron - 2) // 2  # 2 for both (2,0) and (1,1)
    ops1 = build_operators(mol, mf, c_act, ncore, na=2, nb=0)
    ops0 = build_operators(mol, mf, c_act, ncore, na=1, nb=1)
    return mol, mf, ops1, ops0


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
    import numpy as np

    i = 1j
    lx = np.array([[0, 0, 0], [0, 0, -i], [0, i, 0]])
    ly = np.array([[0, 0, i], [0, 0, 0], [-i, 0, 0]])
    lz = np.array([[0, -i, 0], [i, 0, 0], [0, 0, 0]])
    return {"x": lx, "y": ly, "z": lz}
