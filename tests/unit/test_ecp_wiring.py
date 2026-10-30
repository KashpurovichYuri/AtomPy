"""Unit tier of the ECP mode (ticket 18): ccECP pairing and the
valence-electron bookkeeping -- mol-only, no SCF."""

import pytest

from atomterm.system import AtomConfig, build_mol, ncore_of, sector_na_nb
from atomterm.windows import select_shell_ao

pytestmark = pytest.mark.fast

C_ECP = AtomConfig(symbol="C", basis="ccecpccpvdz", ecp="ccecp", spin=2)
SI_ECP = AtomConfig(symbol="Si", basis="ccecpccpvdz", ecp="ccecp", spin=2)
SC_ECP = AtomConfig(symbol="Sc", basis="ccecpccpvdz", ecp="ccecp", spin=1)


def test_ecp_wiring_valence_bookkeeping():
    """mol.nelectron is ECP-honest (valence only) and the sector geometry
    derives from it: C/Si -> 4 e (ncore=1, the explicit ns^2), Sc -> 11 e
    (ncore=4, the explicit 3s^2 3p^6)."""
    for cfg, pat, nelec, ncore in [(C_ECP, "2p", 4, 1), (SI_ECP, "3p", 4, 1)]:
        mol = build_mol(cfg)
        assert mol.nelectron == nelec
        assert not mol.cart
        assert len(select_shell_ao(mol, pat)) == 3  # complete m-set
        assert ncore_of(mol.nelectron, 2) == ncore
        assert sector_na_nb(2, cfg.spin) == (2, 0)
    mol = build_mol(SC_ECP)
    assert mol.nelectron == 11  # Ne-core ECP: 3s 3p 3d 4s explicit
    assert ncore_of(mol.nelectron, 3) == 4
