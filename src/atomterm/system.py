"""System configuration: the physical system under study (a single atom)."""

from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np

from pyscf import gto


@dataclass(frozen=True)
class AtomConfig:
    """Configuration of a single-atom calculation.

    Attributes
    ----------
    symbol : atomic symbol, e.g. ``"C"``.
    basis : PySCF basis set name (must resolve to *spherical* functions;
        cc-pVnZ family qualifies).
    spin : 2S, the spin multiplicity minus 1 (PySCF convention). Carbon
        ground state: spin=2 (triplet, Ms=1 for ROHF).
    charge : integer charge.
    coord : nuclear position in Bohr. Must be the origin: the angular
        momentum operators are defined about r=0 (see atomterm.angular).
    verbose : PySCF output level.
    """

    symbol: str = "C"
    basis: str = "ccpvdz"
    spin: int = 2
    charge: int = 0
    coord: tuple[float, float, float] = (0.0, 0.0, 0.0)
    verbose: int = 0
    extra: dict = field(default_factory=dict)

    @property
    def label(self) -> str:
        return f"{self.symbol}_{self.basis}_s{self.spin}"

    def gto_string(self) -> str:
        x, y, z = self.coord
        return f"{self.symbol} {x} {y} {z}"


def sector_occupations(nelec_act: int, spin: int) -> tuple[int, int]:
    """(na, nb) of the primary Ms sector from the active-electron count
    and PySCF's spin = 2S: ((n + 2S)/2, (n - 2S)/2)."""
    return (nelec_act + spin) // 2, (nelec_act - spin) // 2


def ncore_of(nelectron: int, nelec_act: int) -> int:
    """Doubly-occupied core count: (nelectron - nelec_act) // 2.

    Pass ``mol.nelectron`` — never recompute it from nuclear charge: an
    ECP removes core electrons, so Z - charge overcounts (the shared
    definition keeps pipeline and ladder honest under future ECPs)."""
    return (nelectron - nelec_act) // 2


def build_mol(cfg: AtomConfig) -> gto.Mole:
    """Build and validate a PySCF :class:`Mole` for the configured atom."""
    mol = gto.Mole()
    mol.atom = cfg.gto_string()
    mol.basis = cfg.basis
    mol.charge = cfg.charge
    mol.spin = cfg.spin
    mol.verbose = cfg.verbose
    mol.unit = "Bohr"
    mol.build()

    if mol.cart:
        raise ValueError(
            "L-quantum-number analysis requires spherical (pure) basis "
            "functions; got a cartesian basis."
        )
    if mol.natm != 1:
        raise ValueError("AtomTerm analyses a single atom.")
    if not np.allclose(mol.atom_coords()[0], 0.0, atol=1e-10):
        raise ValueError(
            "The atom must sit at the origin: Lz / L2 operators are defined "
            "about r = 0."
        )
    return mol
