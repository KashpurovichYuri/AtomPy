"""System configuration: the physical system under study (a single atom)."""

from __future__ import annotations

from dataclasses import dataclass

from pyscf import gto

__all__ = ["AtomConfig", "build_mol", "sector_na_nb", "ncore_of"]


@dataclass(frozen=True)
class AtomConfig:
    """Configuration of one single-atom calculation.

    Attributes
    ----------
    symbol : str
        Atomic symbol, e.g. ``"C"``. The defaults are the campaign's
        reference substrate (C / cc-pVDZ / triplet).
    basis : str | dict
        Basis-set name or raw PySCF spec (dict / gto string). Must be
        spherical -- cartesian functions carry no m quantum numbers
        (build_mol rejects them). A raw spec has no name: set ``tag``
        so artifact paths stay readable.
    spin : int
        2S = Nalpha - Nbeta (PySCF convention): 0 singlet, 1 doublet,
        2 triplet, 3 quartet. A doublet atom is spin=1, not 0.5.
    charge : int
        Integer charge.
    ecp : str | None
        Pseudopotential name or raw ECP block; None = all-electron.
        Pairs with a valence basis ("ccecp" / "ccecpccpvdz"). Under an
        ECP mol.nelectron counts valence electrons only.
    tag : str | None
        Artifact stem (chk/report names) overriding the derived label
        -- set it when basis or ecp is a raw spec with no path-safe
        name.
    verbose : int
        PySCF output level; 0 = silent.
    """

    symbol: str = "C"
    basis: str | dict = "ccpvdz"
    spin: int = 2
    charge: int = 0
    ecp: str | None = None
    tag: str | None = None
    verbose: int = 0

    @property
    def label(self) -> str:
        """Artifact stem: tag, else ``symbol_basis[_ecp][_qC]_s{spin}``."""
        if self.tag is not None:
            return self.tag
        # raw (multi-line) ECP strings are not path-safe: use tag
        ecp = f"_{self.ecp}" if self.ecp and "\n" not in self.ecp else ""
        q = f"_q{self.charge}" if self.charge else ""
        return f"{self.symbol}_{self.basis}{ecp}{q}_s{self.spin}"


def sector_na_nb(nelec_act: int, spin: int) -> tuple[int, int]:
    """(na, nb) occupations of the primary Ms sector.

    Parameters
    ----------
    nelec_act : int
        Active-window electron count.
    spin : int
        PySCF 2S = Nalpha - Nbeta.

    Returns
    -------
    tuple[int, int]
        ``((nelec_act + spin) // 2, (nelec_act - spin) // 2)``.

    Raises
    ------
    ValueError
        If ``nelec_act < spin`` or ``nelec_act - spin`` is odd: no
        such sector exists.

    Examples
    --------
    >>> sector_na_nb(2, 2)
    (2, 0)
    """
    if nelec_act < spin:
        raise ValueError(
            f"no Ms sector: spin {spin} exceeds {nelec_act} active electrons")
    if (nelec_act - spin) % 2:
        raise ValueError(
            f"no Ms sector: {nelec_act} active electrons with spin {spin}")
    return (nelec_act + spin) // 2, (nelec_act - spin) // 2


def ncore_of(nelectron: int, nelec_act: int) -> int:
    """Number of doubly-occupied core orbitals.

    Parameters
    ----------
    nelectron : int
        Pass ``mol.nelectron``, never ``Z - charge``: under an ECP it
        counts valence electrons only.
    nelec_act : int
        Active-window electron count.

    Returns
    -------
    int
        ``(nelectron - nelec_act) // 2``.

    Raises
    ------
    ValueError
        If ``nelec_act > nelectron`` or the difference is odd: the
        core would not be closed.
    """
    if nelec_act > nelectron:
        raise ValueError(
            f"core not closed: {nelec_act} active electrons exceed "
            f"{nelectron} total")
    if (nelectron - nelec_act) % 2:
        raise ValueError(
            f"core not closed: {nelectron} - {nelec_act} is odd")
    return (nelectron - nelec_act) // 2


def build_mol(cfg: AtomConfig) -> gto.Mole:
    """Build and validate a PySCF :class:`Mole` for the configured atom.

    The atom sits at the origin by construction: the angular momentum
    operators are defined about r = 0.

    Raises
    ------
    ValueError
        On a cartesian basis (cartesian functions carry no m quantum
        numbers) or a symbol that parses to more than one atom.
    """
    mol = gto.Mole()
    mol.atom = f"{cfg.symbol} 0.0 0.0 0.0"
    mol.basis = cfg.basis
    mol.charge = cfg.charge
    mol.spin = cfg.spin
    mol.ecp = cfg.ecp
    mol.verbose = cfg.verbose
    mol.build()

    if mol.cart:
        raise ValueError(
            "L-quantum-number analysis requires spherical (pure) basis "
            "functions; got a cartesian basis."
        )
    if mol.natm != 1:
        raise ValueError("AtomTerm analyses a single atom.")
    return mol
