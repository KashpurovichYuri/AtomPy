"""The cell layer: one key type addressing every calculation, its store
contract, and the CI ordering bridge between PySCF and the determinant
basis.

A *cell* is THE unit of computation and analysis (one grid point of the
manifest): the run, its checkpoint, its payload, and its table row are
all addressed by one :class:`CellKey`, serialized as the cell
*address* -- ``ti_3d4s_3-1_sa_no_free_pur_avg``. Axes resolve through
their OWNING modules (:mod:`atomterm.windows`, :mod:`atomterm.core`,
:mod:`atomterm.casscf`, :mod:`atomterm.stages`) -- cells never
re-implement physics.

Axis semantics (two of them describe the CORE and are easily conflated):

* ``substrate`` -- the core SPAN (which core columns exist): ``pol`` =
  the SCF's canonical core MOs; ``pur`` = the exactly L-closed
  reconstruction (:func:`atomterm.core.purified_core_orbitals`).
* ``fold`` -- the core FIELD convention of the analysis: ``raw`` = the
  SCF density as converged (the anisotropy detector), ``avg`` = its
  spherical average (the ruler), ``frac`` = a different reference SCF
  entirely (fractional occupations).

  Composition: (raw, pol) detector; (avg, pol) ruler; (avg, pur) the
  certified fold; (frac, pol) own reference. On ``pur`` the density is
  exactly m-scalar (fact 20), so ``raw`` degenerates with ``avg`` --
  the key CANONICALIZES ``fold="raw"`` to ``"avg"`` there, and
  ``frac`` on ``pur`` is rejected. ``substrate`` is the span, ``fold``
  is the field.

Address grammar (manifest reading order: where / the physics box / the
objective / the optimizer freedom / the core)::

    <atom>_<window>_<na>-<nb>_<method>_<start>_<mask|->_<substrate>_<fold>

with SINGULAR axis values throughout (``ao``/``mo``/``no``), the
sector as ``na-nb`` (fixes Ms AND nelec), and ``-`` for no mask
(free optimization; a mask exists only for optimizer methods --
``virtmix`` is a MASK -- :func:`atomterm.casscf.no_core_window_exchange_mask
-- not a start).

The store contract (one chkfile per cell,
``elements/<Atom>/rungs/<address>.chk``): everything ``load_cell`` needs
under ``atomterm/<address>/`` -- ``energies``, ``ci``, ``core_energy``,
``mo_coeff`` (SA runs: per-root energies under ``<address>_roots``).
The PySCF ``mcscf`` group is NOT load-bearing: this PySCF's CASCI
kernel never dumps a chkfile, so the cell layer owns a complete,
method-uniform store.

The CI ordering bridge: PySCF FCI vectors are indexed by
``(alpha-string address, beta-string address)`` with strings in
INTEGER-ASCENDING bit order, while :func:`atomterm.detbasis.determinant_basis`
enumerates lexicographic combinations -- these coincide only for
trivial sectors (norb=3, or nelec<=1 per block) and DIFFER everywhere
else (first at norb=4, na=2). Every consumer that indexes or contracts
a PySCF CI vector against detbasis-order objects must map through
:func:`pyscf_ci_permutation` first; the bridge is certified by the
eigen-equation test (generic integrals admit no permutation
automorphism).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from atomterm.detbasis import determinant_basis

__all__ = [
    "CellKey",
    "CellData",
    "cell_chk",
    "load_cell",
    "pyscf_ci_permutation",
    "ci_in_detbasis_order",
]

METHODS = ("casci", "casscf", "sa")
STARTS = ("ao", "mo", "no")       # window constructions (the span is the window)
MASKS = ("radial", "frozen", "virtmix")   # None = free optimization
SUBSTRATES = ("pol", "pur")
FOLDS = ("raw", "avg", "frac")

# the nominal window label: shell tokens like "2p", "3d4s" -- lowercase
# alphanumerics, no separators (the address splits on "_")
_WINDOW_RE = re.compile(r"[0-9][0-9a-z]*\Z")


@dataclass(frozen=True)
class CellKey:
    """One grid point of the manifest: the unit of run, load, and report.

    Every axis has an owner; the key only names the choice:

    atom : element label, e.g. ``"Ti"`` (the ``elements/<Atom>/`` dir).
    window : NOMINAL span label, e.g. ``"3d4s"`` -- names the table row;
        the construction is the ``start`` axis (window = the span,
        start = the columns defining it).
    sector : ``(na, nb)`` -- alpha/beta active electrons; fixes Ms and
        nelec together.
    method : ``"casci"`` | ``"casscf"`` | ``"sa"``.
    start : ``"ao"`` (shell AOs) | ``"mo"`` (canonical SCF slice) |
        ``"no"`` (natural-orbital groups; needs the parent-prefix
        policy, lands with the run queue).
    mask : optimizer restriction -- ``None`` (free) | ``"radial"`` |
        ``"frozen"`` | ``"virtmix"``; only meaningful for optimizer
        methods (rejected for ``casci``).
    substrate / fold : the core axes -- see the module docstring.
    """

    atom: str
    window: str
    sector: tuple[int, int]
    method: str = "casci"
    start: str = "ao"
    mask: str | None = None
    substrate: str = "pol"
    fold: str = "avg"

    def __post_init__(self):
        if not (isinstance(self.atom, str) and self.atom.isalpha()
                and len(self.atom) <= 2):
            raise ValueError(f"atom {self.atom!r} must be an element label")
        object.__setattr__(self, "atom", self.atom.capitalize())
        if not (isinstance(self.window, str) and _WINDOW_RE.match(self.window)):
            raise ValueError(
                f"window {self.window!r} must be a nominal shell label "
                "(lowercase alphanumerics, e.g. '3d4s'; no separators)")
        na, nb = self.sector
        if not (isinstance(na, int) and isinstance(nb, int)
                and na >= nb >= 0 and na + nb > 0):
            raise ValueError(f"sector {self.sector!r} must be (na >= nb >= 0, n > 0)")
        object.__setattr__(self, "sector", (na, nb))
        for axis, value, vocab in (
                ("method", self.method, METHODS),
                ("start", self.start, STARTS),
                ("substrate", self.substrate, SUBSTRATES),
                ("fold", self.fold, FOLDS)):
            if value not in vocab:
                raise ValueError(
                    f"{axis} {value!r} not in {vocab}")
        if self.mask is not None:
            if self.mask not in MASKS:
                raise ValueError(f"mask {self.mask!r} not in {MASKS} (None = free)")
            if self.method == "casci":
                raise ValueError("a mask needs an optimizer method (casscf/sa)")
        # the pur-fold degeneracy rule: on the purified substrate the
        # core density is exactly m-scalar, so raw == avg (fact 20) --
        # canonicalize instead of storing two addresses for one cell;
        # frac is a different REFERENCE, meaningless on pur
        if self.substrate == "pur":
            if self.fold == "frac":
                raise ValueError("fold='frac' is a different reference SCF; "
                                 "it cannot live on the purified substrate")
            if self.fold == "raw":
                object.__setattr__(self, "fold", "avg")

    @property
    def address(self) -> str:
        """The derived serialization: chk name, payload key, table-row id."""
        na, nb = self.sector
        return "_".join((self.atom.lower(), self.window, f"{na}-{nb}",
                         self.method, self.start, self.mask or "-",
                         self.substrate, self.fold))

    @classmethod
    def from_address(cls, address: str) -> "CellKey":
        """Parse an address (inverse of :attr:`address`)."""
        parts = address.split("_")
        if len(parts) != 8:
            raise ValueError(
                f"address {address!r} does not have 8 fields: "
                "atom_window_na-nb_method_start_mask_substrate_fold")
        atom, window, sector, method, start, mask, substrate, fold = parts
        na, nb = sector.split("-")
        return cls(atom, window, (int(na), int(nb)), method, start,
                   None if mask == "-" else mask, substrate, fold)

    @property
    def shell_patterns(self) -> list[str]:
        """The window label as shell selectors: ``"3d4s"`` -> ``["3d", "4s"]``."""
        return re.findall(r"\d+[spdfgh]", self.window)


@dataclass
class CellData:
    """What a stored cell holds -- pure data, no live mol/mf needed.

    ``ci`` stays in PYSCF FCI ordering over the active orbitals
    (norb = ncas): map through :func:`ci_in_detbasis_order` before any
    detbasis-order use (indexing, contracting with our operator
    matrices, RDMs).
    """

    key: CellKey
    mo_coeff: np.ndarray          # the run's full MO matrix (nao, nao)
    energies: np.ndarray          # total energies per root
    ci: list                      # per-root arrays, pyscf ordering
    core_energy: float            # E_tot - active part (the fold constant)


def cell_chk(key: CellKey, workdir: str | Path = "elements") -> Path:
    """The chkfile of a cell: ``<workdir>/<Atom>/rungs/<address>.chk``."""
    return Path(workdir) / key.atom / "rungs" / f"{key.address}.chk"


def _read_store(chk: Path, name: str) -> tuple | None:
    """Read one ``atomterm/<name>/`` store: (energies, ci, core_energy).

    ``None`` = store absent (caller tries the next name).
    """
    from pyscf import lib

    e = lib.chkfile.load(str(chk), f"atomterm/{name}/energies")
    if e is None:
        return None
    ci = lib.chkfile.load(str(chk), f"atomterm/{name}/ci")
    core = lib.chkfile.load(str(chk), f"atomterm/{name}/core_energy")
    if ci is None or core is None:
        raise KeyError(
            f"cell store {chk}::{name} is incomplete (energies without "
            "ci/core_energy); regenerate the cell")
    return np.asarray(e), list(np.asarray(ci, dtype=object)), float(core)


def load_cell(key: CellKey, workdir: str | Path = "elements") -> CellData:
    """Load a stored cell -- pure data, no live mol/mf required.

    The re-analysis entry point (:func:`analyze_cell`) builds on this;
    the point is re-reading a run (post-hoc properties, tables) without
    re-running the optimizer.
    """
    chk = cell_chk(key, workdir)
    if not chk.exists():
        raise FileNotFoundError(
            f"no stored cell {key.address} at {chk}; run it first "
            "(atomterm.cells.run_cell)")
    store = None
    for name in (f"{key.address}_roots", key.address):  # SA: roots pass wins
        store = _read_store(chk, name)
        if store is not None:
            break
    if store is None:
        raise KeyError(f"no atomterm/{key.address}* store in {chk}")
    from pyscf import lib

    mo = lib.chkfile.load(str(chk), f"atomterm/{key.address}/mo_coeff")
    if mo is None:
        raise KeyError(
            f"{chk} has no atomterm/{key.address}/mo_coeff (a pre-cell "
            "chkfile); regenerate the cell")
    energies, ci, core_energy = store
    return CellData(key, np.asarray(mo), energies, ci, core_energy)


# ---------------------------------------------------------------------
# The CI ordering bridge: pyscf FCI ordering <-> detbasis ordering.
# ---------------------------------------------------------------------


def pyscf_ci_permutation(norb: int, na: int, nb: int) -> np.ndarray:
    """Map ravelled PySCF CI amplitudes onto the detbasis ordering.

    PySCF's direct_spin1 vector ``v`` (shape ``(nalpha_strings,
    nbeta_strings)``, ravelled alpha-major) indexes determinants by
    string ADDRESSES (integer-ascending bit strings); detbasis
    enumerates lexicographic combinations. These orders differ for
    norb >= 4 once a block holds >= 2 electrons -- coinciding only in
    trivial sectors (carbon's norb=3 window is why the difference went
    unnoticed).

    Returns ``perm`` with ``v[perm[k]] == w_detbasis[k]``-semantics:
    the pyscf amplitude at ravelled index ``k`` belongs to determinant
    ``determinant_basis(norb, na, nb)[perm[k]]``. Scatter with
    ``w[perm] = v`` (or use :func:`ci_in_detbasis_order`).
    """
    from pyscf.fci import cistring

    dets = determinant_basis(norb, na, nb)
    index = {d: i for i, d in enumerate(dets)}
    a_strs = cistring.make_strings(range(norb), na)  # index == address
    b_strs = cistring.make_strings(range(norb), nb)
    perm = np.empty(len(dets), dtype=np.int64)
    k = 0
    for sa in a_strs:
        occ_a = [i for i in range(norb) if sa >> i & 1]  # orbital i = bit i
        for sb in b_strs:
            occ_b = [i for i in range(norb) if sb >> i & 1]
            det = tuple(sorted(occ_a)) + tuple(o + norb for o in occ_b)
            perm[k] = index[det]
            k += 1
    return perm


def ci_in_detbasis_order(v: np.ndarray, norb: int, na: int, nb: int) -> np.ndarray:
    """One ravelled PySCF CI vector, reordered to detbasis indexing."""
    v = np.asarray(v).ravel()
    w = np.empty_like(v)
    w[pyscf_ci_permutation(norb, na, nb)] = v
    return w
