"""The cell layer: one key type addressing every calculation, the
run/load/analyze seam over the owning modules, and the self-owned
store.

A *cell* is THE unit of computation and analysis (one grid point of the
manifest): the run, its checkpoint, its payload, and its table row are
all addressed by one :class:`CellKey`, serialized as the cell
*address* (the grammar below). Each axis resolves through its OWNING
module -- cells never re-implements physics.

Division of labor
-----------------
With the declaration layer:
:mod:`atomterm.cases` declares WHAT is studied (one spec per atom --
windows, references, sectors); this module addresses HOW each single
calculation runs -- one frozen key per combination of the axes, the
declared key-sets (:func:`window_grid` and its filters) derived from
those specs. A new axis value arrives here (vocabulary + one
resolution branch through its owner); a new atom arrives there (one
registry row).

Axis semantics
--------------
Two axes describe the CORE and are easily conflated:

* ``substrate`` -- the core SPAN (which core columns exist): ``pol`` =
  the SCF's canonical core MOs; ``pur`` = the exactly L-closed
  reconstruction (:func:`atomterm.core.purified_core_orbitals`).
* ``fold`` -- the core FIELD convention of the analysis: ``raw`` = the
  SCF density as converged (the anisotropy detector), ``avg`` = its
  spherical average (the ruler), ``frac`` = a different reference SCF
  entirely (fractional occupations).

  Composition: (raw, pol) detector; (avg, pol) ruler; (avg, pur) the
  certified fold; (frac, pol) own reference. On ``pur`` the density is
  exactly m-scalar, so ``raw`` degenerates with ``avg`` --
  the key CANONICALIZES ``fold="raw"`` to ``"avg"`` there, and
  ``frac`` on ``pur`` is rejected. ``substrate`` is the span, ``fold``
  is the field.

Address grammar
---------------
Manifest reading order: where / the physics box /
the objective / the optimizer freedom / the core::

    <atom>_<window>_<na>-<nb>_<method>_<start>_<mask|->_<substrate>_<fold>

with SINGULAR axis values throughout (``ao``/``mo``/``no``), the
sector as ``na-nb`` (fixes Ms AND nelec), and ``-`` for no mask
(free optimization; a mask exists only for optimizer methods --
``virtmix`` is a MASK -- :func:`atomterm.casscf.mask_no_core_window_mix
-- not a start). NO windows append a prefix index to the label:
``<source>no<k>`` (e.g. ``3d4sno2``) = the k-th occupation-group
prefix of the source's term-averaged density -- the prefix is part of
the window identity, so differently-sized NO cells of one source get
different addresses.

One real address, decoded slot by slot::

    ti_3d4sno2_3-1_sa-casscf_no_-_pol_avg
    │  │        │   │         │  │ │   └ fold=avg   (the averaged core field)
    │  │        │   │         │  │ └─── substrate=pol (the SCF core span)
    │  │        │   │         │  └───── mask=-      (free optimization)
    │  │        │   │         └──────── start=no    (natural-orbital columns)
    │  │        │   └────────────────── method=sa-casscf
    │  │        └────────────────────── sector=(3, 1) (na-nb; Ms and nelec)
    │  └─────────────────────────────── window=3d4sno2 (the 2-group NO
    │                                   prefix of the 3d4s source)
    └────────────────────────────────── atom=ti

The store contract
------------------
One chkfile per cell,
``elements/<Atom>/rungs/<address>.chk``; everything ``load_cell`` needs
under ``atomterm/<address>/`` -- ``energies``, ``ci_stack`` (a plain
``(nroots, ndet)`` array in PYSCF FCI ordering; PySCF's object-array
serialization does not round-trip for 2-D roots -- map through
:func:`atomterm.operators.ci_in_detbasis_order` before any
detbasis-order use), ``core_energy``, ``mo_coeff`` (SA runs: per-root
energies under ``<address>_roots``). The PySCF ``mcscf`` group is NOT
load-bearing: this PySCF's CASCI kernel never dumps a chkfile, so the
cell layer owns a complete, method-uniform store.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from atomterm.detbasis import sector_dimension

if TYPE_CHECKING:  # the declaration groups consume the intent manifest;
    from atomterm.cases import CaseSpec  # annotation-only (no cycle)

__all__ = [
    "CellKey",
    "CellData",
    "cell_chk",
    "cell_payload",
    "load_cell",
    "run_cell",
    "analyze_cell",
    "spine_keys",
    "window_grid",
    "grid_keys",
    "diagnostics",
    "planned",
    "lattice_keys",
]

METHODS = ("casci", "casscf", "sa-casscf")
STARTS = ("ao", "mo", "no")       # window constructions (the span is the window)
MASKS = ("radial", "frozen", "virtmix")   # None = free optimization
SUBSTRATES = ("pol", "pur")
FOLDS = ("raw", "avg", "frac")

# the nominal window label: shell tokens like "2p", "3d4s" -- lowercase
# alphanumerics, no separators (the address splits on "_")
_WINDOW_RE = re.compile(r"[0-9][0-9a-z]*\Z")


# the cell key -------------------------------------------------------

@dataclass(frozen=True)
class CellKey:
    """One grid point of the manifest: the unit of run, load, and report.

    Every axis has an owner; the key only names the choice.

    Attributes
    ----------
    atom : element label, e.g. ``"Ti"`` (the ``elements/<Atom>/`` dir).
    window : NOMINAL span label, e.g. ``"3d4s"`` -- names the table
        row; the construction is the ``start`` axis (window = the
        span, start = the columns defining it).
    sector : ``(na, nb)`` -- alpha/beta active electrons; fixes Ms
        and nelec together.
    method : ``"casci"`` | ``"casscf"`` | ``"sa-casscf"``.
    start : ``"ao"`` (shell AOs) | ``"mo"`` (canonical SCF slice) |
        ``"no"`` (natural-orbital groups; the label grammar is
        ``<source>no<k>``, see :attr:`no_source`).
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
        """Validate and canonicalize: axis vocabularies, void
        combinations, and the pur fold rule (raw -> avg)."""
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
        if (self.start == "no") != (self.no_source is not None):
            raise ValueError(
                f"start='no' and the window label must agree: NO windows "
                "are '<source>no<k>' (e.g. '3d4sno2'), other starts take "
                "a plain shell label")
        if self.mask is not None:
            if self.mask not in MASKS:
                raise ValueError(f"mask {self.mask!r} not in {MASKS} (None = free)")
            if self.method == "casci":
                raise ValueError("a mask needs an optimizer method "
                                 "(casscf/sa-casscf)")
            if self.mask == "radial" and self.start != "ao":
                raise ValueError(
                    "the radial mask rotates member-labeled AO slots; it "
                    "needs start='ao' (MO and NO columns carry no member "
                    "labels)")
        # (mo, pur) is valid: the canonical columns are
        # re-orthogonalized against the pur core span (run_cell's mo
        # branch).
        # the pur-fold degeneracy rule: on the purified substrate the
        # core density is exactly m-scalar, so raw == avg --
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

    @property
    def no_source(self) -> tuple[str, int] | None:
        """For NO windows: ``(source label, group prefix k)``. The label
        grammar is ``<source>no<k>`` (e.g. ``"3d4sno2"`` = the second
        occupation-group prefix of the 3d4s source's term-averaged
        density) -- shell letters cannot contain "no", so the tail is
        unambiguous. The prefix IS part of the window identity: two NO
        cells of one source differ in k and get different addresses.
        """
        m = re.fullmatch(r"(.+?)no([1-9]\d*)", self.window)
        return (m.group(1), int(m.group(2))) if m else None


@dataclass
class CellData:
    """What a stored cell holds -- pure data, no live mol/mf needed.

    ``ci`` stays in PYSCF FCI ordering over the active orbitals
    (norb = ncas): map through
    :func:`atomterm.operators.ci_in_detbasis_order` before any
    detbasis-order use (indexing, contracting with our operator
    matrices, RDMs).

    Attributes
    ----------
    key : the :class:`CellKey` this data belongs to.
    mo_coeff : the run's full MO matrix (nao, nao).
    energies : total energies per root.
    ci : per-root CI arrays, PySCF ordering.
    core_energy : ``E_tot`` minus the active part (the fold constant).
    ncas : window width -- needed exactly when it is not derivable
        from the label (NO windows); 0 = derive from the shell label.
    """

    key: CellKey
    mo_coeff: np.ndarray
    energies: np.ndarray
    ci: list
    core_energy: float
    ncas: int = 0


# the store ----------------------------------------------------------

def cell_chk(key: CellKey, workdir: str | Path = "elements") -> Path:
    """The chkfile of a cell: ``<workdir>/<Atom>/rungs/<address>.chk``."""
    return Path(workdir) / key.atom / "rungs" / f"{key.address}.chk"


def _read_store(chk: Path, name: str) -> tuple | None:
    """Read one ``atomterm/<name>/`` store: (energies, ci, core_energy).

    ``ci`` comes from the cell layer's ``ci_stack`` (nroots, ndet).
    ``None`` = store absent (caller tries the next name).
    """
    from pyscf import lib

    e = lib.chkfile.load(str(chk), f"atomterm/{name}/energies")
    if e is None:
        return None
    stack = lib.chkfile.load(str(chk), f"atomterm/{name}/ci_stack")
    core = lib.chkfile.load(str(chk), f"atomterm/{name}/core_energy")
    if stack is None or core is None:
        raise KeyError(
            f"cell store {chk}::{name} is incomplete (energies without "
            "ci_stack/core_energy); regenerate the cell")
    return np.asarray(e), list(np.asarray(stack)), float(core)


def load_cell(key: CellKey, workdir: str | Path = "elements") -> CellData:
    """Load a stored cell -- pure data, no live mol/mf required.

    The re-analysis entry point (:func:`analyze_cell`) builds on this;
    the point is re-reading a run (post-hoc properties, tables)
    without re-running the optimizer.

    Raises
    ------
    FileNotFoundError
        No chkfile for the address (run the cell first).
    KeyError
        The store is incomplete or the chkfile predates the store
        contract -- regenerate the cell.
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
            f"{chk} has no atomterm/{key.address}/mo_coeff (the "
            "chkfile predates the store contract); regenerate the cell")
    energies, ci, core_energy = store
    ncas = lib.chkfile.load(str(chk), f"atomterm/{key.address}/ncas")
    return CellData(key, np.asarray(mo), energies, ci, core_energy,
                    int(ncas) if ncas is not None else 0)


# ------------------------------------------------------------- the execution


def _pick_no_window(wins, k: int):
    """The k-GROUP-PREFIX window from no_group_windows's COMPACT list.

    k counts occupation groups (window identity: ``3d4sno2`` = the
    2-group prefix), but the returned list skips prefixes too small to
    open the sector ({3d,4s} over (3,1): the 1-column s prefix) -- so
    labels can be non-contiguous. Select by label, never by index."""
    for w in wins:
        if w.label == f"no{k}":
            return w
    raise ValueError(
        f"the source's usable group prefixes are {[w.label for w in wins]}; "
        f"no{k} is not among them (absent, or skipped as too small to "
        "open the sector)")


def run_cell(key: CellKey, mf, mol=None, *, frac_mf=None, budgets: dict = None,
             fresh: bool = False, verbose: int = 0,
             workdir: str | Path = "elements") -> "StageFacts":
    """Execute one cell: resolve every axis through its OWNER, run,
    store, analyze. The unit of the run queue.

    Parameters
    ----------
    key : the cell to run.
    mf : the atom's converged reference (REQUIRED -- reference
        construction, including character-forced occupations and the
        basin gate, belongs to cases/reference, not to cells). On
        ``fold="frac"`` this is REPLACED by ``frac_mf``.
    mol : the analysis mol; defaults to ``mf.mol``. Only a
        ``fold="frac"`` run needs it stated separately (the run basis
        is ``frac_mf.mol``).
    frac_mf : the fractional-occupation reference; required iff
        ``key.fold == "frac"`` (that fold is a different reference SCF
        -- the run AND the fold use it).
    budgets : per-run knobs that are NOT cell identity:
        ``mc_opts``/``fci_opts`` (optimizer/solver), ``sa_weights``,
        ``verbose``, ``no_cap`` (the NO-source cost stop).
        ``sa_weights`` default: uniform over the MEASURED ground block
        of the start window (:func:`atomterm.operators.ground_block`,
        the SA-anchoring rule). The default presupposes a DEGENERATE
        start window -- the block measurement needs the ground term's
        M_L set intact; on a defective window (the input columns'
        core-orthogonalization split) it anchors
        ONE state: single-root SA, which collapses on TMs. Declare
        ``sa_weights`` explicitly there.
    fresh : ignore a stored cell and re-run (a cache hit otherwise
        re-analyzes the stored run -- budgets are then ignored).
    workdir : the store root (``<workdir>/<Atom>/rungs/``).

    Returns
    -------
    The stage's :class:`~atomterm.stages.StageFacts`; persist with
    :func:`cell_payload` (which embeds the address).
    """
    budgets = dict(budgets or {})
    if key.fold == "frac":
        if frac_mf is None:
            raise ValueError("fold='frac' needs the fractional-occupation "
                             "reference (frac_mf=...)")
        mf = frac_mf
        mol = mol if mol is not None else mf.mol
    mol = mf.mol if mol is None else mol

    chk = cell_chk(key, workdir)
    if chk.exists() and not fresh:
        return analyze_cell(key, mol, mf, workdir=workdir)
    chk.parent.mkdir(parents=True, exist_ok=True)

    from atomterm.casscf import (constrained_casscf, mask_frozen_window,
                                 mask_no_core_window_mix, run_casci)
    from atomterm.core import purified_core_orbitals
    from atomterm.system import ncore_of
    from atomterm.windows import (build_full_mo, occupied_core_idx,
                                  radial_basis, mask_radial_rotations, select_valence_ao)

    na, nb = key.sector
    ncore = ncore_of(mol.nelectron, na + nb)
    act_idx = select_valence_ao(mol, key.shell_patterns)
    ncas = len(act_idx)
    c_pur = (purified_core_orbitals(mol, mf, ncore)[0]
             if key.substrate == "pur" else None)

    # start x mask -> the substrate MO matrix + the optimizer mask.
    # frozen/virtmix are START-INDEPENDENT (they pin the window columns
    # whatever substrate they came from); radial needs the
    # member-labeled AO substrate and is ao-only (validated at the key).
    mask = None
    if key.start == "no":
        # the parent-prefix policy: the label names the SOURCE window,
        # the k-th occupation-group prefix of its term-averaged density
        # is the window (never splits an m-set; exactly closed when the
        # density was m-scalar -- natural.no_group_windows verifies)
        from atomterm.natural import no_group_windows

        mo_src = build_full_mo(mf, occupied_core_idx(mf, ncore), act_idx,
                               c_core=c_pur)
        c_src = mo_src[:, ncore:ncore + ncas]
        wins = no_group_windows(mol, mf, ncore, c_src, na, nb,
                                c_core=c_pur, cap=budgets.get("no_cap", 12000))
        win = _pick_no_window(wins, key.no_source[1])
        ncas = win.c_win.shape[1]
        mo_full = win.columns(mf, ncore, c_core=c_pur)
    elif key.start == "ao":
        if key.mask == "radial":
            mo_full, labels = radial_basis(
                mf, occupied_core_idx(mf, ncore), act_idx, c_core=c_pur)
            mask = mask_radial_rotations(labels, ncore, ncas)
        else:
            mo_full = build_full_mo(mf, occupied_core_idx(mf, ncore),
                                    act_idx, c_core=c_pur)
    else:  # start == "mo": PySCF's default start (canonical MOs,
        # canonical window). On pol the canonical MO matrix IS the
        # substrate; on pur the canonical columns are re-orthogonalized
        # against the pur core span (they are S-orthonormal against
        # the core they were SOLVED in, not the pur one -- exactly the
        # input-window defect build_full_mo's projection removes at the
        # source). Honest
        # expectation: the window's own m-delocalization survives
        # (measured: raw-MO spans never rotation-close) -- these cells measure that, they do not fix it.
        if key.substrate == "pur":
            mo_full = build_full_mo(mf, occupied_core_idx(mf, ncore),
                                    mf.mo_coeff[:, ncore:], c_core=c_pur)
        else:
            mo_full = mf.mo_coeff
    if key.mask == "frozen":
        mask = mask_frozen_window
    elif key.mask == "virtmix":
        mask = mask_no_core_window_mix

    verbose = budgets.get("verbose", verbose)
    common = dict(mc_opts=budgets.get("mc_opts"),
                  fci_opts=budgets.get("fci_opts"), verbose=verbose)
    if key.method == "casci":
        ndet = sector_dimension(ncas, na, nb)
        cas = run_casci(mf, ncas, (na, nb), mo_full, str(chk),
                        nroots=min(9, ndet), name=key.address,
                        verbose=verbose, fci_opts=common["fci_opts"])
    else:
        weights = (budgets.get("sa_weights")
                   if key.method == "sa-casscf" else None)
        if key.method == "sa-casscf" and weights is None:
            from atomterm.operators import build_operators, ground_block

            c_act0 = mo_full[:, ncore:ncore + ncas]
            ops0 = build_operators(mol, mf, c_act0, ncore, na, nb,
                                   c_core=c_pur)
            _, _, weights, _ = ground_block(ops0)  # the SA-anchoring rule
            weights = list(weights)
        cas = constrained_casscf(mf, ncas, (na, nb), str(chk),
                                 mo_full=mo_full, mask=mask, weights=weights,
                                 name=key.address, **common)

    _store_cell(chk, key, cas, ncas,
                  store_name=f"{key.address}_roots"
                            if key.method == "sa-casscf"
                  else None)
    return _facts(key, mol, mf, cas, ncore, ncas, na, nb, c_pur=c_pur)


# analysis + persistence ---------------------------------------------

def analyze_cell(key: CellKey, mol, mf, workdir: str | Path = "elements"):
    """Re-analyze a stored cell (load + one facts pass) -- no mol/mf-bound
    optimizer re-run, no re-solving. The two-function lean: loading is
    pure data (:func:`load_cell`), analysis needs the live reference.

    ``converged`` defaults True (a trusted cached stage).
    """
    from atomterm.casscf import CASResult
    from atomterm.core import purified_core_orbitals
    from atomterm.system import ncore_of
    from atomterm.windows import select_valence_ao

    data = load_cell(key, workdir)
    cas = CASResult(data.energies, data.ci, data.mo_coeff,
                    None, data.core_energy)
    na, nb = key.sector
    ncore = ncore_of(mol.nelectron, na + nb)
    ncas = data.ncas or len(select_valence_ao(mol, key.shell_patterns))
    c_pur = (purified_core_orbitals(mol, mf, ncore)[0]
             if key.substrate == "pur" else None)
    return _facts(key, mol, mf, cas, ncore, ncas, na, nb, c_pur=c_pur)


def cell_payload(key: CellKey, facts) -> dict:
    """The persisted payload of a cell: the stage facts + the address."""
    out = facts.payload(key.method)
    out["cell"] = key.address
    return out


def _store_cell(chk: Path, key: CellKey, cas, ncas: int,
                store_name: str | None = None) -> None:
    """Complete the cell store (the writer half of the load_cell pair):
    energies were stored by the run; add the pieces only the cell
    layer owns. ``ci_stack`` is a plain
    (nroots, ndet) float array -- PySCF's chkfile serializes object
    arrays as nested ``__from_list__`` groups that do not round-trip
    for 2-D roots, so the cell layer owns its own CI serialization.

    ``store_name``: the kernel store the returned ``cas`` corresponds
    to. SA runs create TWO (the ensemble pass under ``<address>``, the
    per-root pass under ``<address>_roots``) and return the roots one
    -- completing ``<address>`` instead would pair the ensemble's
    energies with the roots' CI. :func:`load_cell` prefers the roots
    store for SA; ``mo_coeff``/``ncas`` live under ``<address>``."""
    from pyscf import lib

    name = store_name or key.address
    lib.chkfile.save(str(chk), f"atomterm/{name}/core_energy",
                     float(cas.core_energy))
    lib.chkfile.save(str(chk), f"atomterm/{name}/ci_stack",
                     np.asarray([np.asarray(c).ravel() for c in cas.ci]))
    lib.chkfile.save(str(chk), f"atomterm/{key.address}/mo_coeff",
                     np.asarray(cas.mo_coeff))
    lib.chkfile.save(str(chk), f"atomterm/{key.address}/ncas", int(ncas))


def _facts(key: CellKey, mol, mf, cas, ncore: int, ncas: int,
           na: int, nb: int, c_pur=None):
    """The cell layer's one facts call: the fold's spherical convention
    and the CASCI anisotropy report are the only cell-specific
    choices; everything else is :func:`atomterm.stages.stage_facts`."""
    from atomterm.stages import stage_facts

    c_act = cas.mo_coeff[:, ncore:ncore + ncas]
    return stage_facts(
        mol, mf, cas, c_act, ncore, na, nb, tag=key.address,
        spherical_core=key.fold != "raw",
        # every CASCI cell reports its reference core's anisotropy (the
        # raw fold's is the detector; the pur fold's is the ~0
        # certificate; row parity across the fold block)
        anisotropy_ref=(mf, c_pur) if key.method == "casci" else None,
    )


# the declared key-sets ----------------------------------------------

def spine_keys(atom: str, window: str, sector: tuple[int, int]) -> list[CellKey]:
    """The default cell subset per atom (README · Methods, "the spine"):
    the three folds; frozen x purified under both objectives (the
    certified class); one free row per substrate (the basin detectors).
    Diagnostics beyond these run as needed. A SUBSET view of the grid
    (every spine key is in :func:`window_grid`'s product); pinned here
    because README quotes it verbatim.
    """
    sec = tuple(sector)
    return [
        CellKey(atom, window, sec, fold="raw"),
        CellKey(atom, window, sec, fold="avg"),
        CellKey(atom, window, sec, fold="frac"),
        CellKey(atom, window, sec, method="casscf", mask="frozen",
                substrate="pur"),
        CellKey(atom, window, sec, method="sa-casscf", mask="frozen",
                substrate="pur"),
        CellKey(atom, window, sec, method="casscf"),
        CellKey(atom, window, sec, method="sa-casscf",
                substrate="pur"),
    ]


def _no_label(window: str) -> str:
    """The full-span occupation-group prefix label of a window: k = its
    shell count ({3d,4s} -> the [s,d] groups -> ``3d4sno2``; each shell's
    m-set is one occupation group of the term-averaged density)."""
    return f"{window}no{len(re.findall(r'\d+[spdfgh]', window))}"


def window_grid(spec, window: str | None = None) -> list[CellKey]:
    """One window's full factor grid, in canonical render order.

    Parameters
    ----------
    spec : the :class:`~atomterm.cases.CaseSpec` (duck-typed: symbol,
        shell_patterns, nelec_act, spin are read).
    window : the window label; None = the case's primary window.

    The ORDER is presentation and lives HERE, once -- renderers receive
    this list ready-made and format (branching on key identity only):
    the CASCI fold block first (identity: method=casci, start=ao; the
    pur fold third), then the optimizer product method x mask x
    substrate x start with the virtual-mix diagnostic after the CASSCF
    free/polarized starts, and the NO start addressing the full-span
    group prefix (:func:`_no_label`). mo x pur IS part of the grid
    (the mo branch re-projects the canonical columns on the pur
    span); radial stays ao-only (validated at the key).
    """
    from atomterm.system import sector_na_nb

    w = window or "".join(spec.shell_patterns)
    atom, sec = spec.symbol, sector_na_nb(spec.nelec_act, spec.spin)
    keys = [
        CellKey(atom, w, sec, fold="raw"),
        CellKey(atom, w, sec, fold="avg"),
        CellKey(atom, w, sec, substrate="pur"),   # the CASCI pur fold
        CellKey(atom, w, sec, fold="frac"),
    ]
    for m in ("casscf", "sa-casscf"):
        for mask in (None, "radial", "frozen"):
            for sub in ("pol", "pur"):
                starts = ("ao", "mo", "no") if mask != "radial" else ("ao",)
                for st in starts:
                    keys.append(CellKey(atom, _no_label(w) if st == "no"
                                        else w, sec, method=m, start=st,
                                        mask=mask, substrate=sub))
                if m == "casscf" and mask is None and sub == "pol":
                    keys.append(CellKey(atom, w, sec, method="casscf",
                                        mask="virtmix", substrate=sub))
    return keys


def grid_keys(spec) -> list[CellKey]:
    """The case's runnable set: the primary window's full grid
    (:func:`window_grid`) -- the spine folds, the optimizer product, the
    CASCI pur fold, and the virtual-mix diagnostic; every key distinct
    by construction (the spine is a subset view). What regen executes;
    planned windows are :func:`planned`."""
    return window_grid(spec)


def diagnostics(spec) -> list[CellKey]:
    """The virtual-mix diagnostic rows (no core<->window exchange) --
    the subset of the grid tagged mask='virtmix'."""
    return [k for k in grid_keys(spec) if k.mask == "virtmix"]


def planned(spec) -> list[CellKey]:
    """The planned windows' full grids (declared, not yet run): each of
    ``spec.planned_windows`` gets the same factor grid as the primary
    window (:func:`window_grid`) -- CASCI anchors first when the runs
    come, then the optimizer cells."""
    return [k for w in spec.planned_windows for k in window_grid(spec, w)]


def lattice_keys(spec) -> list[CellKey]:
    """The declared free-ion d2 lattice (``spec.lattice_windows``): the
    CAS(2,·) series -- the ns2 pair closed in the core, sector (2,0) =
    Ms=1 -- over growing m-complete shell spans.

    Three constructions per span: the raw shell AOs (exactly closed),
    the canonical MO slice (the m-delocalized control), and the
    full-span group-prefix NOs of the source's own term-averaged
    density (:func:`_no_label`; NO == AOs by construction -- the
    machinery certificate, not a methodology).

    Parameters
    ----------
    spec : the declaring :class:`~atomterm.cases.CaseSpec` (empty
        ``lattice_windows`` -> empty list).

    Examples
    --------
    Ti (ncore=10): 3d -> CAS(2,5), 3d4p -> CAS(2,8), 3d4p4d ->
    CAS(2,13); the NO labels carry the group count (3d4pno2).
    """
    sec = (2, 0)
    keys = []
    for w in spec.lattice_windows:
        keys.append(CellKey(spec.symbol, w, sec, method="casci",
                            fold="avg"))
        keys.append(CellKey(spec.symbol, w, sec, method="casci",
                            start="mo", fold="avg"))
    for w in spec.lattice_windows:
        keys.append(CellKey(spec.symbol, _no_label(w), sec,
                            method="casci", start="no", fold="avg"))
    return keys
