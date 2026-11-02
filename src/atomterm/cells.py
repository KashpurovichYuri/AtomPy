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
different addresses.  A LANE window carries its own fold as a tag on
the label (``3d3p@semi``, ``3s3p3d@1s2s2p4s`` --
:func:`atomterm.windows.fold_shells` the tag grammar): the tag
declares the cell's core partition AND marks the pinned-lane cell
(:func:`lane_keys`, :func:`run_lane_cell`) -- a tagged window never
collides with a free-route label.

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
    "lane_keys",
    "run_lane_cell",
    "load_lane_cell",
    "void_reason",
]

METHODS = ("casci", "casscf", "sa-casscf")
STARTS = ("ao", "mo", "no")       # window constructions (the span is the window)
MASKS = ("radial", "frozen", "virtmix")   # None = free optimization
SUBSTRATES = ("pol", "pur")
FOLDS = ("raw", "avg", "frac")

# the nominal window label: shell tokens like "2p", "3d4s" -- lowercase
# alphanumerics, no separators (the address splits on "_"); a lane
# window appends its fold tag: "3d3p@semi", "3s3p3d@1s2s2p4s"
_WINDOW_RE = re.compile(r"[0-9][0-9a-z]*(?:@[0-9a-z]+)?\Z")


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
        span, start = the columns defining it).  A LANE window
        carries its own fold as a tag (``"3d3p@semi"`` --
        :meth:`fold_shells` resolves it): the tag declares the cell's
        core partition AND marks the lane (a tagged window never
        collides with a free-route label; NO windows take none).
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
        Lane cells carry ``(pur, avg)`` -- the certified-fold class:
        an exactly-closed shell-span fold whose field is exactly
        m-scalar (the uniform construction of :func:`run_lane_cell`).
        For fold-tagged windows the substrate axis NAMES THE ROUTE:
        ``pur`` = the pinned lane class (run_lane_cell), ``pol`` = the
        free driver (:func:`run_cell`'s fold branch -- the twins).
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
                "(lowercase alphanumerics, e.g. '3d4s', optionally with "
                "a fold tag '3d3p@semi'; no separators)")
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
        if self.no_source is not None and self.fold_tag is not None:
            raise ValueError(
                f"window {self.window!r}: NO windows carry no fold tag "
                "(their fold is the source window's own)")
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
        """The window's span as shell selectors: ``"3d4s"`` ->
        ``["3d", "4s"]``; a lane window's fold tag is not part of the
        span."""
        return re.findall(r"\d+[spdfgh]", self.span_label)

    @property
    def span_label(self) -> str:
        """The window label without its fold tag (the bare span)."""
        return self.window.partition("@")[0]

    @property
    def fold_tag(self) -> str | None:
        """The window's ``@fold`` tag (None = the case's core prefix;
        :meth:`fold_shells` resolves a tag to shells)."""
        span, sep, tag = self.window.partition("@")
        return tag if sep else None

    def fold_shells(self) -> tuple[str, ...]:
        """The key's fold tag resolved to shells
        (:func:`atomterm.windows.fold_shells`; the class rule's open d
        = the LOWEST-n d shell of the window span -- every lane window
        contains the valence d).

        Raises
        ------
        ValueError
            the window carries no tag (its fold is the case's prefix,
            not the key's to resolve), or the tag is not a fold name.
        """
        from atomterm.windows import fold_shells

        tag = self.fold_tag
        if tag is None:
            raise ValueError(
                f"window {self.window!r} carries no fold tag; its fold "
                "is the case's core prefix")
        ds = [p for p in self.shell_patterns if p.endswith("d")]
        open_d = min(ds, key=lambda p: int(p[:-1])) if ds else None
        return fold_shells(tag, open_d=open_d)

    @property
    def no_source(self) -> tuple[str, int] | None:
        """For NO windows: ``(source label, group prefix k)``. The label
        grammar is ``<source>no<k>`` (e.g. ``"3d4sno2"`` = the second
        occupation-group prefix of the 3d4s source's term-averaged
        density) -- shell letters cannot contain "no", so the tail is
        unambiguous. The prefix IS part of the window identity: two NO
        cells of one source differ in k and get different addresses.
        """
        m = re.fullmatch(r"(.+?)no([1-9]\d*)", self.span_label)
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

# the dense-meter budget: sectors past this determinant count skip the
# exact [H, L] commutators (O(ndet^2) memory) -- the ladder's em-dash
_COMM_CAP = 6000


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

    # the route-class grammar first, before any store access: a
    # fold-tagged window on the free route is ao/pol/avg (the pur
    # class is run_lane_cell's; raw/frac of a tagged window are
    # meaningless -- the tag is the partition, the fold exactly
    # m-scalar)
    if key.fold_tag is not None and (key.substrate, key.fold,
                                     key.start) != ("pol", "avg", "ao"):
        raise ValueError(
            f"{key.address}: a free-route ladder window is ao/pol/avg "
            "(the pur class is run_lane_cell's; the tag is the "
            "partition -- raw/frac do not apply)")

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

    # a fold-tagged window: the tag names the core partition (the
    # route-class grammar validated above); the fold's raw-shell
    # Loewdin span IS the substrate.  The construction validations
    # mirror the pinned lane's.
    c_run = None
    if key.fold_tag is not None:
        from atomterm.windows import (closed_electrons, lowdin,
                                      select_shell_ao)
        fold = key.fold_shells()
        overlap = sorted(set(key.shell_patterns) & set(fold))
        if overlap:
            raise ValueError(
                f"{key.address}: {overlap} sit in both the window and "
                "its fold")
        nelec = mol.nelectron - closed_electrons(fold)
        if na + nb != nelec:
            raise ValueError(
                f"{key.address}: the sector holds {na + nb} active "
                f"electrons; the fold {list(fold)} leaves {nelec}")
        fold_rows = [r for sh in fold for r in select_shell_ao(mol, sh)]
        if not fold_rows:
            raise ValueError(
                f"{key.address}: the basis lacks a fold shell")
        c_run = lowdin(np.eye(mol.nao_nr())[:, fold_rows], mf.get_ovlp())
        ncore = c_run.shape[1]

    c_pur = (None if c_run is not None
             else purified_core_orbitals(mol, mf, ncore)[0]
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
        core_idx = (list(range(ncore)) if c_run is not None
                    else occupied_core_idx(mf, ncore))
        if key.mask == "radial":
            mo_full, labels = radial_basis(
                mf, core_idx, act_idx, c_core=c_run if c_run is not None
                else c_pur)
            mask = mask_radial_rotations(labels, ncore, ncas)
        else:
            mo_full = build_full_mo(mf, core_idx, act_idx,
                                    c_core=c_run if c_run is not None
                                    else c_pur)
    else:  # start == "mo": PySCF's default start (canonical MOs,
        # canonical window). On pol the canonical MO matrix IS the
        # substrate; on pur the canonical columns are re-orthogonalized
        # against the pur core span (they are S-orthonormal against
        # the core they were SOLVED in, not the pur one -- exactly the
        # input-window defect build_full_mo's projection removes at the
        # source). Honest
        # expectation: the window's own m-delocalization survives
        # (measured: raw-MO spans never rotation-close) -- these
        # cells measure that, they do not fix it.
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
                                   c_core=c_run if c_run is not None
                                   else c_pur)
            _, _, weights, _ = ground_block(ops0)  # the SA-anchoring rule
            weights = list(weights)
        cas = constrained_casscf(mf, ncas, (na, nb), str(chk),
                                 mo_full=mo_full, mask=mask, weights=weights,
                                 name=key.address, **common)

    _store_cell(chk, key, cas, ncas,
                  store_name=f"{key.address}_roots"
                            if key.method == "sa-casscf"
                  else None)
    # the dense stage analysis (operator matrices over the sector) is
    # budgeted like the lane's meters: past _COMM_CAP the cell stores
    # KERNEL-ONLY (the run's own energies; the meters read em-dash) --
    # None is the caller's signal (:func:`kernel_cell_payload`)
    if sector_dimension(ncas, na, nb) > _COMM_CAP:
        return None
    return _facts(key, mol, mf, cas, ncore, ncas, na, nb,
                  c_pur=c_run if c_run is not None else c_pur)


# analysis + persistence ---------------------------------------------

def analyze_cell(key: CellKey, mol, mf, workdir: str | Path = "elements"):
    """Re-analyze a stored cell (load + one facts pass) -- no mol/mf-bound
    optimizer re-run, no re-solving. The two-function lean: loading is
    pure data (:func:`load_cell`), analysis needs the live reference.

    Fold-tagged windows analyze against their tag's raw-shell span
    (the core size is the fold's column count, not the prefix); sectors
    past the dense-meter budget return None -- the kernel-only class
    (:func:`kernel_cell_payload`).  ``converged`` defaults True (a
    trusted cached stage).
    """
    from atomterm.casscf import CASResult
    from atomterm.core import purified_core_orbitals
    from atomterm.system import ncore_of
    from atomterm.windows import (closed_electrons, lowdin,
                                  select_shell_ao, select_valence_ao)

    data = load_cell(key, workdir)
    cas = CASResult(data.energies, data.ci, data.mo_coeff,
                    None, data.core_energy)
    na, nb = key.sector
    ncore = ncore_of(mol.nelectron, na + nb)
    ncas = data.ncas or len(select_valence_ao(mol, key.shell_patterns))
    c_run = None
    if key.fold_tag is not None:
        fold = key.fold_shells()
        fold_rows = [r for sh in fold for r in select_shell_ao(mol, sh)]
        c_run = lowdin(np.eye(mol.nao_nr())[:, fold_rows], mf.get_ovlp())
        ncore = c_run.shape[1]
    if sector_dimension(ncas, na, nb) > _COMM_CAP:
        return None
    c_pur = (None if c_run is not None
             else purified_core_orbitals(mol, mf, ncore)[0]
             if key.substrate == "pur" else None)
    return _facts(key, mol, mf, cas, ncore, ncas, na, nb,
                  c_pur=c_run if c_run is not None else c_pur)


def cell_payload(key: CellKey, facts) -> dict:
    """The persisted payload of a cell: the stage facts + the address."""
    out = facts.payload(key.method)
    out["cell"] = key.address
    return out


def kernel_cell_payload(key: CellKey,
                        workdir: str | Path = "elements") -> dict:
    """The kernel-only payload of a big-sector cell (run_cell returns
    None past the dense-meter budget): the run's own energies from the
    store, every dense-meter field NaN -- the tables' em-dash class.
    The kernel totals are the values that compare; the meters are
    structural zeros for the certified folds."""
    data = load_cell(key, workdir)
    return {
        "cell": key.address,
        "kind": key.method,
        "n_det": int(sector_dimension(data.ncas, *key.sector)),
        "energies_pyscf": [float(x) for x in np.atleast_1d(data.energies)],
        "ground_l2": float("nan"),
        "comm_norm_h_lz": float("nan"),
        "comm_norm_h_l2": float("nan"),
        "converged": True,
    }


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


# the pinned lane ------------------------------------------------------


def _root_vectors(ci):
    """A solver's return as a list of root vectors (nroots = 1 returns
    a bare array, nroots > 1 a stacked one)."""
    ci = np.asarray(ci)
    return [ci] if ci.ndim == 2 else list(ci)


def run_lane_cell(key: CellKey, mol, mf, term, conv_tol: float = 1e-12,
                  anchor: bool = False, fresh: bool = False, verbose: int = 0,
                  workdir: str | Path = "elements") -> dict:
    """Run one pinned-lane cell: the CAS ladder's UNIFORM construction.
    The fold's raw AO-shell span, Löwdin-orthonormalized; the window
    placed by :func:`atomterm.windows.build_full_mo` against it; one
    ``wfnsym``-pinned solver per branch component
    (:func:`atomterm.irreps.branch_d2h`); meters per irrep block
    (L2/S2/H are totally symmetric, so the block restriction is
    exact); the fold-independent bookkeeping ruler
    (:func:`atomterm.core.closed_shell_core_energy`).  Deterministic:
    raw AO spans, draw-free on both sides of the fold.

    Parameters
    ----------
    key : the lane cell (a :func:`lane_keys` key: ``casci/ao/-/pur/avg``
        with a tagged window).  Its sector must hold exactly the
        electrons its fold leaves active.
    mol, mf : the PINNED substrate -- ``mol`` a converged-reference
        D2h molecule, ``mf`` its ROHF.
    term : the target ``(L, parity)`` (:attr:`~atomterm.cases.CaseSpec.term`);
        parity is the configuration's total inversion parity.
    conv_tol : the pinned solvers' CI tolerance.
    anchor : also eigen-solve each block with OUR H and record the
        max deviation from the pinned solver (the conventions
        certificate; O(block^3) -- keep off large blocks).
    fresh : ignore a stored payload and re-run.
    verbose : print each block as it solves.
    workdir : the store root.

    Returns
    -------
    the payload dict (also stored under the key's address, reloadable
    through :func:`load_lane_cell`): per-component ``energies`` /
    ``l2`` / ``s2`` in branch order, the TERM ``total`` (component
    mean + the closed-shell core energy), the component spread
    ``ptp``, the sector size ``n_det``, ``anchor_dev``, and the dense
    assembled-treatment commutators ``comm_norm_h_lz`` /
    ``comm_norm_h_l2`` (NaN above :data:`_COMM_CAP` -- sectors past
    the dense-meter budget; the certified fold class makes them
    structural there).

    Raises
    ------
    ValueError
        the key is not a lane cell; the window and its fold overlap;
        the sector disagrees with the fold-derived electron count; the
        basis lacks a fold shell.
    """
    from pyscf import ao2mo, fci, lib

    from atomterm.core import closed_shell_core_energy
    from atomterm.irreps import (IRREP_IDS, branch_d2h, det_irreps,
                                 orb_irreps)
    from atomterm.operators import (active_space_integrals,
                                    build_operator_parts,
                                    ci_in_detbasis_order)
    from atomterm.detbasis import sector_dimension, slater_condon_matrix
    from atomterm.windows import (build_full_mo, closed_electrons, lowdin,
                                  select_shell_ao, select_valence_ao)

    if (key.method, key.start, key.mask, key.substrate,
            key.fold) not in (("casci", "ao", None, "pur", "avg"),
                              ("sa-casscf", "ao", None, "pur", "avg")):
        raise ValueError(
            f"a lane cell is casci/sa-casscf at ao/-/pur/avg (the "
            f"certified fold class); {key.address} is not one")
    chk = cell_chk(key, workdir)
    if chk.exists() and not fresh:
        return load_lane_cell(key, workdir)
    chk.parent.mkdir(parents=True, exist_ok=True)

    L, parity = term
    fold = key.fold_shells()
    span = key.shell_patterns
    overlap = sorted(set(span) & set(fold))
    if overlap:
        raise ValueError(
            f"{key.address}: {overlap} sit in both the window and its fold")
    na, nb = key.sector
    nelec = mol.nelectron - closed_electrons(fold)
    if na + nb != nelec:
        raise ValueError(
            f"{key.address}: the sector holds {na + nb} active electrons; "
            f"the fold {list(fold)} leaves {nelec}")
    fold_rows: list[int] = []
    for sh in fold:
        rows = select_shell_ao(mol, sh)
        if not rows:
            raise ValueError(
                f"{key.address}: the basis has no {sh} shell to fold")
        fold_rows.extend(rows)
    s = mf.get_ovlp()
    n = mol.nao_nr()
    c_fold = lowdin(np.eye(n)[:, fold_rows], s)
    nf = c_fold.shape[1]
    rows = select_valence_ao(mol, span)
    mo = build_full_mo(mf, list(range(nf)), rows, c_core=c_fold)
    c = mo[:, nf:nf + len(rows)]
    ncas = c.shape[1]
    h1, _ = active_space_integrals(mol, mf, c, nf, c_core=c_fold)
    eri = ao2mo.kernel(mol, c)
    orbsyms = orb_irreps(mol, c)
    parts = build_operator_parts(mol, mf, c, nf, na, nb, c_core=c_fold)
    dets = parts["dets"]
    ids = det_irreps(orbsyms, na, nb)

    es, l2s, s2s = [], [], []
    anchor_dev = 0.0
    ensemble = None
    if key.method == "sa-casscf":
        # the pinned SA: the component ensemble optimized
        # (state_average_mix_ over the pinned solvers), then the exact
        # per-component roots re-solved in the frozen converged window
        # -- the two-pass pattern; the meters read the READOUT roots
        from atomterm.casscf import run_irrep_sa

        out = run_irrep_sa(mol, mf, c, nf, na, nb, L, parity,
                           conv_tol=conv_tol, c_core=c_fold)
        ensemble = float(out["e_tot"])
        # the term TOTAL is the kernel ensemble (the optimizer may
        # rotate same-irrep core/window columns -- the readout's
        # absolute energies are of the run's own partition; its
        # COMPONENT SPREAD and meters are the diagnostics)
        for irrep, mult in branch_d2h(L, parity).items():
            e_arr, ci = out["roots"][irrep]
            e_arr = np.atleast_1d(np.asarray(e_arr, dtype=float))
            es.extend(float(x) for x in e_arr[:mult])
            sel = ids == IRREP_IDS[irrep]
            block = [d for d, m in zip(dets, sel) if m]
            for meters, op in ((l2s, "L2"), (s2s, "S2")):
                mb = slater_condon_matrix(block, one_body=parts[op][0],
                                          two_body=parts[op][1])
                for vec in _root_vectors(ci):
                    w = ci_in_detbasis_order(vec, ncas, na, nb)[sel]
                    meters.append(float(np.real(w @ mb @ w)))
                del mb
    else:
        for irrep, mult in branch_d2h(L, parity).items():
            sol = fci.direct_spin1_symm.FCI(mol)
            sol.orbsym, sol.wfnsym = orbsyms, IRREP_IDS[irrep]
            sol.nroots, sol.conv_tol = mult, conv_tol
            if verbose:
                print(f"  lane {key.address} block {irrep} nroots={mult}")
            e, ci = sol.kernel(h1, eri, ncas, (na, nb))
            e = np.atleast_1d(e)[:mult]
            es.extend(float(x) for x in e)
            sel = ids == IRREP_IDS[irrep]
            block = [d for d, m in zip(dets, sel) if m]
            w_blk = []
            for vec in _root_vectors(ci):
                w = ci_in_detbasis_order(vec, ncas, na, nb)[sel]
                # the pinned vector lives entirely in its block
                assert abs(np.linalg.norm(w) - 1.0) < 1e-8, irrep
                w_blk.append(w)
            if anchor:
                hb = slater_condon_matrix(block,
                                          one_body=parts["H"][0],
                                          two_body=parts["H"][1])
                ours = np.linalg.eigvalsh(hb)[:mult]
                anchor_dev = max(anchor_dev,
                                 float(np.max(np.abs(np.sort(e) - ours))))
                del hb
            for meters, op in ((l2s, "L2"), (s2s, "S2")):
                mb = slater_condon_matrix(block,
                                          one_body=parts[op][0],
                                          two_body=parts[op][1])
                for w in w_blk:
                    meters.append(float(np.real(w @ mb @ w)))
                del mb

    es = np.asarray(es)
    # the assembled-treatment meters: exact dense commutators where
    # the sector fits memory; above the cap the certified fold class
    # makes [H, L] ~ 0 structural (the component spread ptp is the
    # energy-space certificate) and the columns read em-dash
    comm_lz = comm_l2 = float("nan")
    if len(dets) <= _COMM_CAP:
        from atomterm.detbasis import commutator, frobenius_norm

        h_m = slater_condon_matrix(dets, one_body=parts["H"][0],
                                   two_body=parts["H"][1])
        lz_m = slater_condon_matrix(dets, one_body=parts["Lz"][0],
                                    two_body=parts["Lz"][1])
        l2_m = slater_condon_matrix(dets, one_body=parts["L2"][0],
                                    two_body=parts["L2"][1])
        comm_lz = float(frobenius_norm(commutator(h_m, lz_m)))
        comm_l2 = float(frobenius_norm(commutator(h_m, l2_m)))
        del h_m, lz_m, l2_m
    if key.method == "sa-casscf":
        # the kernel ensemble is the total (the run's own partition
        # may have moved under same-irrep rotations); the readout's
        # degenerate components carry ptp and the meters
        total = ensemble
    else:
        total = float(es.mean() + closed_shell_core_energy(mf, c_fold))
    payload = {
        "cell": key.address,
        "term": [int(L), str(parity)],
        "energies": [float(x) for x in es],
        "total": float(total),
        "ptp": float(np.ptp(es)),
        "l2": l2s,
        "s2": s2s,
        "n_det": int(sector_dimension(ncas, na, nb)),
        "anchor_dev": float(anchor_dev),
        "comm_norm_h_lz": comm_lz,
        "comm_norm_h_l2": comm_l2,
        "fold": list(fold),
    }
    if ensemble is not None:
        payload["ensemble"] = ensemble
    for field_name, value in payload.items():
        lib.chkfile.save(str(chk), f"atomterm/{key.address}/{field_name}",
                         value)
    return payload


def load_lane_cell(key: CellKey, workdir: str | Path = "elements") -> dict:
    """Load a stored lane cell's payload (:func:`run_lane_cell`'s
    reader half) -- pure data, no mol/mf needed, normalized to the
    stored types (HDF5 hands strings back as bytes).

    Raises
    ------
    FileNotFoundError / KeyError
        no chkfile for the address / an incomplete store (regenerate
        the cell).
    """
    from pyscf import lib

    def _stored(v):
        if isinstance(v, bytes):
            return v.decode()
        if isinstance(v, np.generic):
            return v.item()
        if isinstance(v, (list, tuple, np.ndarray)):
            return [_stored(x) for x in v]
        return v

    chk = cell_chk(key, workdir)
    if not chk.exists():
        raise FileNotFoundError(
            f"no stored lane cell {key.address} at {chk}; run it first "
            "(atomterm.cells.run_lane_cell)")
    out = {"cell": key.address}
    for field_name in ("term", "energies", "total", "ptp", "l2", "s2",
                       "n_det", "anchor_dev", "comm_norm_h_lz",
                       "comm_norm_h_l2", "fold"):
        value = lib.chkfile.load(str(chk),
                                 f"atomterm/{key.address}/{field_name}")
        if value is None:
            raise KeyError(
                f"lane store {chk} is incomplete (no {field_name}); "
                "regenerate the cell")
        out[field_name] = _stored(value)
    # optional fields (method-dependent): the SA ensembles
    value = lib.chkfile.load(str(chk),
                             f"atomterm/{key.address}/ensemble")
    if value is not None:
        out["ensemble"] = _stored(value)
    return out


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
    """The declared lattice (``spec.lattice_windows`` at
    ``spec.lattice_sector``): CASCI over growing m-complete shell
    spans -- WHAT the series studies is the declared sector (the
    free-ion series with the ns2 pair folded into the core, or the
    full valence).

    Two constructions per span: the raw shell AOs (exactly closed)
    and the canonical MO slice (the m-delocalized control).  NO
    variants are NOT a lattice arm: the own-span group-prefix label
    presupposes a term-averaged source density, which is a
    free-ion-series property (the ground block is the whole multiplet)
    -- on the valence series the ground block is a split +-pair whose
    m-unequal density groups by noise (measured, 2026-11-02); the NO
    window variant lives in the methods grid, where its donor is
    declared per key.

    Parameters
    ----------
    spec : the declaring :class:`~atomterm.cases.CaseSpec` (empty
        ``lattice_windows`` -> empty list).

    Examples
    --------
    Ti (free-ion series, sector (2,0)): 3d -> CAS(2,5), 3d4p ->
    CAS(2,8), 3d4p4d -> CAS(2,13).  Sc (valence series, sector
    (2,1)): 3d4s -> CAS(3,6).
    """
    sec = spec.lattice_sector
    keys = []
    for w in spec.lattice_windows:
        keys.append(CellKey(spec.symbol, w, sec, method="casci",
                            fold="avg"))
        keys.append(CellKey(spec.symbol, w, sec, method="casci",
                            start="mo", fold="avg"))
    return keys


def lane_keys(spec) -> list[CellKey]:
    """The declared pinned-lane legs (``spec.lane_windows``): one
    CI-level lane cell per leg, run by :func:`run_lane_cell` (the legs
    need ``spec.term`` declared -- validated at the CaseSpec).

    The SECTOR derives from the leg's own fold: activating a shell
    moves its electrons into the window (the 3d-metal semi fold leaves
    ``nelec_act + 4`` active), so a lane key never reuses the case's
    default sector blindly.

    Examples
    --------
    Ti (nelec_act=4, Ar prefix): ``3d4s@ar`` keeps sector (3, 1);
    ``3d3p@semi`` grows to (5, 3); the M-shell folds to (6, 4).
    """
    from dataclasses import replace

    from atomterm.system import sector_na_nb
    from atomterm.windows import closed_electrons

    base = sector_na_nb(spec.nelec_act, spec.spin)
    keys = []
    for label in spec.lane_windows:
        k = CellKey(spec.symbol, label, base, method="casci",
                    substrate="pur", fold="avg")
        nelec = (spec.nelec_act + 2 * spec.ncore
                 - closed_electrons(k.fold_shells()))
        keys.append(replace(k, sector=sector_na_nb(nelec, spec.spin)))
    return keys


def route_keys(spec) -> list[CellKey]:
    """The ROUTE AXIS over the ladder windows: every pinned CI leg's
    twins -- the free CI (a survey row; same window, same fold, the
    free driver), the free CASSCF and free SA-CASSCF, and the pinned
    SA-CASSCF (the two-pass ensemble through the lane runner).  The
    substrate axis names the route for tagged windows (pur = pinned,
    pol = free).  The declared plan; unrun members render pending.

    Examples
    --------
    The ``3d3p4d@semi`` leg spawns ``casci/casscf/sa-casscf`` at pol
    (free) and ``sa-casscf`` at pur (pinned) -- four twins beside its
    pinned CI.
    """
    from dataclasses import replace

    out = list(lane_keys(spec))
    for leg in lane_keys(spec):
        for method in ("casci", "casscf", "sa-casscf"):
            out.append(replace(leg, method=method, substrate="pol"))
        out.append(replace(leg, method="sa-casscf"))
    return out


def void_reason(address, scenarios) -> str | None:
    """The RUN GATE: why the cell at ``address`` (a
    :class:`CellKey` or its address string) is void under the declared
    scenarios, or None if it runs.  A void cell stays addressable -- a
    full regen materializes it, which re-tests the scenario's dated
    claim.  First enabled match wins (scenario order: campaign
    defaults, then the row's).
    """
    a = address.address if isinstance(address, CellKey) else address
    for s in scenarios:
        if s.enabled and s.gates(a):
            return f"{s.name}: {s.reason}"
    return None
