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
* ``core`` -- the core FIELD convention of the analysis: ``raw`` = the
  SCF density as converged (the anisotropy detector), ``avg`` = its
  spherical average (the ruler), ``frac`` = a different reference SCF
  entirely (fractional occupations).

  Composition: (raw, pol) detector; (avg, pol) ruler; (avg, pur) the
  certified core; (frac, pol) own reference. On ``pur`` the density is
  exactly m-scalar, so ``raw`` degenerates with ``avg`` --
  the key CANONICALIZES ``core="raw"`` to ``"avg"`` there, and
  ``frac`` on ``pur`` is rejected. ``substrate`` is the span, ``core``
  is the field.

Address grammar
---------------
Manifest reading order: where / the physics box /
the objective / the optimizer freedom / the core::

    <atom>_<window>_<na>-<nb>_<method>_<start>_<mask|->_<substrate>_<core>

with SINGULAR axis values throughout (``ao``/``mo``/``moc``/``no``), the
sector as ``na-nb`` (fixes Ms AND nelec), and ``-`` for no mask
(free optimization; a mask exists only for optimizer methods --
``virtmix`` is a MASK -- :func:`atomterm.casscf.mask_no_core_window_mix
-- not a start). NO windows append a prefix index to the label:
``<source>no<k>`` (e.g. ``3d4sno2``) = the k-th occupation-group
prefix of the source's term-averaged density -- the prefix is part of
the window identity, so differently-sized NO cells of one source get
different addresses.  A LANE window carries its own core as a tag on
the label (``3d3p@semi``, ``3s3p3d@1s2s2p4s`` --
:func:`atomterm.windows.core_shells` the tag grammar): the tag
declares the cell's core partition AND marks the pinned-lane cell
(:func:`lane_keys`, :func:`run_lane_cell`) -- a tagged window never
collides with a free-route label.  A PRIOR-AXIS window carries its
distillation source as the tag instead (``3d4s@s2`` --
:attr:`CellKey.source_tag` the grammar): ``@r<k>`` the head
configuration's density block, ``@s<k>`` the configuration-series
average (:func:`atomterm.basis.series_window`) -- the one ``@``-slot
holds a core or a source tag, never both.

The REFERENCE FIELD sits outside the eight slots: a cell runs on the
neutral-atom SCF (``field="neutral"`` -- the default, and every store
predating the axis) or on the folded reference
(:func:`atomterm.reference.run_folded` -- the confined-core field
with the valence re-relaxed, the lane's reference class).  The field
serializes as an OPTIONAL ninth slot, appended ONLY when non-default
(``..._<substrate>_<core>_folded``): an eight-slot address parses as
``neutral`` and re-serializes unchanged -- the store identities minted
before the axis keep their exact names.  A folded key requires a
core-tagged window: the folded field is defined BY its confined span,
and the tag is that span's name.

One real address, decoded slot by slot::

    ti_3d4sno2_3-1_sa-casscf_no_-_pol_avg
    │  │        │   │         │  │ │   └ core=avg   (the averaged core field)
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
    "kernel_cell_payload",
    "load_cell",
    "no_label",
    "prior_tag_of",
    "run_cell",
    "analyze_cell",
    "spine_keys",
    "window_grid",
    "grid_keys",
    "planned",
    "lattice_keys",
    "lane_keys",
    "route_keys",
    "run_lane_cell",
    "load_lane_cell",
    "radial_refused",
    "void_reason",
]

METHODS = ("casci", "casscf", "sa-casscf")
STARTS = ("ao", "mo", "moc", "no")  # window constructions (the span is the window)
MASKS = ("radial", "frozen", "virtmix", "frozencore")   # None = free optimization
SUBSTRATES = ("pol", "pur")
CORES = ("raw", "avg", "frac")
FIELDS = ("neutral", "folded")  # the reference field (an optional 9th slot)

# the nominal window label: shell tokens like "2p", "3d4s" -- lowercase
# alphanumerics, no separators (the address splits on "_"); a lane
# window appends its core tag: "3d3p@semi", "3s3p3d@1s2s2p4s"; a
# prior-axis window appends its source tag: "3d4s@s2" (the @-slot
# holds ONE tag -- a core or a source, never both)
_WINDOW_RE = re.compile(r"[0-9][0-9a-z]*(?:@[0-9a-z]+)?\Z")

# the prior-axis source tag: @r<k> (the head configuration's density
# block) | @s<k> (the series average) -- letter-first, k >= 1; disjoint
# from the core vocabulary (noble names, "semi", digit-first shell
# lists) by construction
_SOURCE_TAG_RE = re.compile(r"([rs])([1-9]\d*)\Z")


def _source_label_of(window: str) -> str:
    """The window label with any NO tail removed (``"3d4s@s2no2"`` ->
    ``"3d4s@s2"``) -- the part a core/source tag lives in."""
    m = re.fullmatch(r"(.+?)no([1-9]\d*)", window)
    return m.group(1) if m else window


def prior_tag_of(window: str) -> tuple[str, int] | None:
    """The label's ``@``-source tag ``(source, k)``, None when
    untagged -- the prior-axis grammar's ONE home: the key property,
    the case validation, and the renderer's tag display parse
    identically through here.  A NO tail is not part of the tag
    (``"3d4s@s2no2"`` -> ``("s", 2)``)."""
    span, sep, tag = _source_label_of(window).partition("@")
    if not sep:
        return None
    m = _SOURCE_TAG_RE.fullmatch(tag)
    return (m.group(1), int(m.group(2))) if m is not None else None

# the campaign key: a plain alphanumeric token (an element symbol or a
# declared label) or ONE hyphenated variant suffix (the uncontracted
# twins: "Ti-unc") -- no leading/trailing/double hyphen; the hyphen
# rides inside the address's first "_"-slot unharmed
_ATOM_RE = re.compile(r"[A-Za-z][A-Za-z0-9]*(?:-[A-Za-z0-9]+)?\Z")


# the cell key -------------------------------------------------------

@dataclass(frozen=True)
class CellKey:
    """One grid point of the manifest: the unit of run, load, and report.

    Every axis has an owner; the key only names the choice.

    Attributes
    ----------
    atom : campaign key, e.g. ``"Ti"``, ``"U14e"``, ``"Ti-unc"`` (the
        ``elements/<Atom>/`` dir; a plain element symbol, a declared
        ECP-atom label, or a hyphenated basis-variant label -- at most
        one hyphen, no leading/trailing).
    window : NOMINAL span label, e.g. ``"3d4s"`` -- names the table
        row; the construction is the ``start`` axis (window = the
        span, start = the columns defining it).  A LANE window
        carries its own core as a tag (``"3d3p@semi"`` --
        :meth:`core_shells` resolves it): the tag declares the cell's
        core partition AND marks the lane (a tagged window never
        collides with a free-route label; NO windows take none).
    sector : ``(na, nb)`` -- alpha/beta active electrons; fixes Ms
        and nelec together.
    method : ``"casci"`` | ``"casscf"`` | ``"sa-casscf"``.
    start : the window KIND -- where the window's information comes
        from, always expressed as AO-basis columns (the CI consumes
        AO-basis columns; the kind names the SOURCE, never the
        expression basis):
        ``"ao"`` -- the shell rows of the declared realization,
        DISTILLED when the case declares ``shell_radials`` (the
        member-summed radial Gram of the reference's significant
        occupied columns, top-k replicated on the m-set): MO-informed,
        AO-expressed, m-complete by construction -- the production
        kind;
        ``"mo"`` -- the canonical SCF slice (CONTROL: m-delocalized,
        m-incomplete by construction -- the SCF's hybrid-open gauge;
        the certified fingerprint is a scrambled term reading);
        ``"moc"`` -- the character-selected SCF valence MOs
        (:func:`atomterm.windows.select_valence_mo`; CONTROL: same
        disease, a different pick);
        ``"no"`` -- the occupation groups of a term-averaged density
        (natural-orbital kind; m-complete by sphericity -- the
        density sibling of the distilled kind; the label grammar is
        ``<source>no<k>``, see :attr:`no_source`).
    mask : optimizer restriction -- ``None`` (free) | ``"radial"`` |
        ``"frozen"`` | ``"virtmix"`` | ``"frozencore"`` (the core
        immobile, window+virtuals free); only meaningful for optimizer
        methods (rejected for ``casci``).
    substrate / core : the core axes -- see the module docstring.
        Lane cells carry ``(pur, avg)`` -- the certified-core class:
        an exactly-closed shell-span core whose field is exactly
        m-scalar (the uniform construction of :func:`run_lane_cell`).
        For core-tagged windows the substrate axis NAMES THE ROUTE:
        ``pur`` = the pinned lane class (run_lane_cell), ``pol`` = the
        free driver (:func:`run_cell`'s core branch -- the twins).
    field : ``"neutral"`` | ``"folded"`` -- the REFERENCE field the
        cell's substrate lives on (the module docstring).  Default
        ``neutral``, serialized only when ``folded``: an address
        without the tag IS a neutral cell (stores predating the axis
        need no migration).  A folded key requires a core-tagged
        window.
    """

    atom: str
    window: str
    sector: tuple[int, int]
    method: str = "casci"
    start: str = "ao"
    mask: str | None = None
    substrate: str = "pol"
    core: str = "avg"
    field: str = "neutral"

    def __post_init__(self):
        """Validate and canonicalize: axis vocabularies, void
        combinations, and the pur core rule (raw -> avg)."""
        if not (isinstance(self.atom, str)
                and 1 <= len(self.atom) <= 8
                and _ATOM_RE.match(self.atom)):
            raise ValueError(f"atom {self.atom!r} must be a campaign key "
                             "(an element symbol or a declared label "
                             "like 'U14e', 'Ti-unc')")
        object.__setattr__(self, "atom", self.atom.capitalize())
        if not (isinstance(self.window, str) and _WINDOW_RE.match(self.window)):
            raise ValueError(
                f"window {self.window!r} must be a nominal shell label "
                "(lowercase alphanumerics, e.g. '3d4s', optionally with "
                "a core tag '3d3p@semi'; no separators)")
        na, nb = self.sector
        if not (isinstance(na, int) and isinstance(nb, int)
                and na >= nb >= 0 and na + nb > 0):
            raise ValueError(f"sector {self.sector!r} must be (na >= nb >= 0, n > 0)")
        object.__setattr__(self, "sector", (na, nb))
        for axis, value, vocab in (
                ("method", self.method, METHODS),
                ("start", self.start, STARTS),
                ("substrate", self.substrate, SUBSTRATES),
                ("core", self.core, CORES),
                ("field", self.field, FIELDS)):
            if value not in vocab:
                raise ValueError(
                    f"{axis} {value!r} not in {vocab}")
        if self.field == "folded" and self.core_tag is None:
            raise ValueError(
                f"field='folded' needs a core-tagged window (the tag "
                f"names the confined span that defines the field); "
                f"{self.window!r} carries none")
        if (self.start == "no") != (self.no_source is not None):
            raise ValueError(
                f"start='no' and the window label must agree: NO windows "
                "are '<source>no<k>' (e.g. '3d4sno2'), other starts take "
                "a plain shell label")
        if self.no_source is not None and self.core_tag is not None:
            raise ValueError(
                f"window {self.window!r}: NO windows carry no core tag "
                "(their core is the source window's own)")
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
        # the pur-core degeneracy rule: on the purified substrate the
        # core density is exactly m-scalar, so raw == avg --
        # canonicalize instead of storing two addresses for one cell;
        # frac is a different REFERENCE, meaningless on pur
        if self.substrate == "pur":
            if self.core == "frac":
                raise ValueError("core='frac' is a different reference SCF; "
                                 "it cannot live on the purified substrate")
            if self.core == "raw":
                object.__setattr__(self, "core", "avg")

    @property
    def address(self) -> str:
        """The derived serialization: chk name, payload key, table-row
        id.  The field rides as an optional ninth slot, appended only
        when ``folded`` -- an eight-slot address is a neutral cell, so
        every store name minted before the axis is untouched."""
        na, nb = self.sector
        parts = [self.atom.lower(), self.window, f"{na}-{nb}",
                 self.method, self.start, self.mask or "-",
                 self.substrate, self.core]
        if self.field != "neutral":
            parts.append(self.field)
        return "_".join(parts)

    @classmethod
    def from_address(cls, address: str) -> "CellKey":
        """Parse an address (inverse of :attr:`address`).  Eight slots
        = a neutral cell (the pre-axis spelling); nine = the ninth
        names the field."""
        parts = address.split("_")
        if len(parts) == 8:
            parts.append("neutral")
        if len(parts) != 9:
            raise ValueError(
                f"address {address!r} does not have 8 slots (9 with "
                "the field tag): atom_window_na-nb_method_start_mask_"
                "substrate_core[_field]")
        atom, window, sector, method, start, mask, substrate, core, \
            field_name = parts
        na, nb = sector.split("-")
        return cls(atom, window, (int(na), int(nb)), method, start,
                   None if mask == "-" else mask, substrate, core,
                   field=field_name)

    @property
    def shell_patterns(self) -> list[str]:
        """The window's span as shell selectors: ``"3d4s"`` ->
        ``["3d", "4s"]``; a lane window's core tag is not part of the
        span."""
        return re.findall(r"\d+[spdfgh]", self.span_label)

    @property
    def span_label(self) -> str:
        """The window label without its core tag (the bare span)."""
        return self.window.partition("@")[0]

    @property
    def _source_label(self) -> str:
        """The label of the window's SOURCE construction: the bare
        label with any NO tail removed (``"3d4s@s2no2"`` ->
        ``"3d4s@s2"``) -- the part a core/source tag lives in."""
        return _source_label_of(self.window)

    @property
    def core_tag(self) -> str | None:
        """The window's ``@core`` tag (None = the case's core prefix;
        :meth:`core_shells` resolves a tag to shells).  A source tag
        (:attr:`source_tag`) is not a core."""
        span, sep, tag = self._source_label.partition("@")
        return tag if sep and self.source_tag is None else None

    @property
    def source_tag(self) -> tuple[str, int] | None:
        """The window's ``@source`` tag: ``(source, k)`` -- the PRIOR
        AXIS, the basis-preparation window family
        (:func:`atomterm.basis.series_window`): ``@r2`` = the head
        configuration's density block at k=2 (the reference source),
        ``@s1``/``@s2`` = the configuration-series average at k=1/k=2
        (the series source), k uniform over the window's shells.  The
        realization is the case's declared ``prior_bands`` (the
        widened, certificate-measured band lists k>1 needs).  None =
        the incumbent construction (the case's ``shell_radials`` on
        the reference source -- the stored cells' identity, frozen).

        The tag is part of the window identity: differently-sourced
        windows of one span get different addresses, and the basis
        axis of a comparison is thereby always declared, never
        implicit.
        """
        return prior_tag_of(self.window)

    def core_shells(self) -> tuple[str, ...]:
        """The key's core tag resolved to shells
        (:func:`atomterm.windows.core_shells`; the class rule's open d
        = the LOWEST-n d shell of the window span -- every lane window
        contains the valence d).

        Raises
        ------
        ValueError
            the window carries no tag (its core is the case's prefix,
            not the key's to resolve), or the tag is not a core name.
        """
        from atomterm.windows import core_shells

        tag = self.core_tag
        if tag is None:
            raise ValueError(
                f"window {self.window!r} carries no core tag; its core "
                "is the case's core prefix")
        ds = [p for p in self.shell_patterns if p.endswith("d")]
        open_d = min(ds, key=lambda p: int(p[:-1])) if ds else None
        return core_shells(tag, open_d=open_d)

    @property
    def no_source(self) -> tuple[str, int] | None:
        """For NO windows: ``(source label, group prefix k)``. The label
        grammar is ``<source>no<k>`` (e.g. ``"3d4sno2"`` = the second
        occupation-group prefix of the 3d4s source's term-averaged
        density) -- shell letters cannot contain "no", so the tail is
        unambiguous. The prefix IS part of the window identity: two NO
        cells of one source differ in k and get different addresses.
        The source label may itself carry a source tag
        (``"3d4s@s2no2"`` = the NO kind of the series window).
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
    core_energy : ``E_tot`` minus the active part (the core constant).
    ncas : window width -- needed exactly when it is not derivable
        from the label (NO windows); 0 = derive from the shell label.
    cycles : the optimizer's macro-cycle count (cost pricing +
        restart audits); 0 for runs without an optimizer and for
        stores predating the field.
    max_cycle : the DECLARED optimizer cap (``cycles == max_cycle``
        is the cap-hitting flag by identity); ``None`` for runs
        without a declared cap and for stores predating the field.
    conv_dE : the final macro-cycle |delta e_tot| (``None`` where
        ``cycles`` is 0).
    """

    key: CellKey
    mo_coeff: np.ndarray
    energies: np.ndarray
    ci: list
    core_energy: float
    ncas: int = 0
    cycles: int = 0
    max_cycle: int | None = None
    conv_dE: float | None = None


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
    cycles = lib.chkfile.load(str(chk), f"atomterm/{key.address}/cycles")
    max_cycle = lib.chkfile.load(str(chk),
                                 f"atomterm/{key.address}/max_cycle")
    conv_dE = lib.chkfile.load(str(chk), f"atomterm/{key.address}/conv_dE")
    return CellData(key, np.asarray(mo), energies, ci, core_energy,
                    int(ncas) if ncas is not None else 0,
                    int(cycles) if cycles is not None else 0,
                    int(max_cycle) if max_cycle is not None else None,
                    float(conv_dE) if conv_dE is not None else None)


# ------------------------------------------------------------- the execution

# the dense-meter budget: sectors past this determinant count skip the
# exact [H, L] commutators (O(ndet^2) memory) -- the methods table's em-dash
_COMM_CAP = 6000

# the dense-adjudication budget: CI lane blocks up to this determinant
# count get our exact H eigen-solved per block (the conventions
# certificate AND the correction source for the pinned solver's
# near-degenerate root skips; eigenvectors included, O(block^3))
_ANCHOR_CAP = 40000

# how deep the targeted selection scans the dense spectrum for the
# declared term's components (in units of the block's multiplicity --
# the target need not be the block's ground)
_SELECT_DEPTH = 3


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


def _lane_core(key: CellKey, mol, mf,
               shell_bands: dict | None = None
               ) -> tuple[tuple[str, ...], int, np.ndarray, dict]:
    """One home for the core resolution shared by the free core branch
    (:func:`run_cell`), the cache path (:func:`analyze_cell`), and the
    pinned lane (:func:`run_lane_cell`): resolve the tag's shells
    THROUGH THE SHELL LAYER (a chemical shell's AO realization is
    declared, not read off labels --
    :func:`atomterm.windows.shell_realization`), refuse the malformed
    keys, and return the core's shell tuple, its DOUBLY-OCCUPIED PAIR
    count (the chemistry: ``closed_electrons`` halved -- never the
    span's column count), the frozen core columns, and the core's
    certificate report.

    The frozen core is the reference's OWN core -- the occupied core
    columns (by occupation), S-projected onto the core span and
    orthonormalized within it.  The field and the bookkeeping then
    join the neutral anchor BY CONSTRUCTION on any basis.  The
    m-scalar anisotropy rule applies downstream (the averaged core,
    ``active_space_integrals``'s default).  The retired alternative
    -- the span's own closed-shell Roothaan
    (:func:`atomterm.core.confined_core_orbitals`) -- over-relaxes on
    a wide span (the bare core, solved without the valence); its
    measured record lives in the committed artifacts.

    Notes
    -----
    The retirement measurement: the confined-Roothaan core sat
    +124.3/+38.7 mEh off the anchors on the uncontracted twins,
    basis-stable per the contracted control
    (``elements/Ti/logs/folded_review.log``); the repair collapsed
    the offsets to the projected-density bookkeeping form.

    Raises
    ------
    ValueError
        a CHEMICAL shell sitting in both the window and its core (the
        partition is chemical: a pair is frozen or active, never both;
        the shells' BAND realizations may OVERLAP -- shared radial
        primitives are physical, and the window columns are projected
        against the frozen core span, which owns the shared rows);
        the sector disagreeing with the core-derived electron count;
        the basis lacking one of the core's bands.
    """
    from atomterm.windows import (closed_electrons, lowdin,
                                  occupied_core_idx, select_shell_ao,
                                  shell_realization)

    core = key.core_shells()
    overlap = sorted(set(key.shell_patterns) & set(core))
    if overlap:
        raise ValueError(
            f"{key.address}: {overlap} sit in both the window and its core")
    na, nb = key.sector
    nelec = mol.nelectron - closed_electrons(core)
    if na + nb != nelec:
        raise ValueError(
            f"{key.address}: the sector holds {na + nb} active electrons; "
            f"the core {list(core)} leaves {nelec}")
    core_rows: list[int] = []
    for sh in core:
        for band in shell_realization(sh, shell_bands):
            rows = select_shell_ao(mol, band)
            if not rows:
                raise ValueError(
                    f"{key.address}: the basis has no {band} band to freeze"
                    + (f" (the declared realization of shell {sh!r})"
                       if band != sh else ""))
            core_rows.extend(rows)
    # the shells' realizations may SHARE bands (physical, see above) --
    # the span is the UNION of rows: a duplicate row would make the
    # Gram exactly singular
    core_rows = sorted(set(core_rows))
    c_span = lowdin(np.eye(mol.nao_nr())[:, core_rows], mf.get_ovlp())
    n_pairs = closed_electrons(core) // 2
    # the reference's own core: the occupied columns, S-projected onto
    # the core span and orthonormalized within it -- the span hosts
    # the core shells' full family, so the projection is near-identity
    # (the in-span weights are the certificate)
    s = mf.get_ovlp()
    c_scf = mf.mo_coeff[:, occupied_core_idx(mf, n_pairs)]
    proj = c_span @ (c_span.conj().T @ (s @ c_scf))
    in_span = np.einsum("ai,ab,bi->i", c_scf.conj(), s,
                        proj).real.clip(0, 1)
    c_core = lowdin(proj, s)
    report = {"core": "scf", "in_span": in_span.round(6).tolist()}
    return core, n_pairs, c_core, report


def _key_core(key: CellKey, mol, mf,
              shell_bands: dict | None = None
              ) -> tuple[int, np.ndarray | None, np.ndarray | None]:
    """One home for the key->core resolution shared by the run
    (:func:`run_cell`), the cache path (:func:`analyze_cell`), and the
    kernel payload (:func:`kernel_cell_payload`): returns
    ``(ncore, c_core, c_pur)`` -- a tagged core's SCF core IS
    the core (``ncore`` follows the core's PAIR count, not its span
    width), else the pur substrate's purified span when declared,
    else ``None`` freezes the occupied prefix.  ``c_pur`` is separate
    because the NO-start branch and the anisotropy reference need the
    pur substrate even when a core supplied ``c_core``."""
    from atomterm.core import purified_core_orbitals
    from atomterm.system import ncore_of

    na, nb = key.sector
    ncore = ncore_of(mol.nelectron, na + nb)
    c_run = None
    if key.core_tag is not None:
        _, ncore, c_run, _w = _lane_core(key, mol, mf, shell_bands)
    c_pur = (None if c_run is not None
             else purified_core_orbitals(mol, mf, ncore)[0]
             if key.substrate == "pur" else None)
    return ncore, (c_run if c_run is not None else c_pur), c_pur


def _window_columns(key, mol, mf, shell_bands, shell_radials, ncol=None,
                    prior_bands=None, series_mfs=None):
    """The window object ``build_full_mo`` consumes: the raw AO rows,
    the distilled-shell columns when the case declares radials, or the
    PRIOR-AXIS columns when the window carries a source tag -- ONE
    resolution shared by the free run, the pinned lane, and the
    NO-start source window (start-independent by the same rule as the
    frozen masks: the distillation is the case's WINDOW DEFINITION,
    the mo/moc starts express their own gauge within its width).
    Prints the distillation report (the composition log's entry).

    Returns ``(act, ncas, report)``: row indices, the width, and the
    construction's disclosure (``None`` for raw-band windows -- a
    reporting construction is a compressed one:
    :func:`atomterm.windows.window_distilled` -- the radial signal is
    the OCCUPIED reference block; a source tag routes to
    :func:`atomterm.basis.series_window` -- the density-block family,
    ``@r`` the head configuration, ``@s`` the series average).  The
    report rides the payload (``window_report``) -- the k evidence
    (bands, spectrum, discarded fraction, captures) is store data,
    not console text.
    """
    from atomterm.windows import (realized_patterns, select_valence_ao,
                                  window_distilled)

    if key.source_tag is not None:
        src, k = key.source_tag
        if src == "s":
            if series_mfs is None:
                raise ValueError(
                    f"window {key.window!r} needs the configuration "
                    "series (the case's declared SeriesSpec tuple, "
                    "resolved through atomterm.basis.run_series) -- "
                    "pass series_mfs")
            mfs = list(series_mfs.values())
        else:
            mfs = [mf]
        if ncol is None:
            raise ValueError(
                f"window {key.window!r} needs the core boundary ncol "
                "(the case's core pair count -- the series blocks cut "
                "the valence there)")
        from atomterm.basis import series_window
        act, report = series_window(
            mol, mfs, mf.get_ovlp(), key.shell_patterns,
            prior_bands if prior_bands is not None else shell_bands,
            {sh: k for sh in key.shell_patterns}, ncol)
        print("  [window] prior " + "; ".join(
            f"{sh}: k={r['k']} of {len(r['bands'])} bands, "
            f"discarded {r['discarded']:.4f}, "
            f"capture {min(min(c) for c in r['capture']):.4f}"
            f"-{max(max(c) for c in r['capture']):.4f}"
            for sh, r in report.items()), flush=True)
        return act, act.shape[1], report
    act = select_valence_ao(
        mol, realized_patterns(key.shell_patterns, shell_bands))
    if shell_radials is None:
        return act, len(act), None
    c_occ = mf.mo_coeff[:, np.asarray(mf.mo_occ) > 0]
    act, report = window_distilled(
        mol, c_occ, mf.get_ovlp(), key.shell_patterns,
        shell_bands, shell_radials)
    print("  [window] distilled " + "; ".join(
        f"{sh}: k={r['k']} of {len(r['bands'])} bands, "
        f"capture {min(r['capture']):.4f}-{max(r['capture']):.4f}"
        for sh, r in report.items()), flush=True)
    return act, act.shape[1], report


def _window_ncas(key, mol, shell_bands, shell_radials) -> int:
    """The stored-absent ``ncas`` fallback: the distilled width when
    the case declares radials, the prior tag's width when the window
    carries one, else the raw band-row count."""
    from atomterm.windows import (distilled_ncas, realized_patterns,
                                  select_valence_ao)

    if key.source_tag is not None:
        _src, k = key.source_tag
        return distilled_ncas(key.shell_patterns,
                              {sh: k for sh in key.shell_patterns})
    if shell_radials is not None:
        return distilled_ncas(key.shell_patterns, shell_radials)
    return len(select_valence_ao(
        mol, realized_patterns(key.shell_patterns, shell_bands)))


def run_cell(key: CellKey, mf, mol=None, *, frac_mf=None, budgets: dict = None,
             fresh: bool = False, shell_bands: dict | None = None,
             shell_radials: dict | None = None,
             series_mfs: dict | None = None,
             prior_bands: dict | None = None,
             workdir: str | Path = "elements") -> "StageFacts | None":
    """Execute one cell: resolve every axis through its OWNER, run,
    store, analyze. The unit of the run queue.

    Parameters
    ----------
    key : the cell to run.
    mf : the atom's converged reference (REQUIRED -- reference
        construction, including character-forced occupations and the
        basin gate, belongs to cases/reference, not to cells). On
        ``core="frac"`` this is REPLACED by ``frac_mf``.  The
        reference must live on the KEY's field: a folded key takes
        the folded reference (:func:`atomterm.reference.run_folded`'s
        product) in place of the neutral one.
    mol : the analysis mol; defaults to ``mf.mol``. Only a
        ``core="frac"`` run needs it stated separately (the run basis
        is ``frac_mf.mol``).
    frac_mf : the fractional-occupation reference; required iff
        ``key.core == "frac"`` (that core is a different reference SCF
        -- the run AND the core use it).
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
        ``sa_weights`` explicitly there.  A declared block LONGER
        than the sector (the occupation-group windows: a CAS(3,2)
        sector cannot host a term's 2L+1) falls back to the measured
        rule -- window-local by construction.
    fresh : ignore a stored cell and re-run (a cache hit otherwise
        re-analyzes the stored run -- budgets are then ignored).
    shell_bands : the case's declared shell layer
        (:attr:`atomterm.cases.CaseSpec.shell_bands`) -- chemical
        shells resolve to their declared AO bands through it (None =
        the 1:1 identity, the contracted-basis behavior).
    shell_radials : the case's declared WINDOW distillation depth
        (:attr:`atomterm.cases.CaseSpec.shell_radials`; None = the
        raw band rows).  Start-independent: the ao/no starts carry
        the distilled columns, the mo/moc starts their own canonical
        gauge at the distilled width.  Undeclared combinations fail
        fast: the frac core and the radial mask have no distilled
        meaning (the fractional reference's open shell and the
        member-labeled radial substrate are row-space constructs).
    series_mfs : the configuration series' converged mfs
        (:func:`atomterm.basis.run_series`'s product; the HEAD is the
        case reference itself).  REQUIRED by ``@s``-tagged windows
        (the series source); reference construction belongs to the
        caller, never to cells.
    prior_bands : the prior axis's realization
        (:attr:`atomterm.cases.CaseSpec.prior_bands`) -- the widened,
        certificate-declared band lists the ``@r``/``@s`` windows
        distill from (None = the case's ``shell_bands``).
    workdir : the store root (``<workdir>/<Atom>/rungs/``).

    Returns
    -------
    The stage's :class:`~atomterm.stages.StageFacts`; persist with
    :func:`cell_payload` (which embeds the address).  ``None`` past
    the dense-meter budget (:data:`_COMM_CAP` dets): the cell stores
    kernel-only -- :func:`kernel_cell_payload` is the caller's
    continuation.
    """

    budgets = dict(budgets or {})
    if key.core == "frac":
        if frac_mf is None:
            raise ValueError("core='frac' needs the fractional-occupation "
                             "reference (frac_mf=...)")
        mf = frac_mf
        mol = mol if mol is not None else mf.mol
    mol = mf.mol if mol is None else mol

    # the route-class grammar first, before any store access: a
    # core-tagged window on the free route is ao/pol/avg (the pur
    # class is run_lane_cell's; raw/frac of a tagged window are
    # meaningless -- the tag is the partition, the core exactly
    # m-scalar)
    if key.core_tag is not None and (key.substrate, key.core,
                                     key.start) != ("pol", "avg", "ao"):
        raise ValueError(
            f"{key.address}: a free-route lane window is ao/pol/avg "
            "(the pur class is run_lane_cell's; the tag is the "
            "partition -- raw/frac do not apply)")

    chk = cell_chk(key, workdir)
    if chk.exists() and not fresh:
        return analyze_cell(key, mol, mf, shell_bands=shell_bands,
                            workdir=workdir)
    chk.parent.mkdir(parents=True, exist_ok=True)

    from atomterm.casscf import (constrained_casscf, mask_frozen_core,
                                 mask_frozen_window, mask_no_core_window_mix,
                                 run_casci)
    from atomterm.windows import (build_full_mo, canonical_order,
                                  occupied_core_idx, radial_basis,
                                  mask_radial_rotations, realized_patterns,
                                  select_valence_mo)

    na, nb = key.sector
    if (shell_radials is not None or key.source_tag is not None) and (
            key.core == "frac" or key.mask == "radial"):
        raise ValueError(
            f"{key.address}: the distilled/prior window is undeclared "
            f"for core={key.core}, mask={key.mask} -- the fractional "
            "reference's open shell and the member-labeled radial "
            "substrate are row-space constructs")
    ncore, c_core, c_pur = _key_core(key, mol, mf, shell_bands)
    act, ncas, win_report = _window_columns(
        key, mol, mf, shell_bands, shell_radials,
        ncol=ncore, prior_bands=prior_bands, series_mfs=series_mfs)

    # start x mask -> the substrate MO matrix + the optimizer mask.
    # frozen/virtmix are START-INDEPENDENT (they pin the window columns
    # whatever substrate they came from); radial needs the
    # member-labeled AO substrate and is ao-only (validated at the key).
    mask = None
    if key.start == "no":
        # the parent-prefix policy: the label names the SOURCE window,
        # the k-th occupation-group prefix of its term-averaged density
        # is the window (never splits an m-set; exactly closed when the
        # density was m-scalar -- natural.no_group_windows verifies).
        # The SOURCE rides the PURIFIED span whatever the cell's
        # substrate: the density is a measurement device (like the
        # distilled columns, substrate-independent), and only the
        # closed core gives the exact degenerate block an exactly
        # m-scalar average (measured: the canonical-core block splits,
        # and the split states' average is m-scalar only to the
        # admixture level); the cell's own substrate applies at the
        # final assembly below.
        from atomterm.core import purified_core_orbitals
        from atomterm.natural import no_group_windows

        c_src_core = (c_pur if c_pur is not None
                      else purified_core_orbitals(mol, mf, ncore)[0])
        mo_src = build_full_mo(mf, occupied_core_idx(mf, ncore), act,
                               c_core=c_src_core)
        c_src = mo_src[:, ncore:ncore + ncas]
        wins = no_group_windows(mol, mf, ncore, c_src, na, nb,
                                c_core=c_src_core,
                                cap=budgets.get("no_cap", 12000))
        win = _pick_no_window(wins, key.no_source[1])
        ncas = win.c_win.shape[1]
        mo_full = win.columns(mf, ncore, c_core=c_core)
    elif key.start == "ao":
        # a tagged core supplied the core columns (the run substrate's
        # prefix); otherwise the occupied prefix cores
        core_idx = (list(range(ncore)) if key.core_tag is not None
                    else occupied_core_idx(mf, ncore))
        if key.mask == "radial":
            mo_full, labels = radial_basis(
                mf, core_idx, act, c_core=c_core)
            mask = mask_radial_rotations(labels, ncore, ncas)
        else:
            mo_full = build_full_mo(mf, core_idx, act, c_core=c_core)
    elif key.start == "moc":
        # the character-selected valence window
        # (:func:`windows.select_valence_mo`: the SCF's own valence +
        # the virtual complement ranked by the named bands'
        # S-projection weight -- never energy order, the diffuse
        # echo). pol
        # rides the canonical matrix itself (S-orthonormal in the
        # core it was solved in); pur re-projects the window on the
        # pur core span (:func:`windows.build_full_mo`'s honest
        # re-orthogonalization). Honest expectation: the columns'
        # own m-delocalization survives (measured: raw-MO spans never
        # rotation-close) -- these cells measure that, they do not
        # fix it.
        # the character ranking runs on the REALIZED bands (the
        # layer): on uncontracted bases the raw patterns name
        # single primitives, not the declared families
        perm = select_valence_mo(
            mf, realized_patterns(key.shell_patterns, shell_bands),
            ncore)
        c_can = mf.mo_coeff[:, perm]
        if key.substrate == "pur":
            mo_full = build_full_mo(mf, occupied_core_idx(mf, ncore),
                                    c_can[:, ncore:ncore + ncas],
                                    c_core=c_pur)
        else:
            mo_full = c_can
    else:  # start == "mo": the canonical window in the run convention
        # (:func:`atomterm.windows.canonical_order`: closed | open |
        # virtual -- the raw PySCF order is the m-block listing and
        # would freeze virtual junk on uncontracted bases). On pol the
        # reordered canonical matrix IS the substrate; on pur the
        # canonical columns are re-orthogonalized against the pur core
        # span (they are S-orthonormal against the core they were
        # SOLVED in, not the pur one -- exactly the input-window defect
        # build_full_mo's projection removes at the source). Honest
        # expectation: the window's own m-delocalization survives
        # (measured: raw-MO spans never rotation-close) -- these
        # cells measure that, they do not fix it.
        perm = canonical_order(mf, ncore)
        c_can = mf.mo_coeff[:, perm]
        if key.substrate == "pur":
            mo_full = build_full_mo(mf, list(range(ncore)),
                                    c_can[:, ncore:], c_core=c_pur)
        else:
            mo_full = c_can
    if key.mask == "frozen":
        mask = mask_frozen_window
    elif key.mask == "virtmix":
        mask = mask_no_core_window_mix
    elif key.mask == "frozencore":
        mask = mask_frozen_core

    verbose = budgets.get("verbose", 0)
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
        if key.method == "sa-casscf" and weights is not None \
                and len(weights) > sector_dimension(ncas, na, nb):
            # the declared case block cannot fit this window (the
            # occupation-group class: a CAS(3,2) sector hosts 2
            # roots, not a term's 2L+1) -- the MEASURED block anchors
            # instead (the default rule, window-local by construction)
            weights = None
        if key.method == "sa-casscf" and weights is None:
            from atomterm.operators import build_operators, ground_block

            c_act0 = mo_full[:, ncore:ncore + ncas]
            ops0 = build_operators(mol, mf, c_act0, ncore, na, nb,
                                   c_core=c_core)
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
    return _facts(key, mol, mf, cas, ncore, ncas, na, nb, c_core=c_core,
                  window_report=win_report)


# analysis + persistence ---------------------------------------------

def analyze_cell(key: CellKey, mol, mf, shell_bands: dict | None = None,
                 shell_radials: dict | None = None,
                 workdir: str | Path = "elements"):
    """Re-analyze a stored cell (load + one facts pass) -- no mol/mf-bound
    optimizer re-run, no re-solving. The two-function lean: loading is
    pure data (:func:`load_cell`), analysis needs the live reference.

    Core-tagged windows analyze against their tag's SCF core
    (the core size is the core's PAIR count -- the span width on a
    multi-column realization); sectors past the dense-meter budget
    return None -- the kernel-only class
    (:func:`kernel_cell_payload`).  ``converged`` defaults True (a
    trusted cached stage).
    """
    from atomterm.casscf import CASResult

    data = load_cell(key, workdir)
    cas = CASResult(data.energies, data.ci, data.mo_coeff,
                    None, data.core_energy,
                    cycles=data.cycles or None,
                    conv_dE=data.conv_dE)
    na, nb = key.sector
    ncas = data.ncas or _window_ncas(key, mol, shell_bands, shell_radials)
    if sector_dimension(ncas, na, nb) > _COMM_CAP:
        return None
    ncore, c_core, _ = _key_core(key, mol, mf, shell_bands)
    return _facts(key, mol, mf, cas, ncore, ncas, na, nb, c_core=c_core)


def cell_payload(key: CellKey, facts) -> dict:
    """The persisted payload of a cell: the stage facts + the address."""
    out = facts.payload(key.method)
    out["cell"] = key.address
    return out


def _apply_family(mol, mf, c_win, parts, ncore, c_core, ncas, na, nb):
    """One home for the size-independent apply family of a window:
    ``(apply_h, apply_m, apply_l2, apply_s2)`` -- H through the PySCF
    C-contraction boundary (:func:`atomterm.operators.pyscf_apply`),
    Lz as its real form ``M = -i Lz`` (norms unchanged), L2 with its
    two-body part, S2 one-body.  Shared by :func:`kernel_cell_payload`
    and :func:`run_lane_cell`; the callers own their root loops."""
    from pyscf import ao2mo

    from atomterm.detbasis import operator_apply
    from atomterm.operators import active_space_integrals, pyscf_apply

    h1, _ = active_space_integrals(mol, mf, c_win, ncore, c_core=c_core)
    apply_h = pyscf_apply(h1, ao2mo.kernel(mol, c_win), ncas, na, nb)
    apply_m = operator_apply(parts["dets"],
                             one_body=(-1j * parts["Lz"][0]).real)
    apply_l2 = operator_apply(parts["dets"], one_body=parts["L2"][0],
                              two_body=parts["L2"][1])
    apply_s2 = operator_apply(parts["dets"], one_body=parts["S2"][0],
                              two_body=parts["S2"][1])
    return apply_h, apply_m, apply_l2, apply_s2


def kernel_cell_payload(key: CellKey, mol=None, mf=None,
                        shell_bands: dict | None = None,
                        shell_radials: dict | None = None,
                        workdir: str | Path = "elements") -> dict:
    """The kernel-only payload of a big-sector cell (run_cell returns
    None past the dense-meter budget): the run's own energies from the
    store.  With the substrate given (``mol``/``mf``), the payload
    ALSO carries the apply-based meters -- the ground root's
    ``ground_l2`` and the state residual ``res_h_lz``,
    :func:`run_lane_cell`'s size-independent certificates -- filling
    the columns the dense budget leaves NaN (the tables' em-dash
    class).  The meters read the STORED window (the run's own
    converged MO) through the same core resolution :func:`run_cell`
    uses (``shell_bands`` -- the case's declared shell layer); the
    kernel totals are the values that compare."""
    data = load_cell(key, workdir)
    out = {
        "cell": key.address,
        "kind": key.method,
        "n_det": int(sector_dimension(data.ncas, *key.sector)),
        "energies_pyscf": [float(x) for x in np.atleast_1d(data.energies)],
        "ground_l2": float("nan"),
        "comm_norm_h_lz": float("nan"),
        "comm_norm_h_l2": float("nan"),
        "converged": True,
        "cycles": data.cycles or None,
        "max_cycle": data.max_cycle,
        "conv_dE": data.conv_dE,
        "win_lz": float("nan"),
        "win_lperp": float("nan"),
    }
    if mol is None or mf is None or key.core == "frac":
        # frac cells ran under the fractional reference (not this mf);
        # their meters belong to the dense path
        return out
    from atomterm.operators import (build_operator_parts,
                                    ci_in_detbasis_order, window_closure)

    na, nb = key.sector
    ncore, c_core, _ = _key_core(key, mol, mf, shell_bands)
    ncas = data.ncas or _window_ncas(key, mol, shell_bands, shell_radials)
    c_act = np.asarray(data.mo_coeff)[:, ncore:ncore + ncas]
    parts = build_operator_parts(mol, mf, c_act, ncore, na, nb,
                                 c_core=c_core)
    # the one-particle closure certificate of the stored window (the
    # span projector's commutators -- cheap at any sector size)
    closure = window_closure(mol, c_act)
    out["win_lz"] = float(closure["lz"])
    out["win_lperp"] = float(closure["lperp"])
    k = int(np.argmin(np.atleast_1d(data.energies)))
    w = ci_in_detbasis_order(np.asarray(data.ci)[k], ncas, na, nb)
    apply_h, apply_m, apply_l2, _apply_s2 = _apply_family(
        mol, mf, c_act, parts, ncore, c_core, ncas, na, nb)
    out["ground_l2"] = float(np.real(w @ apply_l2(w)))
    out["res_h_lz"] = float(np.linalg.norm(
        apply_h(apply_m(w)) - apply_m(apply_h(w))))
    out["res_h_l2"] = float(np.linalg.norm(
        apply_h(apply_l2(w)) - apply_l2(apply_h(w))))
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
    cycles = getattr(cas, "cycles", None)
    if cycles is not None:
        lib.chkfile.save(str(chk), f"atomterm/{key.address}/cycles",
                         int(cycles))
        # the cap rides the count: the cap-hitting flag compares by
        # identity in the payloads (never a bare literal)
        cap = getattr(getattr(cas, "mc", None), "max_cycle_macro", None)
        if cap is not None:
            lib.chkfile.save(str(chk), f"atomterm/{key.address}/max_cycle",
                             int(cap))
    conv_dE = getattr(cas, "conv_dE", None)
    if conv_dE is not None:
        lib.chkfile.save(str(chk), f"atomterm/{key.address}/conv_dE",
                         float(conv_dE))


def _facts(key: CellKey, mol, mf, cas, ncore: int, ncas: int,
           na: int, nb: int, c_core=None, window_report=None):
    """The cell layer's one facts call: the core's spherical convention
    and the CASCI anisotropy report are the only cell-specific
    choices; everything else is :func:`atomterm.stages.stage_facts`
    (``window_report`` carries the construction's disclosure through
    to the payload)."""
    from atomterm.stages import stage_facts

    c_act = cas.mo_coeff[:, ncore:ncore + ncas]
    return stage_facts(
        mol, mf, cas, c_act, ncore, na, nb, tag=key.address,
        spherical_core=key.core != "raw",
        # every CASCI cell reports its reference core's anisotropy (the
        # raw core's is the detector; the pur core's is the ~0
        # certificate; row parity across the core block)
        anisotropy_ref=(mf, c_core) if key.method == "casci" else None,
        window_report=window_report,
    )


# the pinned lane ------------------------------------------------------


# the payload's required fields, in store order -- the ONE field list
# shared by the writer (:func:`run_lane_cell`) and the reader
# (:func:`load_lane_cell`); method-dependent extras are optional
_LANE_FIELDS = ("term", "energies", "total", "ptp", "l2", "s2",
                "n_det", "anchor_dev", "win_lz", "win_lperp",
                "comm_norm_h_lz", "comm_norm_h_l2", "res_h_lz",
                "res_h_l2", "fold")
# optional fields, loaded when present: the SA ensemble and, on
# optimizer cells, the converged FULL MO matrix (restart +
# orbital-level route comparisons) with the optimizer's cycle count;
# the frozen-core certificate (the SCF core's minimum in-span weight;
# the pre-repair stores carry the retired confined-core Roothaan's
# cycles/comm instead)
_LANE_FIELDS_OPT = ("ensemble", "mo_coeff", "cycles", "conv_dE",
                    "core_pairs", "core_cycles", "core_comm",
                    "core_in_span_min", "anchor_corrected")


def _root_vectors(ci):
    """A solver's return as a list of root vectors (nroots = 1 returns
    a bare array, nroots > 1 a stacked one)."""
    ci = np.asarray(ci)
    return [ci] if ci.ndim == 2 else list(ci)


def run_lane_cell(key: CellKey, mol, mf, term, conv_tol: float = 1e-12,
                  anchor: bool = False, fresh: bool = False, verbose: int = 0,
                  shell_bands: dict | None = None,
                  shell_radials: dict | None = None,
                  workdir: str | Path = "elements") -> dict:
    """Run one pinned-lane cell: the lane family's UNIFORM construction.
    The core resolves through the shell layer to its raw AO-shell span
    (Loewdin-orthonormalized), and the FROZEN CORE is the reference's
    OWN core -- the occupied columns S-projected onto that span
    (:func:`_lane_core`; the field and the bookkeeping join the
    neutral anchor BY CONSTRUCTION; ``ncore`` follows the pair count,
    never the span width); the
    window placed by :func:`atomterm.windows.build_full_mo` against
    it; one ``wfnsym``-pinned solver per branch component
    (:func:`atomterm.irreps.branch_d2h`); meters per irrep block
    (L2/S2/H are totally symmetric, so the block restriction is
    exact); the core-independent bookkeeping ruler
    (:func:`atomterm.core.closed_shell_core_energy` over the frozen
    columns -- the projected-density form).  Deterministic: raw AO
    spans on both sides of the core.

    Parameters
    ----------
    key : the lane cell (a :func:`lane_keys` key: ``casci/ao/-/pur/avg``
        with a tagged window).  Its sector must hold exactly the
        electrons its core leaves active.
    mol, mf : the PINNED substrate -- ``mol`` a converged-reference
        D2h molecule, ``mf`` its ROHF.  The substrate lives on the
        KEY's field: a folded key takes the folded reference
        (:func:`atomterm.reference.run_folded`'s product, its own D2h
        clone) in place of the neutral one.
    term : the target ``(L, parity)`` (:attr:`~atomterm.cases.CaseSpec.term`);
        parity is the configuration's total inversion parity.
    conv_tol : the pinned solvers' CI tolerance.
    anchor : also eigen-solve each block PAST the dense-adjudication
        budget with OUR H and record the max deviation from the
        pinned solver (the conventions certificate; O(block^3) --
        keep deliberate).  Blocks UNDER that budget adjudicate
        automatically, deviation correcting the roots (recorded in
        ``anchor_corrected``).
    fresh : ignore a stored payload and re-run.
    shell_bands : the case's declared shell layer (chemical shell ->
        its AO bands; None = the 1:1 identity).
    shell_radials : the case's declared WINDOW distillation depth
        (:attr:`atomterm.cases.CaseSpec.shell_radials`; None = the raw
        band rows) -- the pinned term readout at the chemically-sized
        window; the CORE stays at full family width either way.
    verbose : print each block as it solves.
    workdir : the store root.

    Returns
    -------
    the payload dict (also stored under the key's address, reloadable
    through :func:`load_lane_cell`): per-component ``energies`` /
    ``l2`` / ``s2`` in branch order, the TERM ``total`` (component
    mean + the closed-shell core energy), the component spread
    ``ptp``, the sector size ``n_det``, ``anchor_dev``,
    the dense-adjudication record (``anchor_corrected`` -- the blocks
    whose roots our exact H replaced or re-selected: the pinned
    Davidson's near-degenerate skips, and the targeted selection when
    the declared term is not the block's ground), the dense
    assembled-treatment commutators ``comm_norm_h_lz`` /
    ``comm_norm_h_l2`` (NaN above :data:`_COMM_CAP` -- sectors past
    the dense-meter budget; the certified core class makes them
    structural there).  SA cells add the ``ensemble``, the converged
    FULL MO ``mo_coeff`` (nested lists in the payload -- the
    JSON-store contract; :func:`load_lane_cell` returns the matrix),
    and the optimizer's ``cycles``.

    Raises
    ------
    ValueError
        the key is not a lane cell; the core construction refuses the
        key (:func:`_lane_core`: window/core overlap, a sector
        disagreeing with the core-derived electron count, an absent
        core shell).
    """
    from pyscf import lib

    from atomterm.core import closed_shell_core_energy
    from atomterm.irreps import (IRREP_IDS, branch_d2h, det_irreps,
                                 orb_irreps)
    from atomterm.operators import (build_operator_parts,
                                    ci_in_detbasis_order,
                                    window_closure)
    from atomterm.detbasis import sector_dimension, slater_condon_matrix
    from atomterm.windows import build_full_mo

    if (key.method, key.start, key.mask, key.substrate,
            key.core) not in (("casci", "ao", None, "pur", "avg"),
                              ("sa-casscf", "ao", None, "pur", "avg")):
        raise ValueError(
            f"a lane cell is casci/sa-casscf at ao/-/pur/avg (the "
            f"certified core class); {key.address} is not one")
    chk = cell_chk(key, workdir)
    if chk.exists() and not fresh:
        return load_lane_cell(key, workdir)
    chk.parent.mkdir(parents=True, exist_ok=True)

    L, parity = term
    na, nb = key.sector
    core, nf, c_core, core_report = _lane_core(key, mol, mf, shell_bands)
    act, ncas, _ = _window_columns(key, mol, mf, shell_bands, shell_radials)
    mo = build_full_mo(mf, list(range(nf)), act, c_core=c_core)
    c = mo[:, nf:nf + ncas]

    # the kernel half: the pinned solvers live in the evaluator (one
    # per branch component, multiplicity as nroots) -- SA returns the
    # ensemble plus the exact roots re-solved in the frozen converged
    # window (the two-pass pattern); CI solves the start window
    es, l2s, s2s = [], [], []
    anchor_dev = 0.0
    ensemble = None
    extra: dict = {}
    if key.method == "sa-casscf":
        from atomterm.casscf import run_irrep_sa

        out = run_irrep_sa(mol, mf, c, nf, na, nb, L, parity,
                           conv_tol=conv_tol, c_core=c_core)
        ensemble = float(out["e_tot"])
        roots = out["roots"]
        # the converged window PERSISTS (restart + orbital-level route
        # comparisons) -- the run's own partition moves under
        # same-irrep rotations, the span certificate survives.  The
        # payload dict is the JSON-store contract (regen dumps it
        # verbatim): mo_coeff travels as nested lists, and
        # load_lane_cell hands the matrix back.
        extra = {"mo_coeff": np.asarray(out["mo"], dtype=float).tolist(),
                 "cycles": int(out["cycles"]),
                 "conv_dE": out["conv_dE"]}
        c_read = out["mo"][:, nf:nf + ncas]
    else:
        from atomterm.casscf import run_irrep_casci

        if verbose:
            print(f"  lane {key.address} pinned blocks "
                  f"{list(branch_d2h(L, parity))}")
        roots = run_irrep_casci(mol, mf, c, nf, na, nb, L, parity,
                                conv_tol=conv_tol, c_core=c_core)
        c_read = c

    # the readout half (shared), on the READOUT window (CI: the start
    # window; SA: the converged one): per component block, the pinned
    # roots bridge into detbasis order and live entirely in their
    # block (the support assert); the meters read them restricted to
    # the block -- L2/S2/H are totally symmetric, so the restriction
    # is exact.  The SA TOTAL is the kernel ensemble (the optimizer
    # may rotate same-irrep core/window columns -- the readout's
    # absolute energies are of the run's own partition; its COMPONENT
    # SPREAD and meters are the diagnostics).
    orbsyms = orb_irreps(mol, c_read)
    parts = build_operator_parts(mol, mf, c_read, nf, na, nb,
                                 c_core=c_core)
    dets = parts["dets"]
    ids = det_irreps(orbsyms, na, nb)
    # the apply family FIRST: the per-root meters (L2/S2 expectations,
    # the state residuals) ride it -- matrix-free and exact at any
    # sector size, where the per-block dense meter matrices would be
    # O(block^2) memory
    apply_h, apply_m, apply_l2, apply_s2 = _apply_family(
        mol, mf, c_read, parts, nf, c_core, ncas, na, nb)
    w_all: list[np.ndarray] = []
    corrected: list[str] = []
    for irrep, mult in branch_d2h(L, parity).items():
        e_arr, ci = roots[irrep]
        e_arr = np.atleast_1d(np.asarray(e_arr, dtype=float))
        sel = ids == IRREP_IDS[irrep]
        block = [d for d, m in zip(dets, sel) if m]
        # the DENSE ADJUDICATION (every CI block the budget fits): the
        # pinned Davidson can deterministically SKIP near-degenerate
        # lowest roots (measured on the Sc-unc lane: the Ag ground
        # PAIR of the 2D term -- the solver's second root landed four
        # eigenvalues up; on Ti-unc ALL FOUR blocks skipped), so our
        # own exact H restricted to the block is the ADJUDICATOR: it
        # certifies the pinned roots (anchor_dev) and, on deviation,
        # REPLACES them (the corrected blocks recorded; the vectors
        # rebuild from the dense eigenvectors).  The external anchor
        # stands -- both numbers are kept; the lane's truth is the
        # exact construction, the solver's deviation the recorded
        # diagnostic.
        ours_vecs = None
        if key.method == "casci" and len(block) <= _ANCHOR_CAP:
            hb = slater_condon_matrix(block,
                                      one_body=parts["H"][0],
                                      two_body=parts["H"][1])
            ev, evec = np.linalg.eigh(hb)
            # the TARGETED selection: the declared term need not be the
            # block's ground (measured on the U14e port: an L=6 quintet
            # sits BELOW the declared L=8 in the raw band window -- the
            # block-lowest readout interleaves the terms).  The dense
            # spectrum holds every root; the term meter (<L2> per
            # eigenvector, cheap through the apply family) SELECTS the
            # block's lowest `mult` eigenvectors of the TARGET L
            # (falling back to the lowest when the target is not
            # found -- the payload's l2 then shows it)
            want = min(len(ev), _SELECT_DEPTH * mult)
            idx = np.flatnonzero(sel)
            l2_of = []
            for k in range(want):
                wv = np.zeros(len(dets))
                wv[idx] = np.real(evec[:, k])
                l2_of.append(float(np.real(wv @ apply_l2(wv))))
            l2_of = np.asarray(l2_of)
            hit = np.flatnonzero(np.abs(l2_of - L * (L + 1)) < 1e-4)
            pick = (hit[:mult] if len(hit) >= mult
                    else np.arange(mult))          # the fallback
            dev = float(np.max(np.abs(np.sort(e_arr)[:mult]
                                      - ev[pick])))
            anchor_dev = max(anchor_dev, dev)
            if dev > 1e-8 or not np.array_equal(pick, np.arange(mult)):
                corrected.append(irrep)
                e_arr = ev[pick].copy()
                ours_vecs = evec[:, pick]
            del hb
        es.extend(float(x) for x in e_arr[:mult])
        w_blk = []
        if ours_vecs is None:
            for vec in _root_vectors(ci):
                w = ci_in_detbasis_order(vec, ncas, na, nb)
                # the pinned vector lives entirely in its block
                assert abs(np.linalg.norm(w[sel]) - 1.0) < 1e-8, irrep
                w_blk.append(w)
        else:
            idx = np.flatnonzero(sel)
            for k in range(ours_vecs.shape[1]):
                w = np.zeros(len(dets))
                w[idx] = np.real(ours_vecs[:, k])
                w_blk.append(w)
        w_all.extend(w_blk)
        if anchor and key.method == "casci" and len(block) > _ANCHOR_CAP:
            # the opt-in anchor for blocks past the dense-adjudication
            # budget (O(block^3) -- keep deliberate)
            hb = slater_condon_matrix(block,
                                      one_body=parts["H"][0],
                                      two_body=parts["H"][1])
            ours = np.linalg.eigvalsh(hb)[:mult]
            anchor_dev = max(anchor_dev,
                             float(np.max(np.abs(np.sort(e_arr) - ours))))
            del hb
        for w in w_blk:
            # the meters read the roots through the apply family
            # (matrix-free, size-independent)
            l2s.append(float(np.real(w @ apply_l2(w))))
            s2s.append(float(np.real(w @ apply_s2(w))))

    # the state residual ||[H, Lz]psi|| over the component roots -- the
    # SIZE-INDEPENDENT symmetry certificate (exact 0 iff a symmetry
    # eigenstate; the operator Frobenius norm it backs up past the
    # dense budget needs every matrix element).  Lz enters as its real
    # form M = -i Lz; H through the C contraction boundary on the
    # readout window's own feed (:func:`_apply_family` is the ONE
    # construction -- built once above, shared with the meters).
    res_lz = 0.0
    for w in w_all:
        r = apply_h(apply_m(w)) - apply_m(apply_h(w))
        res_lz = max(res_lz, float(np.linalg.norm(r)))
    # the L2 residual: L2 carries a two-body part, so its closure pays
    # the pair classification once (minutes-class at the 2e5-det
    # sector, cached after) -- priced by measurement, not declined
    res_l2 = 0.0
    for w in w_all:
        r = apply_h(apply_l2(w)) - apply_l2(apply_h(w))
        res_l2 = max(res_l2, float(np.linalg.norm(r)))

    es = np.asarray(es)
    # the assembled-treatment meters: exact dense commutators where
    # the sector fits memory; above the cap the certified core class
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
        total = float(es.mean() + closed_shell_core_energy(mf, c_core))
    # the one-particle closure certificate of the READOUT window (the
    # span projector's commutators -- gauge-invariant, any size)
    closure = window_closure(mol, c_read)
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
        "win_lz": float(closure["lz"]),
        "win_lperp": float(closure["lperp"]),
        "comm_norm_h_lz": comm_lz,
        "comm_norm_h_l2": comm_l2,
        "res_h_lz": float(res_lz),
        "res_h_l2": float(res_l2),
        "fold": list(core),
        # the frozen-core certificates: the pair count the core
        # freezes, and the SCF core's in-span weights (the span
        # hosts the core shells' full family -- the minimum weight
        # is the restriction certificate)
        "core_pairs": int(nf),
        "core_in_span_min": float(min(core_report["in_span"])),
        # the blocks whose pinned-solver roots our dense H replaced
        # (the Davidson's near-degenerate root skips; empty = agreement)
        "anchor_corrected": corrected,
    }
    if ensemble is not None:
        payload["ensemble"] = ensemble
    payload.update(extra)
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
    for field_name in _LANE_FIELDS:
        value = lib.chkfile.load(str(chk),
                                 f"atomterm/{key.address}/{field_name}")
        if value is None:
            raise KeyError(
                f"lane store {chk} is incomplete (no {field_name}); "
                "regenerate the cell")
        out[field_name] = _stored(value)
    # optional fields (method-dependent): the SA ensemble, the
    # converged MO matrix, the optimizer's cycles and conv_dE --
    # mo_coeff loads as the array it was stored as (consumers do
    # linear algebra)
    for field_name in _LANE_FIELDS_OPT:
        value = lib.chkfile.load(str(chk),
                                 f"atomterm/{key.address}/{field_name}")
        if value is not None:
            out[field_name] = (np.asarray(value, dtype=float)
                               if field_name == "mo_coeff"
                               else _stored(value))
    return out


# the declared key-sets ----------------------------------------------

def spine_keys(atom: str, window: str, sector: tuple[int, int]) -> list[CellKey]:
    """The default cell subset per atom (README · Methods, "the spine"):
    the three cores; frozen x purified under both objectives (the
    certified class); one free row per substrate (the basin detectors).
    Diagnostics beyond these run as needed. A SUBSET view of the grid
    (every spine key is in :func:`window_grid`'s product); pinned here
    because README quotes it verbatim.
    """
    sec = tuple(sector)
    return [
        CellKey(atom, window, sec, core="raw"),
        CellKey(atom, window, sec, core="avg"),
        CellKey(atom, window, sec, core="frac"),
        CellKey(atom, window, sec, method="casscf", mask="frozen",
                substrate="pur"),
        CellKey(atom, window, sec, method="sa-casscf", mask="frozen",
                substrate="pur"),
        CellKey(atom, window, sec, method="casscf"),
        CellKey(atom, window, sec, method="sa-casscf",
                substrate="pur"),
    ]


def no_label(window: str) -> str:
    """The full-span occupation-group prefix label of a window: k = its
    shell count ({3d,4s} -> the [s,d] groups -> ``3d4sno2``; each shell's
    m-set is one occupation group of the term-averaged density)."""
    return f"{window}no{len(re.findall(r'\d+[spdfgh]', window))}"


def window_grid(spec, window: str | None = None) -> list[CellKey]:
    """One window's full factor grid, in canonical render order.

    Parameters
    ----------
    spec : the :class:`~atomterm.cases.CaseSpec` (duck-typed: key,
        shell_patterns, nelec_act, spin, planned_sectors are read).
    window : the window label; None = the case's primary window.

    The ORDER is presentation and lives HERE, once -- renderers receive
    this list ready-made and format (branching on key identity only):
    the CASCI core block first (identity: method=casci, start=ao; the
    pur core third), each ruler core (avg, pur) carrying the COMPLETE
    KIND BLOCK -- mo, moc, no beside ao (the window-source comparison
    at fixed window and core; raw/frac are core-field controls and
    stay ao-only -- an m-incomplete start there would confound two
    diseases in one row) -- then the optimizer product method x mask x
    substrate x start, and the NO start addressing the full-span
    group prefix (:func:`no_label`). mo x pur and moc x pur ARE part
    of the grid (the canonical/character-selected columns are
    re-orthogonalized against the pur span); radial stays ao-only
    (validated at the key).  A SOURCE-TAGGED window drops the frac
    core (structurally: the prior family has no fractional-reference
    class -- it reads the series' own occupations, and a fractional
    reference is a different density altogether).
    """
    from atomterm.system import sector_na_nb

    w = window or "".join(spec.shell_patterns)
    atom = spec.key
    sec = (spec.planned_sectors.get(w)
           or sector_na_nb(spec.nelec_act, spec.spin))
    tagged = _SOURCE_TAG_RE.fullmatch(w.partition("@")[2]) is not None \
        if "@" in w else False
    keys = [
        CellKey(atom, w, sec, core="raw"),
        CellKey(atom, w, sec, core="avg"),
        CellKey(atom, w, sec, substrate="pur"),   # the CASCI pur core
    ]
    if not tagged:
        keys.append(CellKey(atom, w, sec, core="frac"))
    # the COMPLETE KIND BLOCK at each ruler core: the window-source
    # comparison (mo/moc/no beside ao) at fixed window and core --
    # raw/frac are core-field controls and stay ao-only (an
    # m-incomplete start there would confound two diseases in one
    # row)
    keys += [
        CellKey(atom, w, sec, method="casci", start="mo", core="avg"),
        CellKey(atom, w, sec, method="casci", start="moc", core="avg"),
        CellKey(atom, no_label(w), sec, method="casci", start="no",
                core="avg"),
        CellKey(atom, w, sec, method="casci", start="mo", substrate="pur"),
        CellKey(atom, w, sec, method="casci", start="moc",
                substrate="pur"),
        CellKey(atom, no_label(w), sec, method="casci", start="no",
                substrate="pur"),
    ]
    for m in ("casscf", "sa-casscf"):
        for mask in (None, *spec.optimizer_masks):
            for sub in ("pol", "pur"):
                starts = (("ao", "mo", "moc", "no") if mask != "radial"
                          else ("ao",))
                for st in starts:
                    keys.append(CellKey(atom, no_label(w) if st == "no"
                                        else w, sec, method=m, start=st,
                                        mask=mask, substrate=sub))
    return keys


def grid_keys(spec) -> list[CellKey]:
    """The case's runnable set: the primary window's full grid
    (:func:`window_grid`) -- the spine cores, the optimizer product, the
    CASCI pur core, and the virtual-mix diagnostic; every key distinct
    by construction (the spine is a subset view). What regen executes;
    planned windows are :func:`planned`."""
    return window_grid(spec)


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
    -- a split ground block groups by noise and the label misstates
    the window; the NO window variant lives in the methods grid,
    where its donor is declared per key.

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
    from dataclasses import replace
    keys = []
    for w in spec.lattice_windows:
        # ONE grammar home: the full factor grid (:func:`window_grid`)
        # over a sector view -- the lattice spans carry the same
        # structure as the primary (core groups, kinds, treatments),
        # their pending placeholders holding the unrun places
        view = replace(spec, planned_sectors={w: spec.lattice_sector})
        # frac stays out (a different reference SCF -- never a
        # lattice arm)
        keys.extend(k for k in window_grid(view, w) if k.core != "frac")
    return keys


def lane_keys(spec) -> list[CellKey]:
    """The declared pinned-lane legs (``spec.lane_windows``): one
    CI-level lane cell per leg, run by :func:`run_lane_cell` (the legs
    need ``spec.term`` declared -- validated at the CaseSpec).

    The SECTOR derives from the leg's own core: activating a shell
    moves its electrons into the window (the 3d-metal semi core leaves
    ``nelec_act + 4`` active), so a lane key never reuses the case's
    default sector blindly.

    Examples
    --------
    Ti (nelec_act=4, Ar prefix): ``3d4s@ar`` keeps sector (3, 1);
    ``3d3p@semi`` grows to (5, 3); the M-shell core lands at (6, 4).
    """
    from dataclasses import replace

    from atomterm.system import sector_na_nb
    from atomterm.windows import closed_electrons

    base = sector_na_nb(spec.nelec_act, spec.spin)
    keys = []
    for label in spec.lane_windows:
        k = CellKey(spec.key, label, base, method="casci",
                    substrate="pur", core="avg")
        nelec = (spec.nelec_act + 2 * spec.ncore
                 - closed_electrons(k.core_shells()))
        keys.append(replace(k, sector=sector_na_nb(nelec, spec.spin)))
    return keys


def route_keys(spec) -> list[CellKey]:
    """The ROUTE AXIS over the lane windows: every pinned CI leg's
    twins -- the free CI (a survey row; same window, same core, the
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

    legs = lane_keys(spec)
    out = list(legs)
    for leg in legs:
        for method in ("casci", "casscf", "sa-casscf"):
            out.append(replace(leg, method=method, substrate="pol"))
        out.append(replace(leg, method="sa-casscf"))
    return out


def radial_refused(key, spec) -> bool:
    """The construction refusal at the key: the radial mask's pair
    predicate needs member-labeled basis rows, and a prior window
    (``@r``/``@s``) or a distilled substrate (the case's declared
    ``shell_radials``) has no row identity -- the cell can never run
    on this case, whatever its scenario standing.  One home for the
    predicate: regen's queue skips on it (nothing stored, no place
    held), the renderer's skeleton consults it, and :func:`run_cell`
    keeps the runtime fail-fast as the deep guard."""
    return (key.mask == "radial"
            and (key.source_tag is not None
                 or spec.shell_radials is not None))


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
        # gates() owns the enabled check -- no second copy here
        if s.gates(a):
            return f"{s.name}: {s.reason}"
    return None
