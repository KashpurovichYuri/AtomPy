"""The declaration layer: one spec per atom case, the campaign registry,
and the launch entry that executes a declaration.

A case is DECLARED, not guessed: a :class:`CaseSpec` fixes the atom, the
reference data, the valence window, and the sectors -- the intent
manifest. :func:`run_case` executes the declaration
(the term tables are the contract, always produced); nothing
auto-generates. ``CAMPAIGN`` is the
per-atom registry the cell layer and the runners derive from --
nothing else re-declares atoms, windows, or the gated references (the
one-declaration rule).

Division of labor
-----------------
With the cell layer: this module declares WHAT is studied and runs
the classic staged pipeline (ROHF -> CASCI per Ms
sector -> CASSCF -> SA-CASSCF, checkpointed under stage tags);
:mod:`atomterm.cells` addresses HOW each single calculation of the
methods grid runs -- one frozen :class:`~atomterm.cells.CellKey` per
axis combination, its grids derived from these specs. The campaign's
two run vocabularies (stage tags here, cell addresses there) are
DELIBERATE: the term tables are a stable contract, the grid an
evolving experiment -- coupling them would put the committed tables'
reproducibility at risk.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from pyscf import gto

from atomterm.detbasis import sector_dimension
from atomterm.system import AtomConfig, ncore_of, sector_na_nb
from atomterm.basis import SeriesSpec
from atomterm.cells import no_label, prior_tag_of
from atomterm.reference import ReferenceSpec, run_rohf
from atomterm.scans import ScanSpec
from atomterm.windows import (build_full_mo, closed_electrons, core_shells,
                              occupied_core_idx, select_valence_ao)
from atomterm.casscf import (
    run_casci,
    constrained_casscf,
)
from atomterm.stages import stage_facts
from atomterm.analysis import format_report

__all__ = ["CaseSpec", "Scenario", "CAMPAIGN", "campaign_spec",
           "scenarios_for", "run_case", "run_atom", "run_carbon"]


# the declaration ----------------------------------------------------

@dataclass(frozen=True)
class Scenario:
    """One declared scenario: a named, dated physics claim that a set
    of cells need not run.

    A scenario is DECLARED DATA, interpreted by exactly one consumer
    -- the run gate (:func:`atomterm.cells.void_reason`: a scenario'd
    cell is void-with-reason at its key, still addressable; a full
    regen materializes it).  The vocabulary is generic on purpose:
    what a scenario MEANS physically lives in its ``reason``.  (The
    render half the original spec declared never came -- presentation
    lives with the emitters.)

    Attributes
    ----------
    name : the scenario's provenance key (campaign defaults and
        per-atom declarations merge on it).
    cells : cell ADDRESS PREFIXES (:attr:`atomterm.cells.CellKey.address`)
        -- a cell is scenario'd iff its address starts with one.
    reason : the dated physics claim the gate carries ("... (2026-11-02)");
        ``--full`` re-tests it.
    enabled : the whole scenario's gate.
    """

    name: str
    cells: tuple[str, ...]
    reason: str
    enabled: bool = True

    def __post_init__(self):
        if not self.name:
            raise ValueError("a scenario needs a name")
        if not self.cells:
            raise ValueError(f"scenario {self.name!r} names no cells")
        if not self.reason:
            raise ValueError(f"scenario {self.name!r} carries no reason")

    def gates(self, address: str) -> bool:
        """Whether ``address`` is this scenario's cell (prefix match)."""
        return self.enabled and address.startswith(self.cells)


# campaign-level scenarios applied to EVERY atom unless its row
# overrides them by name (none today -- the merge is the tested seam;
# a campaign-wide scenario arrives with its first campaign-wide claim)
DEFAULT_SCENARIOS: dict[str, Scenario] = {}


def _ms_label(na: int, nb: int) -> str:
    """Filename/report label of a sector's M_S = (na - nb) / 2.

    Whole numbers stay bare ("1", "0" -- matching the committed chkfile
    names); halves insert a "p" so paths stay clean ("0p5" for a doublet).
    """
    assert na >= nb, "sectors are built with na >= nb"
    two_ms = na - nb
    if two_ms % 2 == 0:
        return str(two_ms // 2)
    return f"{two_ms / 2:.1f}".replace(".", "p")


@dataclass(frozen=True)
class CaseSpec:
    """The declaration of one atom case -- what to run and what not.

    The term tables (T1/T2) are the contract: every executed case
    produces them. The single-root and SA stages are switched by their
    own flags.

    Attributes
    ----------
    key : the campaign key -- the ``elements/<Key>/`` dir name and the
        cell addresses' atom token (defaults to the symbol; a declared
        ECP atom carries its own, e.g. ``"U14e"``).  Must match the
        spec's ``CAMPAIGN`` registry key.
    symbol / basis / ecp / spin : the atom (:class:`~atomterm.system.AtomConfig`
        inputs; ``basis`` is a label like ``"ccpvdz"`` or a dict).
    nelec_act / ncore : the active-electron count and the closed-shell
        count (the core).
    shell_patterns : the valence-shell label prefixes of the standard
        window (e.g. ``["2p"]``; ``["3d", "4s"]``).
    shell_bands : the declared SHELL LAYER -- chemical shell -> the
        AO bands realizing it on this basis (``{"3d": ("12d", "13d")}``);
        every shell/window/core resolution goes through it
        (:func:`atomterm.windows.shell_realization`).  None = the 1:1
        identity (each shell IS its own band -- the contracted-basis
        behavior).  The realization is BASIS-SPECIFIC and certified by
        the measured composition (``scripts/band_composition.py``):
        a chemical shell over multi-zeta contractions or uncontracted
        primitive bands is multi-column, and its CORE freezes
        ``closed_electrons`` pairs over that span
        (:func:`atomterm.core.confined_core_orbitals`) -- never one
        pair per column.
    shell_radials : the declared WINDOW distillation depth (chemical
        shell -> k radial directions; ``None`` = no distillation, the
        raw band rows).  A multi-band realization's chemically-sized
        window: each shell's k radial directions distilled from the
        reference's occupied MOs by the member-summed radial Gram,
        replicated on every m-member
        (:func:`atomterm.windows.window_distilled`) -- the CAS label
        reads the delivered records' shapes (CAS(3,6), not CAS(3,43)).
        Window-side only: the CORE stays at full family width.
    series : the declared CONFIGURATION SERIES
        (:class:`atomterm.basis.SeriesSpec` tuple) -- the prior axis's
        radial-signal sources: the atom's valence configurations as
        forced-occupation, basin-gated references, the HEAD (first
        entry) being the case's own reference configuration.  Required
        by any ``@s``-tagged window (the series source); the mfs
        resolve once per atom through :func:`atomterm.basis.run_series`
        (the caller's job -- reference construction never enters a
        cell run).
    prior_bands : the PRIOR axis's realization (chemical shell -> its
        widened band tuple) -- what ``@r``/``@s`` windows distill
        from, beside the case's incumbent ``shell_bands`` (None =
        the incumbent realization).  The k>1 tags need the width;
        a widening beyond the measured certificate is a DECLARED
        widening (the declaration cites it).
    planned_windows : further window labels of the campaign, declared
        but not yet run (Ti: the (4,10)/(4,14) partitions) -- the cell
        grid derives their full factor grids (:func:`atomterm.cells.planned`).
    planned_sectors : per-window sector overrides for planned windows
        (window label -> (na, nb)) -- a promoted-pair window's core
        count follows its sector (:func:`atomterm.system.ncore_of`:
        the (5,1) key freezes 4 of the 5 core pairs).  Absent = the
        case's default sector.
    lattice_windows : the lattice spans (empty = none) -- CASCI cells
        over growing m-complete shell spans at ``lattice_sector``
        (:func:`atomterm.cells.lattice_keys`).
    lattice_sector : the lattice series' sector ``(na, nb)`` -- WHAT
        the lattice studies: the free-ion series (the ns2 pair folded
        into the core) or the full valence.  Required when
        ``lattice_windows`` is non-empty.
    scan : the declared frozen-scan anatomy
        (:class:`~atomterm.scans.ScanSpec` -- family pair, angle grid,
        fixed extras; None = no scan; the runner is
        ``scripts/radial_scan.py``).
    tag : the workdir/checkpoint label (defaults to the symbol).
    ref : the :class:`~atomterm.reference.ReferenceSpec` (character-forced
        occupations + basin gate) where default SCF filling fails; None
        lets plain aufbau run.
    sa_weights : state-average weights for the primary sector (uniform
        over the ground term's M_L set by convention).
    run_casscf : the single-root stage (off on multi-solution atoms).
    nroots_primary : PySCF roots for the primary CASCI.
    term : the certified ground term ``(L, parity)`` the pinned lane
        targets (parity = the configuration's total inversion parity);
        None = not certified for lane cells.
    lane_windows : the pinned-lane legs -- window labels carrying
        their own core (``"3d3p4d@semi"``, :func:`atomterm.windows.core_shells`
        the tag grammar); each becomes one CI-level lane cell
        (:func:`atomterm.cells.lane_keys`).  The tag marks the lane
        cell: a tagged window NEVER collides with a free-route label.
    scenarios : this atom's :class:`Scenario` declarations (merged
        over :data:`DEFAULT_SCENARIOS` by :func:`scenarios_for`).
    window_names : window label -> chemical display name for the
        renderers (an uncontracted basis counts radial primitives,
        not shells: the stored addresses keep the construction
        labels; the table speaks shells).
    core_chemistry : effective core count -> the core's shell set
        (the noble-guess fallback lies when the ECP leaves
        non-noble closed shells: U14e's core is 6s6p7s, not Ne).
    anchor_note : the ROHF anchor row's filling note, in shells.
    render_primary : the primary window's family renders (False =
        the promoted windows carry the results; the primary rows
        stay in the store).
    render_spans : the review-stage family filter (None = every
        family renders): only the listed spans' families do -- the
        growth path as review proceeds; the store and the runs are
        untouched by it.
    """

    symbol: str
    basis: str | dict
    nelec_act: int
    ncore: int
    shell_patterns: list[str]
    shell_bands: dict[str, tuple[str, ...]] | None = None
    # the WINDOW-side compression depth (chemical shell -> k reference
    # radials x the m-set): a multi-band realization's chemically
    # sized window -- CAS(3,6)-shaped, not CAS(3,43)-shaped.  The CORE
    # ignores it (full family spans; the compression is a window-side
    # choice only); every shell of every window the case declares must
    # carry one (a mixed distilled/raw window is not offered).
    shell_radials: dict[str, int] | None = None
    # the prior axis's declarations: the configuration series (the
    # @s windows' radial-signal sources) and the widened realization
    # the @r/@s windows distill from
    series: tuple | None = None
    prior_bands: dict[str, tuple[str, ...]] | None = None
    ecp: str | None = None
    spin: int = 0
    tag: str | None = None
    key: str | None = None
    ref: ReferenceSpec | None = None
    sa_weights: list[float] | None = None
    run_casscf: bool = True
    nroots_primary: int = 9
    planned_windows: tuple[str, ...] = ()
    planned_sectors: dict[str, tuple[int, int]] | None = None
    lattice_windows: tuple[str, ...] = ()
    lattice_sector: tuple[int, int] | None = None
    scan: ScanSpec | None = None
    term: tuple[int, str] | None = None
    lane_windows: tuple[str, ...] = ()
    scenarios: tuple[Scenario, ...] = ()
    # the chemistry display layer (the renderers' vocabulary; the
    # stored addresses keep their construction labels):
    # window_names: window label -> chemical display name (the
    #   uncontracted bases count radial primitives, not shells).
    # core_chemistry: effective core count -> its shell set (the
    #   noble-guess fallback lies when the ECP leaves non-noble
    #   closed shells).
    # anchor_note: the ROHF row's filling note, in shells.
    # render_primary: the primary window's family renders (False =
    #   the promoted windows carry the results; the primary rows stay
    #   in the store).
    window_names: dict[str, str] | None = None
    core_chemistry: dict[int, str] | None = None
    anchor_note: str | None = None
    render_primary: bool = True
    render_spans: tuple[str, ...] | None = None

    def __post_init__(self):
        """Validate the declaration's own grammar: the term target and
        every lane window (its tag resolves, its core stays out of its
        window).  Fail here, at the seam that owns the grammar."""
        object.__setattr__(self, "key", self.key or self.symbol)
        object.__setattr__(self, "planned_sectors",
                           self.planned_sectors or {})
        if self.lane_windows and self.term is None:
            raise ValueError(
                "lane windows need the term target (term=(L, 'g' | 'u'))")
        if self.lattice_windows and self.lattice_sector is None:
            raise ValueError(
                "lattice windows need the lattice sector "
                "(lattice_sector=(na, nb) -- the series' own Ms)")
        na, nb = self.lattice_sector or (1, 0)
        if not (na >= nb >= 0 and na + nb > 0):
            raise ValueError(
                f"lattice_sector {self.lattice_sector!r} must be "
                "(na >= nb >= 0, n > 0)")
        if self.term is not None:
            L, parity = self.term
            if not (isinstance(L, int) and L >= 0) or parity not in ("g", "u"):
                raise ValueError(
                    f"term {self.term!r} must be (L >= 0, 'g' | 'u')")
        open_d = next((p for p in self.shell_patterns
                       if p.endswith("d")), None)
        for label in self.lane_windows:
            span, sep, tag = label.partition("@")
            shells = re.findall(r"\d+[spdfgh]", span)
            if not sep or not tag or "".join(shells) != span:
                raise ValueError(
                    f"lane window {label!r} must be <shells>@<core> "
                    "(e.g. '3d3p@semi')")
            core = core_shells(tag, open_d=open_d)
            overlap = sorted(set(shells) & set(core))
            if overlap:
                raise ValueError(
                    f"lane window {label!r}: {overlap} sit in both the "
                    "window and its core")
        if self.shell_radials is not None:
            # the tagged planned windows carry their k in the source
            # tag -- only the untagged ones read shell_radials here
            windows = ["".join(self.shell_patterns), *self.lattice_windows,
                       *(w.partition("@")[0]
                         for w in self.lane_windows),
                       *(w for w in self.planned_windows if "@" not in w)]
            shells = {s for w in windows
                      for s in re.findall(r"\d+[spdfgh]", w)}
            bare = sorted(shells - set(self.shell_radials))
            if bare:
                raise ValueError(
                    f"shell_radials is undeclared for {bare} -- every "
                    "window shell the case names carries its k (a mixed "
                    "distilled/raw window is not offered)")
        # the prior axis: every @-tagged declared window is a source
        # tag (core tags live on lane windows, which never enter this
        # validation), and the series source needs the series
        tagged = [w for w in (*self.planned_windows,
                              "".join(self.shell_patterns))
                  if "@" in w]
        stray = [w for w in tagged if prior_tag_of(w) is None]
        if stray:
            raise ValueError(
                f"the prior windows {stray} must carry a source tag "
                "<shells>@<r|s><k> (e.g. '3d4s@s2'); core tags belong "
                "to lane windows")
        series = [w for w in tagged
                  if (t := prior_tag_of(w)) is not None and t[0] == "s"]
        if series and self.series is None:
            raise ValueError(
                f"the series windows {series} need the case's declared "
                "series (SeriesSpec tuple -- basis.run_series's "
                "sources)")
        if self.series is not None and self.series[0].spin != self.spin:
            raise ValueError(
                f"the series head {self.series[0].tag!r} must be the "
                f"case's reference configuration (spin {self.spin}, "
                f"not {self.series[0].spin})")
        if self.render_spans is not None:
            declared = {"".join(self.shell_patterns)}
            declared |= {w.partition("@")[0] for w in self.planned_windows}
            declared |= set(self.lattice_windows or ())
            declared |= {w.partition("@")[0]
                         for w in (self.lane_windows or ())}
            stray = [s for s in self.render_spans if s not in declared]
            if stray:
                raise ValueError(
                    f"render_spans names undeclared spans {stray} -- a "
                    "typo would silently render a partial table")

    def atom_config(self) -> AtomConfig:
        """The spec as an :class:`AtomConfig` (``tag`` defaults to the
        symbol)."""
        return AtomConfig(symbol=self.symbol, basis=self.basis, ecp=self.ecp,
                          spin=self.spin, tag=self.tag or self.symbol)


# The declared campaign (one spec per atom) -- THE intent manifest: the
# cell grid (cells.window_grid/grid_keys/planned), the regen runner, and
# the table renderers all derive from these; nothing else re-declares
# atoms, windows, or the gated references (the one-declaration rule).
# New atom = new row here.

# The ECP atom's raw specs live in data/ (raw basis/ECP text carries no
# path-safe name -- hence the explicit tag). Its window tiers are the
# record's exponent-band lattices expressed as their exact shell
# subsets: this basis is its own exponent grid, so a band IS a set of
# whole shells -- m-complete like any shell window.
_DATA = Path(__file__).resolve().parents[2] / "data"
_U_BASIS = gto.parse((_DATA / "u_dyall_cut_basis.txt").read_text())
_U_ECP = (_DATA / "u_stuttgart_rlc_ecp78.txt").read_text()
_U_TIER = {
    "3f": ("9f", "10f", "11f", "13d"),
    "4f": ("8f", "9f", "10f", "11f", "13d"),
    "2f2d": ("9f", "10f", "12d", "13d"),
    "3f2d": ("9f", "10f", "11f", "12d", "13d"),
    "4f2d": ("8f", "9f", "10f", "11f", "12d", "13d"),
}

# the standing protocol for planned windows, enforced: their OPTIMIZER
# cells stay void until their CASCI anchors have been read (a full
# regen materializes on demand).  The prefixes name every label form
# the grid gives the windows (plain and NO-prefix).
_ANCHOR_FIRST = Scenario(
    name="anchor-first",
    cells=("ti_3d4d_3-1_casscf_", "ti_3d4dno2_3-1_casscf_",
           "ti_3d4d_3-1_sa-casscf_", "ti_3d4dno2_3-1_sa-casscf_",
           "ti_3d4s4p4d_3-1_casscf_", "ti_3d4s4p4dno4_3-1_casscf_",
           "ti_3d4s4p4d_3-1_sa-casscf_", "ti_3d4s4p4dno4_3-1_sa-casscf_"),
    reason="planned windows anchor first: no optimizer cell of a window "
           "whose CASCI anchors are unread (2026-11-02)")

# the prior-axis windows (the @r/@s tags), shared by the contracted
# transition metals -- the probe's prior set as declared cells
_PRIOR_WINDOWS = ("3d@r1", "3d@s1", "3d@s2",
                  "3d4s@r1", "3d4s@s1", "3d4s@s2")

# Ti's same-l planned windows (the 3d+4d partitions) wait on
# window_distilled's declared future seam (the member-preserving
# cross-shell orthogonalization): two same-l shells' distilled blocks
# overlap and the assembly refuses -- void until the seam lands
# the twins' @s windows: the series' non-head configurations need
# fresh forced-occupation SCFs, and every fresh run diverges on the
# uncontracted basis (measured 2026-11-11, map Issues) -- the cells
# hold their STRUCTURAL place (pending rows) but never run until
# that repair lands
def _unc_s_block(key: str) -> Scenario:
    """The twins' @s windows' guard: the series' non-head
    configurations need fresh forced-occupation SCFs, and every fresh
    run diverges on the uncontracted basis (measured 2026-11-11, map
    Issues) -- the cells hold their STRUCTURAL place (pending rows)
    but never run until that repair lands.  The guard is the LABEL
    prefix alone (sector- and method-free): the claim is
    window-level, and a sector-bearing prefix silently leaks every
    cell the planned_sectors mirror moves (caught on the twin port,
    2026-11-11)."""
    return Scenario(
        name="series-source-diverges",
        cells=tuple(f"{key}_{lbl}_"
                    for lbl in ("3d@s1", "3d@s1no1", "3d@s2", "3d@s2no1",
                                "3d4s@s1", "3d4s@s1no2",
                                "3d4s@s2", "3d4s@s2no2")),
        reason="the series source (the non-head configurations' fresh "
               "forced-occupation SCFs) diverges on this basis -- the "
               "fresh-SCF disease, measured 2026-11-11 (map Issues); "
               "the @s windows wait on that repair")

def _same_l_windows(windows) -> tuple[str, ...]:
    """The windows whose shell composition REPEATS an l (two d shells
    in one span, e.g. 3d+4d): their distilled member blocks overlap
    and ``window_distilled`` refuses to assemble (the S-orthonormality
    fail-fast).  Derived from the labels -- the class is shell
    composition, never atom identity."""
    import re
    out = []
    for w in windows:
        shells = re.findall(r"\d+[spdf]", w.partition("@")[0])
        ls = [sh[-1] for sh in shells]
        if len(ls) > len(set(ls)):
            out.append(w)
    return tuple(out)


def _same_l_void(key: str, windows) -> Scenario:
    """The same-l windows' guard, sector- and method-free (the claim
    is window-level: the distilled assembly cannot be built,
    whatever cell asks): they wait on the member-preserving
    cross-shell orthogonalization -- window_distilled's declared
    future seam."""
    return Scenario(
        name="same-l-window",
        cells=tuple(f"{key}_{lbl}_"
                    for w in windows
                    for lbl in (w, no_label(w))),
        reason="same-l distilled windows wait on the member-preserving "
               "cross-shell orthogonalization (window_distilled's "
               "declared future seam; void since the distillation wave, "
               "2026-11-10)")





def _prior_anchor_first(atom: str, sectors) -> Scenario:
    """The prior windows' anchor-first gate (every declared sector --
    the case default and any planned_sectors overrides): their
    OPTIMIZER cells stay void until the CASCI anchors run -- the
    prior comparison is CASCI; the treatment rows over a chosen
    basis come after the adoption decision (the radial-CASSCF legs
    re-home there)."""
    cells = []
    for w in _PRIOR_WINDOWS:
        for lbl in (w, no_label(w)):
            cells += [f"{atom.lower()}_{lbl}_{na}-{nb}_{m}_"
                      for na, nb in sectors
                      for m in ("casscf", "sa-casscf")]
    return Scenario(
        name="prior-anchor-first", cells=tuple(cells),
        reason="the prior-axis comparison is CASCI: optimizer cells of "
               "a @r/@s window stay void until the basis decision "
               "(2026-11-10)")

# U14e's gate (Ti's pattern, widened to the whole optimizer class):
# the f-shell optimizer class is the declared multi-solution class
# (``run_casscf=False``'s spirit made explicit) -- the CASCI anchors
# (the primary's four cores, the four tiers' four) run; every
# casscf/sa-casscf cell of every window stays void at the key, one
# ``--full`` away.  The third prefix per window catches the NO-label
# form (``<w>no<k>`` -- start='no' never occurs on a casci cell, so
# the prefix cannot touch an anchor).

# the promoted-pair windows (the restricted optimizer batch,
# user-directed 2026-11-05; shell identities from the band probe --
# logs/probe_bands.log: the valence f/d character sits on 10f/13d,
# the core 6s on 14s+15s, the 6p trio on 14p+15p, the 7s on 18s+19s):
#   10f13d   CAS(4,12)  the record's minimal f+d window
#   10f13d7s CAS(6,13)  the record's six13 -- the 7s pair promoted,
#                       the core shrinks to 4 (the full core minus the
#                       7s pair); it carries its own sector
#   10f13d14p CAS(4,15) the p-entry probe -- the physical 6p trio's
#                       dominant band in-window, its pair folded (the
#                       Ti 3d4s4p4d pattern); the m-complete 6p
#                       PROMOTION would need a character-selected core
#                       (new machinery -- the recorded future step, the
#                       record's own open question)
_U14E_PROMOTED = ("10f13d", "10f13d7s", "10f13d14p")
_U14E_PROMOTED_SECTORS = {"10f13d7s": (5, 1)}

_U14E_WINDOWS = ("".join(_U_TIER["3f2d"]),) + tuple(
    "".join(_U_TIER[t]) for t in ("3f", "4f", "2f2d", "4f2d")) \
    + _U14E_PROMOTED

# THE SHELL LAYER (the port of the multi-column core, the machinery
# proven on Sc-unc/Ti-unc): the measured realizations from the record's
# band probe (elements/U14e/logs/probe_bands.log -- the closed 6s on
# 13s-15s, the 6p trio on 13p-15p, the 7s on 18s-20s; the open 5f on
# 9f-11f, the 6d on 12d-14d).  The BAND-labeled windows above stay
# their own identities (an undeclared pattern passes through the
# layer unchanged); the LANE window is chemical -- 5f6d@6s6p7s, the
# core 5 pairs over the 15-column family span, its frozen core the
# span's own closed-shell Roothaan.
_U14E_BANDS = {
    "5f": ("9f", "10f", "11f"),
    "6d": ("12d", "13d", "14d"),
    "6s": ("13s", "14s", "15s"),
    "6p": ("13p", "14p", "15p"),
    "7s": ("18s", "19s", "20s"),
}


def _u14e_sector(window: str) -> tuple[int, int]:
    """A U14e window's sector: the promoted-pair override, else the
    case's primary sector."""
    return _U14E_PROMOTED_SECTORS.get(window, (4, 0))


_U14E_ANCHOR_FIRST = Scenario(
    name="anchor-first",
    cells=tuple(f"u14e_{w}_{na}-{nb}_{m}_"
                for w in _U14E_WINDOWS
                for na, nb in (_u14e_sector(w),)
                for m in ("casscf", "sa-casscf"))
    + tuple(f"u14e_{w}no{k}_{na}-{nb}_{m}_"
            for w in _U14E_WINDOWS
            for na, nb in (_u14e_sector(w),)
            for k in (2, 3, 4, 5)
            for m in ("casscf", "sa-casscf")),
    reason="the optimizer class is multi-solution on the f shell "
           "(run_casscf=False's spirit made explicit): the CASCI "
           "anchors run first (2026-11-05; the kind block's casci "
           "cells -- mo/moc/no -- run: amended 2026-11-09)")
# the frac cores' gate: the fractional-occupation reference has no
# converged basin on the ECP atom (plain DIIS diverges E ~ +137;
# Newton-from-minao converges to a collapsed -12.74, 37 Eh above the
# aufbau RHF -49.04; Newton-from-ROHF does not leave the seed) -- the
# MOM-anchored block selection (the runner docstring's prescribed fix)
# is the recorded future step.
_U14E_FRAC_VOID = Scenario(
    name="frac-reference-unconverged",
    cells=tuple(f"u14e_{w}_{na}-{nb}_casci_ao_-_pol_frac"
                for w in _U14E_WINDOWS
                for na, nb in (_u14e_sector(w),)),
    reason="the fractional-occupation reference has no converged basin "
           "on the ECP atom: DIIS diverges, Newton collapses or does "
           "not move; MOM-anchored selection is the recorded fix "
           "(2026-11-05)")

# The uncontracted validation twins (the shell-layer wave): the same
# chemical shells and term targets as the certified rows, on the
# primitive split of the same basis (data/<el>_ccpvdz_unc.txt, the
# make_uncontracted.py split).  The SHELL LAYER is declared as
# CONTIGUOUS BAND RANGES partitioning each l-family (disjoint across
# shells; the boundaries from the measured composition,
# elements/<Key>/logs/band_composition.log -- the greedy capture
# tables locate each shell's MOs; the RANGES are the generous
# declaration).  The core spans are family-covering because the
# confined core is the span's OWN closed-shell Roothaan
# (core.confined_core_orbitals): its pairs fill the lowest orbitals
# of the span, so the core's realization must cover its shells'
# radial space DOWN TO THE CUSP (a norm-capture minimal list misses
# the tight primitives -- measured: 99.99% capture with +13 Eh
# core-energy error) and must EXCLUDE every skipped or active shell's
# bands (else the Roothaan fills them).  The references are
# trusted-cache chks: Sc-unc's own plain aufbau lands the deep basin;
# Ti-unc's plain aufbau lands a SHALLOW same-filling basin 11 mEh high
# and NO fresh start reaches the deepest one (the forced hooks
# diverge under every policy from every start) -- its rungs/rohf.chk
# is the CENSUSED deepest basin (logs/census.log: the seed-refinement
# DIIS attractor), trusted as every cached chk is.
_TM_UNC_BANDS = {
    "1s": tuple(f"{k}s" for k in range(1, 9)),
    "2s": tuple(f"{k}s" for k in range(9, 15)),
    "3s": tuple(f"{k}s" for k in range(15, 18)),
    "4s": ("18s", "19s", "20s"),
    "2p": tuple(f"{k}p" for k in range(2, 12)),
    "3p": tuple(f"{k}p" for k in range(12, 18)),
    "3d": tuple(f"{k}d" for k in range(3, 11)),
}
# the twins' WINDOW distillation depth (k=1: the single reference
# radial per shell, replicated on the m-set): the CAS labels mirror
# the certified contracted records' shapes (CAS(3,6)/CAS(4,6)) at
# minutes per cell.  The CORE families stay full -- the cusp is the
# core's business; the distillation is a window-side choice only.
_TM_UNC_RADIALS = {"3d": 1, "4s": 1}
_SC_UNC_BASIS = gto.parse((_DATA / "sc_ccpvdz_unc.txt").read_text())
_TI_UNC_BASIS = gto.parse((_DATA / "ti_ccpvdz_unc.txt").read_text())

# the validation-wave gates: the optimizer class stays void on the
# uncontracted bases (not declared this wave -- the CASCI anchors and
# the lane legs are the deliverable), and the frac cores wait for a
# measured fractional-occupation basin (the same conservatism the ECP
# atom's record set).
def _unc_gates(key: str, sec: tuple[int, int]) -> tuple[Scenario, ...]:
    w = "3d4s"
    return (
        Scenario(
            name="anchor-first",
            cells=tuple(f"{key}_{win}_{na}-{nb}_{m}_"
                        for win in (w, w + "no2", w + "@ar")
                        for m in ("casscf", "sa-casscf")
                        for na, nb in (sec,)),
            reason="the optimizer class is not declared on the "
            "uncontracted bases this wave: the CASCI anchors and the "
            "multi-column lane legs run first (2026-11-06; the kind "
            "block's casci cells -- mo/moc/no -- run: amended "
            "2026-11-09)"),
        Scenario(
            name="frac-reference-unmeasured",
            cells=tuple(f"{key}_{w}_{na}-{nb}_casci_ao_-_pol_frac"
                        for na, nb in (sec,)),
            reason="no measured fractional-occupation basin on the "
            "uncontracted basis yet (2026-11-06)"),
    )


CAMPAIGN: dict[str, CaseSpec] = {
    "C": CaseSpec("C", "ccpvdz", nelec_act=2, ncore=2,
                  shell_patterns=["2p"], spin=2, sa_weights=[1 / 3] * 3,
                  term=(1, "g"), lane_windows=("2p@1s2s",),
                  render_spans=("2p",),
                  # the prior axis (the p-shell port of the Sc
                  # grammar): the head = the trusted anchor; the
                  # promotion member 2s1 2p3 (quintet S, censused
                  # 2026-11-12: three arms, one basin, all
                  # occupation-safe) -- the widened chain 2p -> 3p
                  series=(SeriesSpec("p2", ("1s", "2s"), ("2p",),
                                     2, -37.68241788145875),
                          SeriesSpec("s1p3", ("1s",), ("2s", "2p"),
                                     4, -37.59068260)),
                  prior_bands={"2p": ("2p", "3p")},
                  planned_windows=("2p@r1", "2p@s1", "2p@s2"),
                  scenarios=(_prior_anchor_first("C", [(2, 0)]),),
                  scan=ScanSpec(fam=("2p", "3p"),
                                grid=(0.0, 0.8, 17))),
    "Si": CaseSpec("Si", "ccpvdz", nelec_act=2, ncore=6,
                   shell_patterns=["3p"], spin=2, sa_weights=[1 / 3] * 3,
                   term=(1, "g"), lane_windows=("3p@1s2s2p3s",),
                   render_spans=("3p",),
                   # the prior axis, mirrored (the 3p -> 4p chain; the
                   # promotion 3s1 3p3 censused 2026-11-12: one basin,
                   # all arms occupation-safe)
                   series=(SeriesSpec("p2", ("1s", "2s", "2p", "3s"),
                                      ("3p",), 2, -288.8464369230717),
                           SeriesSpec("s1p3", ("1s", "2s", "2p"),
                                      ("3s", "3p"), 4, -288.75397878)),
                   prior_bands={"3p": ("3p", "4p")},
                   planned_windows=("3p@r1", "3p@s1", "3p@s2"),
                   scenarios=(_prior_anchor_first("Si", [(2, 0)]),),
                   scan=ScanSpec(fam=("3p", "4p"),
                                 grid=(0.0, 1.1, 23))),
    "Sc": CaseSpec("Sc", "ccpvdz", nelec_act=3, ncore=9,
                   shell_patterns=["3d", "4s"], spin=1,
                   shell_radials={"3d": 1, "4s": 1, "3p": 1, "3s": 1,
                                  "4d": 1, "4f": 1, "4p": 1, "5d": 1},
                   sa_weights=[1 / 5] * 5,
                   term=(2, "g"),
                   lane_windows=("3d4s@ar", "3d@ar",
                                 "3d3p@semi", "3d3p4d@semi",
                                 "3d3p4s@1s2s2p3s",
                                 "3s3p3d@1s2s2p4s"),
                   lattice_windows=("3d", "3d4s", "3d4d", "3d4s4p",
                                    "3d4s4d", "3d4s4p4d",
                                    "3d4s4p4d5d", "3d4s4p4d4f"),
                   lattice_sector=(2, 1),
                   scan=ScanSpec(fam=("3d", "4d"), grid=(0.0, 0.4, 9),
                                 # 4s fixed: the 4s->5s direction is
                                 # harmful and would mask the d structure
                                 extra_fixed="4s"),
                   # the prior axis: the series (gates = the CENSUS
                   # occupation-safe basins, 2026-11-10 -- the probe's
                   # draws were 155-177 mEh SHALLOW garbage) and the
                   # widened realization (Sc's certificate: 3d ->
                   # (3d, 4d), 4s -> (4s, 6s))
                   series=(SeriesSpec("d1s2", ("1s", "2s", "2p", "3s",
                                               "3p", "4s"), ("3d",),
                                      1, -759.73625510),
                           SeriesSpec("d2s1", ("1s", "2s", "2p", "3s",
                                               "3p"), ("3d", "4s"),
                                      3, -759.69899090),
                           SeriesSpec("d3", ("1s", "2s", "2p", "3s",
                                             "3p"), ("3d",),
                                      3, -759.57064777)),
                   prior_bands={"3d": ("3d", "4d"), "4s": ("4s", "6s")},
                   # the prior windows: @r1 (the head's own block -- the
                   # probe's "neutral"), @s1/@s2 (the series mean, the
                   # mean+spread); @r2 does not exist on these atoms (a
                   # single-configuration block is rank-1 -- series_window
                   # refuses the noise second radial)
                   planned_windows=_PRIOR_WINDOWS,
                   # the review stage: the 3d4s window (CAS(3,6) +
                   # the @s2 CAS(3,12)) and the 3d section (CAS(3,5))
                   # -- opened with the === section rule
                   render_spans=("3d4s", "3d",),
                   scenarios=(_prior_anchor_first("Sc", [(2, 1)]),),
                   ref=ReferenceSpec(
                       occ_closed=("1s", "2s", "2p", "3s", "3p", "4s"),
                       occ_open=("3d",), occ_e_ref=-759.73625510,
                       occ_e_ref_uhf=-759.73789123)),
    "Ti": CaseSpec("Ti", "ccpvdz", nelec_act=4, ncore=9,
                   shell_patterns=["3d", "4s"], spin=2,
                   shell_radials={"3d": 1, "4s": 1, "4p": 1, "4d": 1,
                                  "3p": 1, "3s": 1},
                   sa_weights=[1 / 7] * 7,
                   term=(3, "g"),
                   lane_windows=("3d4s@ar", "3d@ar", "3d4d@ar",
                                 "3d4s4p4d@ar",
                                 "3d3p@semi", "3d3p4d@semi",
                                 "3s3p3d@1s2s2p4s",
                                 "3d3p4s@1s2s2p3s"),
                   render_spans=("3d4s", "3d"),
                   scenarios=(_ANCHOR_FIRST,
                              _prior_anchor_first("Ti", [(3, 1), (2, 0)])),
                   # the prior axis: the series (the HEAD gate = the
                   # TRUSTED reference energy; the further gates = the
                   # CENSUS occupation-safe basins 2026-11-10 -- the
                   # probe's d4 draw was 94 mEh shallow, d3s1 4.7 mEh);
                   # the realization WIDENS 3d beyond Ti's certificate
                   # (3d -> (3d), capture 0.9991 -- a DECLARED widening:
                   # k=2 needs the 2-wide radial space)
                   series=(SeriesSpec("d2s2", ("1s", "2s", "2p", "3s",
                                               "3p", "4s"), ("3d",),
                                      2, -848.40648875),
                           SeriesSpec("d3s1", ("1s", "2s", "2p", "3s",
                                               "3p"), ("3d", "4s"),
                                      4, -848.38650790),
                           SeriesSpec("d4", ("1s", "2s", "2p", "3s",
                                             "3p"), ("3d",),
                                      4, -848.24887522)),
                   prior_bands={"3d": ("3d", "4d"), "4s": ("4s", "6s")},
                   planned_windows=("3d4d", "3d4s4p4d", *_PRIOR_WINDOWS),
                   # the 3d-only prior windows compare against the
                   # LATTICE's d2 sector (ti_3d_2-0), not the case
                   # default (3,1 -- a different term class with no
                   # incumbent twin)
                   planned_sectors={f"3d@{t}": (2, 0)
                                    for t in ("r1", "s1", "s2")},
                   ref=ReferenceSpec(
                       occ_closed=("1s", "2s", "3s", "4s", "2p", "3p"),
                       occ_open=("3d",), occ_e_ref=-848.40648875,
                       occ_e_ref_uhf=-848.40991227),
                   lattice_windows=("3d", "3d4p", "3d4p4d"),
                   lattice_sector=(2, 0)),
    # The uncontracted twins: the multi-column core's validation atoms
    # (the ar core = 9 pairs over a 27-column declared span; the 3d4s
    # window COMPRESSED to its reference radials -- CAS(3,6)/CAS(4,6),
    # the certified records' shapes).  CASCI-class only this wave
    # (run_casscf=False: the optimizer class is not declared on the
    # uncontracted bases yet).
    "Sc-unc": CaseSpec("Sc", _SC_UNC_BASIS, nelec_act=3, ncore=9,
                        shell_patterns=["3d", "4s"], spin=1,
                        shell_bands=_TM_UNC_BANDS,
                        shell_radials=_TM_UNC_RADIALS,
                        # the prior axis, @r1 ONLY: the head's own block
                        # (the head IS the trusted cache -- no fresh SCF;
                        # the @s windows need the series' other
                        # configurations, whose fresh forced-occupation
                        # runs diverge on this basis -- measured
                        # 2026-11-11, the twins' fresh-SCF disease; the
                        # port's @s half waits on that repair)
                        series=(SeriesSpec("d1s2", ("1s", "2s", "2p", "3s",
                                                    "3p", "4s"), ("3d",),
                                           1, -759.7373848622801),),
                        prior_bands={"3d": _TM_UNC_BANDS["3d"],
                                     "4s": _TM_UNC_BANDS["4s"]},
                        planned_windows=("3d@r1", "3d@s1", "3d@s2",
                                         "3d4s@r1", "3d4s@s1",
                                         "3d4s@s2"),
                        sa_weights=[1 / 5] * 5,
                        run_casscf=False, key="Sc-unc",
                        render_spans=("3d4s", "3d"),
                        scenarios=_unc_gates("sc-unc", (2, 1))
                        + (_prior_anchor_first("sc-unc", [(2, 1)]),
                           _unc_s_block("sc-unc")),
                        term=(2, "g"),
                        lane_windows=("3d4s@ar", "3d@ar"),
                        lattice_windows=("3d",), lattice_sector=(2, 1)),
    "Ti-unc": CaseSpec("Ti", _TI_UNC_BASIS, nelec_act=4, ncore=9,
                        shell_patterns=["3d", "4s"], spin=2,
                        shell_bands=_TM_UNC_BANDS,
                        shell_radials=_TM_UNC_RADIALS,
                        # the prior axis, @r1-only (the head IS the
                        # trusted cache; the @s windows guarded -- the
                        # fresh-SCF disease, map Issues)
                        series=(SeriesSpec("d2s2", ("1s", "2s", "2p", "3s",
                                                    "3p", "4s"), ("3d",),
                                           2, -848.4074434287603),),
                        prior_bands={"3d": _TM_UNC_BANDS["3d"],
                                     "4s": _TM_UNC_BANDS["4s"]},
                        planned_windows=("3d@r1", "3d@s1", "3d@s2",
                                         "3d4s@r1", "3d4s@s1",
                                         "3d4s@s2"),
                        # the prior sectors mirror Ti's: the 3d-span
                        # priors at the lattice's d2 sector
                        planned_sectors={f"3d@{t}": (2, 0)
                                         for t in ("r1", "s1", "s2")},
                        sa_weights=[1 / 7] * 7,
                        run_casscf=False, key="Ti-unc",
                        render_spans=("3d4s", "3d"),
                        scenarios=_unc_gates("ti-unc", (3, 1))
                        + (_prior_anchor_first("ti-unc", [(3, 1), (2, 0)]),
                           _unc_s_block("ti-unc")),
                        term=(3, "g"),
                        lane_windows=("3d4s@ar", "3d@ar"),
                        lattice_windows=("3d",), lattice_sector=(2, 0)),
    # CAMPAIGN key `U14e`: the ECP-atom record's element dir
    # (``elements/U14e/``); the SYMBOL stays `U` (PySCF builds the
    # element; the ECP removes its 78 core electrons) and the `tag`
    # stays `U_stuttgart78_dyall` (the basis/psp label on chk names).
    "U14e": CaseSpec("U", _U_BASIS, nelec_act=4, ncore=5,
                   shell_patterns=list(_U_TIER["3f2d"]), spin=4,
                   ecp=_U_ECP, tag="U_stuttgart78_dyall", key="U14e",
                   shell_bands=_U14E_BANDS,
                   # every window shell carries k=1 -- the IDENTITY
                   # for the single-primitive band shells (a 1x1
                   # radial Gram; the stored band-window data stays
                   # valid) and the actual DISTILLATION for the 7s
                   # (its realization is the 3-band split -- one
                   # radial delivers the chemical width, CAS(6,13)).
                   # Corollary: the case is distilled -- the radial
                   # mask is construction-refused on its windows
                   shell_radials={"5f": 1, "6d": 1, "7s": 1,
                                  "9f": 1, "10f": 1, "11f": 1,
                                  "12d": 1, "13d": 1, "14p": 1},
                   sa_weights=[1 / 17] * 17,  # the ground term's M_L set
                   run_casscf=False,          # the multi-solution class
                   # the lane class RETURNS at the CHEMICAL width
                   # (the case now distills -- the lane window rides
                   # _window_columns like every other, delivering the
                   # CAS(4,12)/(6,13) shapes instead of the retired
                   # 36-column band span): one lane per section at its
                   # ruler core.  The (4,15) section carries NONE: its
                   # window shares the 6p rows with the 5-pair core
                   # and _lane_core refuses the overlap (a declared
                   # non-disjoint partition is a grid axis, not a
                   # lane one)
                   lane_windows=("5f6d@6s6p7s", "5f6d7s@6s6p"),
                   nroots_primary=48,
                   # the tier CAS(4,>20) families: RETIRED for now
                   # (user-directed 2026-11-06 -- "most of CASCI
                   # results (4,>20) can be deleted for now": the
                   # band-lattice anchor rows left the table with the
                   # declaration; the windows stay in git history and
                   # the store, re-declarable by re-adding the tuple)
                   # + tuple("".join(_U_TIER[t])
                   #         for t in ("3f", "4f", "2f2d", "4f2d")),
                   planned_windows=_U14E_PROMOTED,
                   planned_sectors=_U14E_PROMOTED_SECTORS,
                   render_spans=_U14E_PROMOTED + ("5f6d", "5f6d7s"),
                   scenarios=(_U14E_ANCHOR_FIRST, _U14E_FRAC_VOID),
                   window_names={"9f10f11f12d13d": "5f6d (3f2d radials)",
                                 "10f13d": "5f6d",
                                 "10f13d7s": "5f6d7s",
                                 "10f13d14p": "5f6d6p"},
                   core_chemistry={5: "Rn·7s", 4: "Rn"},
                   anchor_note="forced filling: 6s²6p⁶7s² closed, "
                               "5f³6d¹ open",
                   render_primary=False,
                   term=(8, "u"),
                   ref=ReferenceSpec(
                       occ_closed=("l0", "l1"), occ_open=("l2", "l3"),
                       occ_e_ref=-50.38185233,
                       occ_e_ref_uhf=-50.38881343,
                       conv_accelerator="newton")),
}


def scenarios_for(spec: CaseSpec) -> tuple[Scenario, ...]:
    """The spec's effective scenarios: the campaign defaults
    (:data:`DEFAULT_SCENARIOS`) merged with the row's own -- the row
    wins by name.  The ONE evaluated policy object the run gate
    reads; disabled scenarios stay in the tuple (their ``enabled``
    flag is the consumer's check).
    """
    merged = {**DEFAULT_SCENARIOS, **{s.name: s for s in spec.scenarios}}
    # the same-l guard DERIVES from the declared windows (ladder +
    # planned): the class is shell composition, never a per-atom
    # bookkeeping item -- a row that declares its own by name wins
    same_l = _same_l_windows((*spec.lattice_windows, *spec.planned_windows))
    if same_l and "same-l-window" not in merged:
        merged["same-l-window"] = _same_l_void(spec.key.lower(), same_l)
    return tuple(merged.values())


def campaign_spec(atom: str) -> CaseSpec:
    """The declared spec by campaign key -- the ONE lookup, owning the
    ``key == registry key`` contract: every element dir, store, and
    chkfile path derives from the key, so a spec whose ``key`` field
    disagrees with its ``CAMPAIGN`` entry refuses here instead of
    splitting the artifacts across two dirs.

    Raises
    ------
    KeyError
        the atom is not declared.
    ValueError
        the spec's ``key`` field disagrees with its registry key.
    """
    spec = CAMPAIGN[atom]
    if spec.key != atom:
        raise ValueError(
            f"CAMPAIGN key {atom!r} != spec.key {spec.key!r} -- the "
            "element dir would split between the two; fix the row")
    return spec


# the executor -------------------------------------------------------

def run_case(spec: CaseSpec, workdir: str = "elements") -> dict:
    """Execute a declared case: the staged pipeline's launch entry for
    any atom, new or old (an atom is runnable iff declared here).

    Parameters
    ----------
    spec : the :class:`CaseSpec` to execute.
    workdir : where the checkpoints and reports land (generated inserts
        go to ``elements/<Atom>/``).

    Returns
    -------
    the :func:`run_atom` stage dict (``rohf_energy`` + one entry per
    computed stage tag).
    """
    res = run_atom(
        spec.atom_config(),
        list(spec.shell_patterns),
        spec.nelec_act,
        workdir=workdir,
        sa_weights=list(spec.sa_weights) if spec.sa_weights else None,
        run_casscf=spec.run_casscf,
        nroots_primary=spec.nroots_primary,
        ref=spec.ref,
    )
    return res


def run_atom(
    cfg: AtomConfig,
    shell_patterns: list[str],
    nelec_act: int,
    workdir: str = "elements",
    sa_weights: list[float] | None = None,
    run_casscf: bool = True,
    nroots_primary: int = 9,
    ref: ReferenceSpec | None = None,
) -> dict:
    """End-to-end pipeline for any single atom.

    Parameters
    ----------
    cfg : the :class:`AtomConfig` of the atom.
    shell_patterns : valence-shell label prefixes selecting the active AO
        space, concatenated in order (e.g. ``["2p"]`` for C 2p^2,
        ``["3d", "4s"]`` for Sc 3d^1 4s^2). Each pattern must match at least
        one AO; each should contribute a complete m-set (rotation closure).
    nelec_act : electrons treated as active. The Ms sectors follow from
        cfg.spin = 2S: the primary sector (na, nb) = ((n+2S)/2, (n-2S)/2),
        plus an Ms=0 sector (n/2, n/2) when n is even and S > 0 (skipped when
        odd -- no Ms=0 exists, e.g. a doublet).
    sa_weights : run SA-CASSCF over these state weights (primary sector).
    run_casscf : run the single-root CASSCF stage (off for transition metals,
        where orbital optimization is multi-solution).
    nroots_primary : PySCF roots for the primary CASCI (clamped to the sector
        size; the Ms=0 CASCI always diagonalizes the full sector so the
        complete multiplet table is produced).
    ref : passed through to :func:`run_rohf` as a :class:`ReferenceSpec`
        -- the character-forced occupation patterns + reference-basin
        energy gate, mandatory where the default ROHF is an aufbau
        runaway (Ti; see :func:`run_rohf`).
    workdir : where checkpoints and reports land.

    Raises
    ------
    ValueError
        ``spin`` negative or of unequal parity with ``nelec_act`` (no
        such Ms sector exists).

    Every stage is checkpointed under ``workdir`` and analyzed with the exact
    L/S operator matrices; see :func:`run_carbon` for the carbon instance.

    Returns
    -------
    dict
        ``"rohf_energy"`` plus one payload per stage tag
        (``casci_ms*``, ``casscf``, ``sa-casscf``), each carrying the
        analysis object, the run's energies, window weights, max|dE|,
        and the particle-number sanity.
    """
    out = Path(workdir)
    out.mkdir(parents=True, exist_ok=True)
    label = cfg.label
    spin = cfg.spin  # 2S (PySCF convention)
    if spin < 0 or (nelec_act + spin) % 2 != 0:
        raise ValueError(
            f"spin 2S={spin} and nelec_act={nelec_act} must be >= 0 and "
            "of equal parity")

    # ---- stage 1: ROHF (checkpointed) -------------------------------------
    rohf_chk = str(out / f"{label}_rohf.chk")
    mol, mf = run_rohf(cfg, rohf_chk, ref=ref)
    print(f"[ROHF      ] E = {mf.e_tot:.8f} Eh   (chk: {rohf_chk})")

    act_idx = select_valence_ao(mol, shell_patterns)
    ncas = len(act_idx)

    # Ms sectors: (na, nb, nroots) for each CASCI diagonalization
    na_p, nb_p = sector_na_nb(nelec_act, spin)
    sector_defs = [(na_p, nb_p, min(nroots_primary, sector_dimension(ncas, na_p, nb_p)))]
    if nelec_act % 2 == 0 and spin > 0:
        na0 = nelec_act // 2
        sector_defs.append((na0, na0, sector_dimension(ncas, na0, na0)))

    results: dict[str, dict] = {"rohf_energy": float(mf.e_tot)}
    reports: list[str] = []

    def stage_analysis(tag: str, cas, c_act, na: int, nb: int):
        """The shared per-stage render: core + analyze through
        :func:`stage_facts`, print the report block, accumulate the
        results dict (one closure per run -- it captures the run's
        mol/mf/reports)."""
        ncore = ncore_of(mol.nelectron, na + nb)
        # build with the stage's OWN core (inside stage_facts: it uses
        # cas.mo_coeff[:, :ncore]): orbital optimization may legally
        # reorganize the inactive span, and a different core than the one
        # PySCF's energies used breaks E_tot = E_core + e_cas. For CASCI
        # stages the passed core IS the mf core (built by build_full_mo).
        ms = (na - nb) / 2  # the true M_S (chkfile tags label M_S too)
        ms_txt = "0" if ms == 0 else f"{ms:+g}"
        sector = f"Ms={ms_txt} (na={na}, nb={nb})"
        facts = stage_facts(mol, mf, cas, c_act, ncore, na, nb, tag=sector)
        title = (f"[{tag}]  total energies vs PySCF: "
                 f"max|dE| = {facts.max_de_vs_pyscf:.2e} Eh")
        text = format_report(facts.analysis, title)
        off_shell = 1.0 - facts.window_rows.max(axis=1)
        text += (f"\nwindow composition (l-weights): "
                 + ", ".join(f"{n}: {v:.3f}" for n, v in
                             facts.window.items() if v > 1e-4)
                 + f"   (mean off-dominant-l weight {off_shell.mean():.2e};"
                   " l-purity is necessary, not sufficient, for closure)")
        print(text)
        reports.append(text)
        results[tag] = {
            "analysis": facts.analysis,
            "energies_pyscf": np.asarray(facts.energies_pyscf).tolist(),
            "window_weights": facts.window,
            "max_de_vs_pyscf": facts.max_de_vs_pyscf,
            # particle-number sanity of the built determinant space: N is
            # diagonal by construction and must equal na + nb everywhere
            "n_dev": facts.n_dev,
        }

    # ---- stage 2: CASCI in every Ms sector (checkpointed) ------------------
    mo_by_sector: dict[tuple[int, int], np.ndarray] = {}
    for na, nb, nroots in sector_defs:
        ncore = ncore_of(mol.nelectron, na + nb)
        mo = build_full_mo(mf, occupied_core_idx(mf, ncore), act_idx)
        mo_by_sector[(na, nb)] = mo
        tag = f"casci_ms{_ms_label(na, nb)}"
        cas = run_casci(mf, ncas, (na, nb), mo, str(out / f"{label}_{tag}.chk"),
                        nroots=nroots, verbose=0)
        print(f"[CASCI Ms={_ms_label(na, nb)}] roots: {np.round(cas.energies, 6)}")
        stage_analysis(tag, cas, cas.mo_coeff[:, ncore:ncore + ncas], na, nb)

    mo_p = mo_by_sector[(na_p, nb_p)]
    ncore_p = ncore_of(mol.nelectron, na_p + nb_p)

    # ---- stage 3: CASSCF (single root) on the primary sector ---------------
    if run_casscf:
        casscf_chk = str(out / f"{label}_casscf.chk")
        cass = constrained_casscf(mf, ncas, (na_p, nb_p), casscf_chk,
                                  mo_full=mo_p, name="casscf", verbose=0)
        print(f"[CASSCF    ] E = {cass.energies[0]:.8f} Eh")
        c_cass = cass.mo_coeff[:, ncore_p : ncore_p + ncas]
        stage_analysis("casscf", cass, c_cass, na_p, nb_p)

    # ---- stage 4: SA-CASSCF over the requested state weights ---------------
    if sa_weights:
        sa_chk = str(out / f"{label}_sa-casscf.chk")
        sa = constrained_casscf(mf, ncas, (na_p, nb_p), sa_chk,
                                mo_full=mo_p, weights=list(sa_weights),
                                name="sa-casscf", verbose=0)
        print(f"[SA-CASSCF ] roots: {np.round(sa.energies, 6)}")
        c_sa = sa.mo_coeff[:, ncore_p : ncore_p + ncas]
        stage_analysis("sa-casscf", sa, c_sa, na_p, nb_p)

    # ---- persist -----------------------------------------------------------
    (out / f"{label}_report.md").write_text(
        f"# Atomic term analysis: {label}\n\n"
        f"ROHF/cc-pVDZ energy: {mf.e_tot:.8f} Eh\n\n" + "\n\n".join(reports),
        encoding="utf-8",
    )
    serial = {
        tag: {
            "energies_pyscf": payload["energies_pyscf"],
            "states": [vars(s) for s in payload["analysis"].states],
            "lz_spectrum": payload["analysis"].lz_spectrum,
            "sz_spectrum": payload["analysis"].sz_spectrum,
            "comm_norm_h_lz": payload["analysis"].comm_norm_h_lz,
            "comm_norm_h_l2": payload["analysis"].comm_norm_h_l2,
            "window_weights": payload["window_weights"],
            "max_de_vs_pyscf": payload["max_de_vs_pyscf"],
        }
        for tag, payload in results.items()
        if tag != "rohf_energy"
    }
    (out / f"{label}_report.json").write_text(
        json.dumps({"rohf_energy": results["rohf_energy"], **serial}, indent=2),
        encoding="utf-8",
    )
    print(f"\nreport written to {out}/{label}_report.md and .json")
    return results


def run_carbon(cfg: AtomConfig | None = None, workdir: str = "elements") -> dict:
    """Carbon pilot (2p^2, CAS(2,3)): the generic pipeline driven by
    the DECLARED row -- the anatomy (window, spin, SA weights) reads
    off :data:`CAMPAIGN`\\ ['C'], never re-declared here."""
    row = CAMPAIGN["C"]
    return run_atom(
        cfg or row.atom_config(),
        row.shell_patterns,
        row.nelec_act,
        workdir=workdir,
        sa_weights=row.sa_weights,
    )


if __name__ == "__main__":  # pragma: no cover
    run_carbon()
