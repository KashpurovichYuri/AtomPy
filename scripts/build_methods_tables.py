"""Render the methods-table section of an element's TABLE.md from its
cells.json payloads -- the key-addressed provenance (the render
contract: the table and its source cannot diverge).

Usage: python scripts/build_methods_tables.py [atom]

THE ONE METHODS TABLE (the element file IS the table): every run cell
of the atom a row, ROHF anchor first, the UHF anchor beside it when
the store carries the record; families in declaration order,
separated by the rule row -- the lattice spans (CI only), the
primary window (each core group's CASCI BESIDE its optimizer
variants -- CI interleaved, the row-collapse policy is DELETED),
then every lane window whose route twins ran (the pinned and free
CI, the free CASSCF/SA, the pinned SA).  NO-start cells are not
methods rows.  The PLANNED windows follow as pointers when gated
(declared, not yet run -- one regen away).

No verdict column: the classes and their row-level discriminators
are the reading key in README -- the one home.

Notes
-----
Column conventions and the reading key live in README (Methods and
Results) -- the one home; this renderer formats, it does not
re-state them. Two operational rules it does own:

- the anchor line in each caption is COMPUTED from the store's
  ``rohf_energy`` (Eh + eV), never typed;
- rows branch on key IDENTITY only (method/start/mask/substrate/
  core) and on the ROW POLICY; order and grouping live once in the
  declaration builder (cells.window_grid / cells.route_keys) -- this
  renderer receives ready-made key lists and formats.

Applicability: any declared atom (``cases.CAMPAIGN``; windows and
cell addresses derive from the declaration -- this file owns
presentation only: headings, notes, value formatting).

methods-table edit protocol lives here: sections are rebuilt from constructed
blocks, structural assertions checked (each heading once, row counts,
no pending rows), temp file + verified move.
"""

import json
import re
import sys

import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atomterm.cases import campaign_spec  # noqa: E402
from atomterm.cells import (CellKey, band_tokens, grid_keys,  # noqa: E402
                            lattice_keys, prior_tag_of, radial_refused,
                            route_keys, window_grid)
from _tables_common import (DASH, EV, cas_label,  # noqa: E402
                            comm_cell, conv_cell, core_col,
                            l2_text, minus, payload_of, replace_block)
BOUNDARY = "## Studies (frozen history; ticket 28)"

# the methods-table grammar (one shape with the CAS table); the
# family separator row renders as a horizontal rule mid-table
HEAD = ("| method | route | window | CAS | core | dets | E, Eh | "
        "ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | "
        "‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |")
SEP = "|" + "---|" * 18
RULE = "|" + " --- |" * 18
# the section rule: a family OUTSIDE the primary span opens a new
# SECTION (the added spans), set apart from the primary's families
HEAVY = "|" + " === |" * 18
BEGIN, END = "<!-- methods:begin -->", "<!-- methods:end -->"

# identity -> display (the ONLY branching the renderer does)
METHOD_NAMES = {"casci": "CASCI", "casscf": "CASSCF",
                "sa-casscf": "SA-CASSCF"}
# the distilled kind renders D (its construction: the distilled windows)
START_NAMES = {"ao": "D", "mo": "MO", "moc": "MOC", "no": "NO"}
MASK_NAMES = {None: "free", "radial": "radial", "frozen": "frozen",
              "virtmix": "virtual-mix", "frozencore": "frozen core"}


CSV_HEADER = ["method", "route", "window", "CAS", "core", "dets", "E_Eh",
              "dE_mEh", "E_eV", "dE_eV", "conv", "L2", "ptp_Eh", "PLz",
              "PLperp", "HLz", "HL2", "note"]


def csv_of(body: str) -> str:
    """The CSV mirror of the methods table: the same data rows (the
    rule rows and the markdown furniture dropped), ASCII header
    names, csv-quoted (the braced core names carry commas)."""
    import csv
    import io
    out = io.StringIO()
    w = csv.writer(out, lineterminator="\n")
    w.writerow(CSV_HEADER)
    for line in body.splitlines():
        if not line.startswith("|"):
            continue
        if set(line.replace("|", "").replace(" ", "")) == {"-"}:
            continue                      # the SEP and RULE rows
        cells = [c.strip() for c in line.split("|")[1:-1]]
        if cells[0] == "method":
            continue                      # the markdown header row
        w.writerow(cells)
    return out.getvalue()


def _method(k) -> str:
    """The method column: the objective + its restriction only (the
    start rides in the window cell, the substrate in the core
    column; an unrestricted optimizer names nothing beyond its
    method)."""
    name = METHOD_NAMES[k.method]
    if k.method != "casci" and k.mask is not None:
        name += f", {MASK_NAMES[k.mask]}"
    return name


def _e0(p: dict) -> float:
    """A payload's as-run comparison energy (the collapse ruler):
    the lane payloads' TERM total or the free cells' run total."""
    return float(p["total"]) if "total" in p \
        else float(p["energies_pyscf"][0])


def _num(v, fmt) -> str:
    """A numeric cell, or the em-dash when absent/NaN."""
    return format(v, fmt) if isinstance(v, float) and v == v else DASH


def _route_of(k):
    """The route column: the pinned lane class vs the free route."""
    return ("pinned" if k.core_tag is not None and k.substrate == "pur"
            else "free")


def _row(k, p, anchor, spec, sa_floor) -> str:
    """One run cell in the methods table's columns -- every run cell a row,
    fully identified (method, route, window/start, CAS, core)."""
    route = _route_of(k)
    core = core_col(k, spec)
    if k.core_tag is not None:
        core += "·lane"     # the raw-shell Loewdin core, not the prefix
    e = _e0(p)
    de = e - anchor
    ptp = p.get("ptp")
    if ptp is None and len(p.get("energies_pyscf", [])) > 1:
        ptp = float(np.ptp(p["energies_pyscf"]))   # the SA root set
    cyc, cap = p.get("cycles"), p.get("max_cycle")
    # the cap-hitting flag by IDENTITY on payloads carrying the cap
    # (the store's max_cycle field -- the run's actual mc cap);
    # payloads predating the field keep the legacy == 200 key (the
    # cap the measurement runs declared; PySCF's default is 50)
    note = "cycle cap" if cyc and cyc == (cap if cap is not None else 200) \
        else ""
    # an UNCONVERGED optimizer renders flagged, never clean (the row
    # stays -- the store is the record -- but the reading must know)
    if p.get("converged") is False:
        note = "not converged" + ("; " if note else "") + note
    # the swap-basin disclosure: the optimizer rotated the core away
    # from the reference core and the analysis read the raw field
    # (stages._CORE_ANCHOR_TOL) -- the row's energy is the swap
    # solution's own, never comparable to the declared window's
    if p.get("swap_basin"):
        note += ("; " if note else "") + "swap"
    if note:
        note += "; "
    # the kind vocabulary's verdict language: the SCF-MO kinds are
    # CONTROLS (m-incomplete by construction -- the SCF's hybrid-open
    # gauge; the scrambled term reading is the certified fingerprint,
    # an explained expectation, never a construction failure); the
    # pinned rows carry the bookkeeping-form disclosure (the residual
    # vs the anchor is the projected-density ruler's form, mEh-class)
    if k.start in ("mo", "moc"):
        note += "m-incomplete"
    elif route == "pinned":
        note += "bookkeeping form"
    # the single-root floor: a CASSCF row above its span's SA ground
    # root is a converged LOCAL minimum (1 mEh guards the cross-core
    # bookkeeping)
    if (k.method == "casscf" and k.mask is None and sa_floor is not None
            and (fl := sa_floor.get(
                (k.span_label, k.sector, p.get("n_det")))) is not None
            and e > fl + 1e-3):
        note += ("; " if note else "") + "local minimum"
    # a prior window renders under its full label (the tag is part of
    # the window's identity -- the address grammar verbatim)
    if k.source_tag is not None:
        wname = k.window
        note += ("; " if note else "") + f"bands {k.span_label}"
    else:
        wname = (spec.window_names or {}).get(k.span_label, k.span_label)
        # a CHEMICAL display name: the basis-band counters demote to
        # the note column (the provenance the user asked to keep)
        if wname != k.span_label:
            bands = "+".join(band_tokens(k.span_label))
            note += ("; " if note else "") + f"bands {bands}"
    note = note.rstrip("; ")    # no dangling separator on the last part
    cells = [
        _method(k), route, f"{wname} / {START_NAMES[k.start]}",
        cas_label(k, spec), core, f"{p['n_det']:,}",
        f"{e:.4f}", f"{de * 1e3:+.1f}", f"{e * EV:.4f}",
        f"{de * EV:+.4f}", conv_cell(p), l2_text(p),
        _num(ptp, ".0e"),
        _num(p.get("win_lz"), ".0e"),
        _num(p.get("win_lperp"), ".0e"),
        comm_cell(p, "comm_norm_h_lz", "res_h_lz"),
        comm_cell(p, "comm_norm_h_l2", "res_h_l2"),
        note,
    ]
    # sign minus on the ENERGY cells only (names keep their hyphens)
    for i in (6, 7, 8, 9):
        cells[i] = minus(cells[i])
    return "| " + " | ".join(cells) + " |"


def _method_row(k) -> bool:
    """The methods-row policy, ONE home: the table renders the
    COMPARISON classes -- the distilled block (the incumbent and the
    prior windows), the SCF controls, the treatments.  A NO-start
    cell never renders: its window is its SOURCE's span re-gauged
    into the term-averaged density's natural orbitals -- the CI is
    degenerate with the source's D row and the optimizer reaches the
    same basins (measured, Sc 2026-11-11: free to 5 decimals), a
    CERTIFICATE, not a comparison row; every filter in the table's
    assembly routes through here."""
    return k.start != "no"


def _sa_floor(cells: dict) -> dict[tuple, float]:
    """The single-root floor: the SA ground-root energy per
    (span, sector, determinant count).  A single-root CASSCF
    minimizes one state unconstrained by averaging, so it can never
    legitimately sit ABOVE its span's SA ground -- a row that does
    is a converged LOCAL minimum (the variational argument), and
    :func:`_row` flags it."""
    out: dict[tuple, float] = {}
    for addr, p in cells.items():
        if p.get("kind") != "sa-casscf" or "n_det" not in p:
            continue
        try:
            k = CellKey.from_address(addr)
        except (ValueError, IndexError):
            continue
        e = p.get("energies_pyscf", [p.get("e_ground")])[0]
        if e is None:
            continue
        key3 = (k.span_label, k.sector, p["n_det"])
        out[key3] = min(out.get(key3, float("inf")), float(e))
    return out


def _rows(keys, cells, anchor, spec, floor) -> list[str]:
    """Every declared cell of the key set (the :func:`_method_row`
    policy), in the given order: a payload renders its run, a
    payload-less runnable cell renders its PENDING placeholder (the
    structure), and the construction-refused class holds no place
    (:func:`_pendingable`)."""
    out = []
    for k in keys:
        if not _method_row(k):
            continue
        p = payload_of(cells, k)
        if p is not None and ("e_ground" in p or "energies_pyscf" in p):
            out.append(_row(k, p, anchor, spec, floor))
        elif _pendingable(k, spec):
            # DATA always renders (the store is the record -- legacy
            # rows included); the construction-refused class gates
            # only the PLACEHOLDER
            out.append(_pending_row(k, spec))
    return out


def _anchor_row(anchor: float, spec) -> str:
    """The ROHF anchor row: ΔE's zero (the caption states it); the
    filling note in the declared chemistry when the row declares it."""
    cells = ["ROHF (anchor)", "—", "—", "—", "—", "—",
             f"{anchor:.4f}", "0", f"{anchor * EV:.4f}", "+0.0000",
             DASH, DASH, DASH, DASH, DASH, DASH, DASH,
             spec.anchor_note or "forced filling"]
    for i in (6, 8):
        cells[i] = minus(cells[i])
    return "| " + " | ".join(cells) + " |"


def _uhf_anchor_row(anchor: float, uhf: dict) -> str:
    """The UHF anchor row: the spin-broken reference beside the ROHF
    one -- ΔE against the same zero, the meter reading (⟨S²⟩, the
    weakest corresponding-orbital pair) in the note."""
    e = uhf["e_tot"]
    de = (e - anchor) * 1000.0
    # a diverged record still renders (the measurement IS the verdict)
    # -- the note says so beside the meter reading
    note = ("diverged, " if uhf.get("converged") is False else "")
    cells = ["UHF (anchor)", "—", "—", "—", "—", "—",
             f"{e:.4f}", f"{de:+.1f}", f"{e * EV:.4f}",
             f"{de * EV / 1000.0:+.4f}",
             DASH, DASH, DASH, DASH, DASH, DASH, DASH,
             f"{note}⟨S²⟩={uhf['s2']:.4f}, "
             f"σ_min={uhf['sigma_min']:.4f}"]
    # sign minus on the ENERGY cells only (names keep their hyphens)
    for i in (6, 7, 8, 9):
        cells[i] = minus(cells[i])
    return "| " + " | ".join(cells) + " |"


def _uhf_pending_row() -> str:
    """The UHF anchor row's place when the store carries no record:
    dashed, noted pending (on the twins the class diverges -- map
    Issues -- the row states the missing measurement, never hides
    it)."""
    return "| " + " | ".join(
        ["UHF (anchor)", *[DASH] * 16, "pending"]) + " |"


def _pending_row(k, spec) -> str:
    """A declared, runnable cell that has no run yet: the row holds
    its place in the structure (label, CAS, core) with every value
    dashed and the note 'pending' -- the table doubles as the run
    plan; the runs land on approval."""
    wname = (k.window if k.source_tag is not None
             else (spec.window_names or {}).get(k.span_label, k.span_label))
    core = core_col(k, spec)
    if k.core_tag is not None:
        core += "·lane"     # the raw-shell Loewin core, not the prefix
    cells = [_method(k),
             _route_of(k),
             f"{wname} / {START_NAMES[k.start]}",
             cas_label(k, spec), core,
             *[DASH] * 12, "pending"]
    return "| " + " | ".join(cells) + " |"


def _pendingable(k, spec) -> bool:
    """The skeleton excludes the CONSTRUCTION-REFUSED class
    (:func:`atomterm.cells.radial_refused` -- the one home): those
    rows can never run, so they hold no place."""
    return not radial_refused(k, spec)


def _comparison_key(k):
    """The family's row order -- THE COMPARISON READING: every CI row
    before the optimizer rows; within CI, the distilled block by
    (k, source) -- the incumbent, the head prior ``r``, the series
    prior ``s`` at k=1, the double-budget ``s2`` -- the MO/MOC
    controls last; optimizer variants keep their declared order after
    the CI block."""
    t = k.source_tag
    k_rank, src = (t[1], t[0]) if t is not None else (1, "")
    # both tiers, ONE block structure: the D-family by (k, source),
    # then the MO/MOC controls last; within a start the OPTIMIZER block reads
    # CASSCF-then-SA, masks in the declaration order (free first)
    block = 0 if k.start == "ao" else 2
    tier = 0 if k.method == "casci" else 1
    m_rank = 1 if k.method == "sa-casscf" else 0
    mask_rank = {"radial": 1, "frozen": 2, "frozencore": 3,
                 "virtmix": 4}.get(k.mask, 0)
    return (tier, block, k_rank, src, k.start, m_rank, mask_rank)


def methods_section(spec, cells, anchor, uhf: dict | None = None) -> list[str]:
    """THE METHODS TABLE: one table for the whole atom.  Families in
    declaration order, separated by the rule row -- the lattice
    spans (the distilled CI row per span), the primary window when
    ``render_primary`` (each core group's CI block beside its
    optimizer variants), then every lane window whose route twins
    ran.  The PRIOR windows join their span's family (the same-CAS
    comparison: at a given CAS, every approach adjacent -- the
    incumbent's CI row, the ``r``/``s`` priors beside it, the
    double-budget block, then the optimizer rows over the chosen
    basis), in :func:`_comparison_key` order; a prior span with no
    incumbent family renders as its own family.
    A lattice span equal to the primary's renders once -- the
    primary's core groups own it.  The rendering DIET (store data
    stays; rows do not repeat it): the prior windows contribute
    their distilled rows in every core group where they ran (their
    controls repeat the incumbent's at equal width and stay out),
    the lattice renders one row per
    span (the m-incomplete control story renders once, at the
    primary's controls), the NO start never renders (its window is
    its source's span re-gauged -- a certificate, not a row), and
    the frac reference class does not render (a
    different reference SCF -- its home is the reference study,
    not the methods grid).  The ROHF anchor row
    leads, the UHF anchor row follows when the store carries the
    record; optimizer twins is not this table's story yet."""
    body: list[str] = []
    seen: set[str] = set()

    floor = _sa_floor(cells)
    sections: dict[str, list[list[str]]] = {}

    def _fam(keys) -> None:
        rows = _rows(sorted([k for k in keys
                             if k.address not in seen and _method_row(k)],
                            key=_comparison_key),
                     cells, anchor, spec, floor)
        seen.update(k.address for k in keys)
        if rows:
            # families collect under their CAS; the emission renders
            # each CAS as one section (=== between CAS, --- between
            # core families within) -- whatever the assembly's visit
            # order; the LANE families close their section (is_lane
            # sorts last)
            is_lane = any(k.core_tag is not None for k in keys)
            sections.setdefault(cas_label(keys[0], spec), []).append(
                (is_lane, rows))

    # the prior windows, gathered by span (the sector match guards
    # the join: a prior at another sector has no incumbent twin
    # there and renders as its own family).  The DIET filter: the
    # distilled CI row at the ruler core only -- the prior windows'
    # controls repeat the incumbent's at equal width and their
    # raw/pur rows repeat their certificates
    priors: dict[str, list] = {}
    for w in spec.planned_windows:
        if prior_tag_of(w) is None:
            continue
        priors.setdefault(w.partition("@")[0], []).extend(
            k for k in window_grid(spec, w) if k.start == "ao")
    primary_span = "".join(spec.shell_patterns) if spec.render_primary \
        else None

    # the review-stage family filter (None = every family): only
    # the listed spans' families render -- the growth path back to
    # the full table
    keep = (lambda sp: True) if spec.render_spans is None \
        else (lambda sp: sp in spec.render_spans)

    # the primary family FIRST (the main section): core groups in
    # declaration order, each group's CI block (the incumbent's row
    # beside the prior windows' -- the comparison) followed by its
    # optimizer variants; the frac reference class does not render
    # (the diet); suppressed by the declaration (render_primary=
    # False -- the promoted windows carry the results; the primary
    # rows stay in the store)
    if spec.render_primary and keep(primary_span):
        grid = [k for k in grid_keys(spec)
                if _method_row(k) and k.core != "frac"]
        join = [k for k in priors.get(primary_span, [])
                if k.sector == grid[0].sector]
        groups: dict[str, list] = {}
        for k in [*grid, *join]:
            groups.setdefault((core_col(k, spec), cas_label(k, spec)),
                              []).append(k)
        for keys in groups.values():
            _fam(keys)

    # the route/lane families of the primary span close its section
    routes = route_keys(spec)
    for leg in [k for k in routes
                if k.substrate == "pur" and k.method == "sa-casscf"]:
        if not keep(leg.window.partition("@")[0]):
            continue
        twins = [k for k in routes
                 if k.window == leg.window and k.core != "frac"]
        # the family renders whatever ran -- unrun twins hold their
        # pending places (the structural era; the any-twin-ran gate
        # retired with it)
        _fam([k for k in twins if _method_row(k)])

    # the other spans' sections follow (opened with the === rule)
    for span in dict.fromkeys(k.window for k in lattice_keys(spec)):
        if span == primary_span:
            continue    # the primary's core groups own this span
        if not keep(span):
            continue
        # the SAME core-group structure as the primary section: each
        # group's block (the incumbent's rows beside the prior
        # windows' at the lattice sector), pending placeholders
        # holding the unrun places
        keys = [k for k in lattice_keys(spec)
                if k.window == span and k.core != "frac"]
        keys += [k for k in priors.get(span, [])
                 if k.sector == spec.lattice_sector]
        groups: dict[tuple, list] = {}
        for k in keys:
            groups.setdefault((core_col(k, spec), cas_label(k, spec)),
                              []).append(k)
        for g in groups.values():
            _fam(g)

    # the planned windows that did not join a family above: the
    # untagged tiers (their own families, as declared) and any prior
    # span without an incumbent twin (its own comparison family) --
    # the SAME core-group grammar as the primary/lattice sections
    # (raw | avg | pur behind --- rules), the same row order
    for w in spec.planned_windows:
        if not keep(w.partition("@")[0]):
            continue
        if prior_tag_of(w) is not None:
            keys = priors.pop(w.partition("@")[0], [])
            if not keys:
                continue
        else:
            keys = [k for k in window_grid(spec, w)
                    if _method_row(k) and k.core != "frac"]
        groups: dict[tuple, list] = {}
        for k in keys:
            groups.setdefault((core_col(k, spec), cas_label(k, spec)),
                              []).append(k)
        for g in groups.values():
            _fam(g)

    # sections emission: === between spans, --- between families
    emitted: list[str] = []
    prev_span = None
    for span, families in sections.items():
        for _, rows in sorted(families, key=lambda f: f[0]):
            if emitted:
                emitted.append(RULE if span == prev_span else HEAVY)
            emitted.extend(rows)
            prev_span = span
    out = [BEGIN, HEAD, SEP, _anchor_row(anchor, spec)]
    out.append(_uhf_anchor_row(anchor, uhf) if uhf is not None
               else _uhf_pending_row())
    out.append(HEAVY)     # the anchor block is its own SECTION
    out.extend([*emitted, END])

    return out


def main(atom: str) -> None:
    """Rebuild ``<atom>``'s methods block (the one table) from
    cells.json under the methods-table edit protocol (the marker pair;
    structural checks; verified move) and write the CSV mirror."""
    src = Path(f"elements/{atom}/TABLE.md")
    if not src.exists():
        sys.exit(f"{src}: not a one-table atom (no TABLE.md) -- the "
                 "atom's pass migrates it")
    data = json.loads(Path(f"elements/{atom}/cells.json").read_text())
    cells, anchor = data["cells"], data["rohf_energy"]
    spec = campaign_spec(atom)

    body = "\n".join(methods_section(spec, cells, anchor,
                                    uhf=data.get("uhf")))

    # structural assertions (the table protocol): every pending row
    # is well-formed (the pending placeholders are VALID rows -- the
    # table doubles as the run plan), no starred gauge numbers,
    # every table row well-formed
    assert not re.search(r"\d\*", body)      # no starred gauge numbers
    for line in body.splitlines():
        if line.startswith("|"):
            assert line.count("|") == HEAD.count("|"), line

    replace_block(src, BEGIN, END, body)
    src.with_suffix(".csv").write_text(csv_of(body))
    n = sum(1 for line in body.splitlines()
            if line.startswith(("| CASCI", "| CASSCF", "| SA-CASSCF")))
    print(f"{src}: methods table rendered ({n} run rows) + TABLE.csv mirror")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Ti")
