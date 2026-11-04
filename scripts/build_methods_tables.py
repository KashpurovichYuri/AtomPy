"""Render the method-ladder section of an element's TABLES.md from its
cells.json payloads -- the key-addressed provenance (the render
contract: the table and its source cannot diverge).

Usage: python scripts/build_methods_tables.py [atom]

THE LADDER SHAPE (one grammar with the CAS ladder): route | method |
window | CAS | fold | dets | E | ΔE | eV | conv | ⟨L²⟩ | ptp |
comms -- one section per CAS window family, every OPTIMIZER run cell
a row (the CI story lives in the survey / the CAS ladder):

- the PRIMARY window: its factor grid's optimizer cells under the ROW
  POLICY (the compression) -- within a (method, mask, substrate,
  fold) group the AO row is the head and MO/NO runs agreeing with it
  within EPS collapse into its treatment name -- disagreement stays
  as its own row (the basin and construction classes live there),
  and the ROHF anchor row leads (ΔE's zero);
- the PLANNED windows: their optimizer cells when they ran, else the
  pointer line (declared, gated -- one regen away);
- every LANE window whose route twins ran: the free CASSCF/SA and
  the pinned SA beside it (the methods axis across CAS sizes).

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
  fold) and on the ROW POLICY; order and grouping live once in the
  declaration builder (cells.window_grid / cells.route_keys) -- this
  renderer receives ready-made key lists and formats.

Applicability: any declared atom (``cases.CAMPAIGN``; windows and
cell addresses derive from the declaration -- this file owns
presentation only: headings, notes, value formatting).

TABLES edit protocol lives here: sections are rebuilt from constructed
blocks, structural assertions checked (each heading once, row counts,
no pending rows), temp file + verified move.
"""

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.cells import (grid_keys, lane_keys, route_keys,  # noqa: E402
                            window_grid)
from _tables_common import (DASH, EV, SHELL_COLS, cas_label,  # noqa: E402
                            comm_cell, conv_cell, fold_col, fold_label,
                            l2_text, minus)
EPS = 1e-4          # Eh -- agreeing starts collapse into the AO head

BOUNDARY = "## Studies (frozen history; ticket 28)"
POINTER = ("Declared, not yet run -- the optimizer cells are void "
           "until their CASCI anchors are read (the anchor-first "
           "scenario). One regen away: `regen_cells.py {atom} "
           "--window {window}`.")

# the ladder grammar (one shape with the CAS ladder table)
HEAD = ("| route | method | window | CAS | fold | dets | E, Eh | "
        "ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | "
        "‖[H,Lz]‖ | ‖[H,L²]‖ |")
SEP = "|" + "---|" * 15

# identity -> display (the ONLY branching the renderer does)
METHOD_NAMES = {"casci": "CASCI", "casscf": "CASSCF",
                "sa-casscf": "SA-CASSCF"}
START_NAMES = {"ao": "", "mo": ", MO", "no": ", NO"}
MASK_NAMES = {None: "free", "radial": "radial", "frozen": "frozen",
              "virtmix": "virtual-mix"}


def _treatment(k, collapsed=()) -> str:
    """The method column: the cell's composed treatment identity
    (mask, substrate, starts) -- the ladder's window/CAS/fold/route
    columns carry the rest; a collapsed head names the starts it
    absorbed."""
    name = METHOD_NAMES[k.method]
    if k.fold_tag is not None:
        return name        # tagged cells: the route column carries it
    name += f", {MASK_NAMES[k.mask]}"
    if k.substrate == "pur":
        name += ", pur"
    name += START_NAMES[k.start]
    if collapsed:
        name += f" (+{','.join(collapsed)})"
    return name


def _run(p: dict | None) -> dict | None:
    """The payload of a RUN cell (None for pending and void cells)."""
    if p is None or "void_reason" in p:
        return None
    return p


def _groups(keys, cells) -> dict[tuple, dict[str, tuple]]:
    """The RUN cells grouped by identity (method, mask, substrate,
    fold) -> {start: (key, payload)} -- the policy's one home
    (rendering and the structural assertions share it)."""
    groups: dict[tuple, dict[str, tuple]] = {}
    for k in keys:
        p = _run(cells.get(k.address))
        if p is None:
            continue
        groups.setdefault(
            (k.method, k.mask, k.substrate, k.fold), {})[k.start] = (k, p)
    return groups


def _e0(p: dict) -> float:
    """A payload's as-run comparison energy (the collapse ruler):
    the lane payloads' TERM total or the free cells' run total."""
    return float(p["total"]) if "total" in p \
        else float(p["energies_pyscf"][0])


def _collapsed(groups) -> dict[tuple, tuple[str, ...]]:
    """Per group: the non-AO starts agreeing with the AO head within
    EPS (they collapse into its start cell). Groups without a RUN AO
    head collapse nothing."""
    out = {}
    for g, members in groups.items():
        head = members.get("ao")
        if head is None:
            continue
        e0 = _e0(head[1])
        out[g] = tuple(
            st for st in ("mo", "no") if st in members
            and abs(_e0(members[st][1]) - e0) <= EPS)
    return out


def _row(k, p, anchor, spec, collapsed=()) -> str:
    """One optimizer run cell in the ladder columns."""
    route = ("pinned" if k.fold_tag is not None and k.substrate == "pur"
             else "free")
    e = _e0(p)
    de = e - anchor
    ptp = p.get("ptp")
    cells = [
        route, _treatment(k, collapsed), k.span_label, cas_label(k),
        fold_col(k, spec), f"{p['n_det']:,}",
        f"{e:.4f}", f"{de * 1e3:+.1f}", f"{e * EV:.4f}",
        f"{de * EV:+.4f}", conv_cell(p), l2_text(p),
        f"{ptp:.0e}" if ptp is not None else DASH,
        comm_cell(p, "comm_norm_h_lz", "res_h_lz"),
        comm_cell(p, "comm_norm_h_l2", "res_h_l2"),
    ]
    # sign minus on the ENERGY cells only (method names keep their
    # ASCII hyphens)
    for i in (6, 7, 8, 9):
        cells[i] = minus(cells[i])
    return "| " + " | ".join(cells) + " |"


def _grid_rows(keys, cells, anchor, spec) -> list[str]:
    """The primary window's optimizer cells under the row policy, in
    declaration order: the AO head per group (naming absorbed
    starts), disagreeing starts as their own rows."""
    groups = _groups(keys, cells)
    collapsed = _collapsed(groups)
    rows = []
    for k in keys:
        if k.method == "casci":
            continue
        p = _run(cells.get(k.address))
        if p is None:
            continue
        g = (k.method, k.mask, k.substrate, k.fold)
        members = groups[g]
        if k.start != "ao" and k.start in collapsed.get(g, ()):
            continue                    # absorbed into the AO head
        coll = collapsed.get(g, ()) if k.start == "ao" else ()
        rows.append(_row(k, p, anchor, spec, coll))
    return rows


def _anchor_row(anchor: float) -> str:
    """The ROHF anchor row: ΔE's zero (the caption states it)."""
    cells = ["—", "ROHF (anchor)", "—", "—", "—", "—",
             f"{anchor:.4f}", "0", f"{anchor * EV:.4f}", "+0.0000",
             DASH, DASH, DASH, DASH, DASH]
    return "| " + " | ".join(minus(v) for v in cells) + " |"


def _caption(anchor: float) -> list[str]:
    return ["Optimizer rows across CAS windows · ΔE vs the ROHF "
            f"anchor ({anchor:.4f} Eh) · routes in rows · source "
            "`cells.json` (as-run)."]


def _table(keys, cells, anchor, spec, lead_anchor: bool) -> list[str]:
    rows = ([_anchor_row(anchor)] if lead_anchor else []) \
        + _grid_rows(keys, cells, anchor, spec)
    return ["", HEAD, SEP, *rows] if rows else []


def _section(heading, keys, cells, anchor, spec,
             lead_anchor=False) -> list[str]:
    """One window family's block: the ladder table when an optimizer
    cell ran (a grid section may lead with the ROHF anchor row);
    pending/void grids render the pointer line -- no pending rows in
    the committed table."""
    has_run = any(k.method != "casci" and _run(cells.get(k.address))
                  for k in keys)
    if not has_run:
        return [heading, "", *_caption(anchor), "",
                POINTER.format(atom=keys[0].atom, window=keys[0].window)]
    return [heading, "", *_caption(anchor), *_table(keys, cells, anchor,
                                                    spec, lead_anchor)]


def _heading(prefix: str, k) -> str:
    """The section heading of a window family from any of its keys
    (the fold names the lane families' partition)."""
    shells = re.findall(r"\d+[spdfgh]", k.span_label)
    display = "{" + ",".join(shells) + "}"
    cas = f"({sum(k.sector)},{sum(SHELL_COLS[s[-1]] for s in shells)})"
    if k.fold_tag is not None:
        return f"{prefix} CAS{cas}/{display} · {fold_label(k.fold_tag)}"
    return f"{prefix} CAS{cas}/{display}"


def sections_of(spec, cells, anchor) -> list[list[str]]:
    """The method ladder's sections, in declaration order: the
    primary window's optimizer grid (leading the ROHF anchor row),
    the planned windows (their grids or the pointer), then every
    lane window whose route twins ran (the methods axis across CAS
    sizes).  Windows with no optimizer run cells render nothing --
    they are not this table's story yet."""
    from atomterm.cases import scenarios_for
    from atomterm.cells import void_reason

    primary = grid_keys(spec)
    heading = _heading("## Method ladder,", primary[0])
    out = [_section(heading, primary, cells, anchor, spec,
                    lead_anchor=True)]

    scens = scenarios_for(spec)
    for w in spec.planned_windows:
        keys = window_grid(spec, w)
        gated = any(void_reason(k.address, scens) for k in keys)
        mark = " — declared, gated" if gated else ""
        out.append(_section(_heading("###", keys[0]) + mark,
                            keys, cells, anchor, spec))
    routes = route_keys(spec)
    for leg in [k for k in routes
                if k.substrate == "pur" and k.method == "sa-casscf"]:
        twins = [k for k in routes
                 if k.window == leg.window and k.method != "casci"]
        ran = any(k.method != "casci" and _run(cells.get(k.address))
                  for k in twins)
        if ran:        # un-run lane windows are not this table's story
            out.append(_section(_heading("###", leg),
                                twins, cells, anchor, spec))
    return out


def main(atom: str) -> None:
    """Rebuild ``<atom>``'s method-ladder section from cells.json under
    the render contract (the declared windows derive from the
    CaseSpec; the TABLES edit protocol's structural checks; verified
    move)."""
    src = Path(f"elements/{atom}/TABLES.md")
    data = json.loads(Path(f"elements/{atom}/cells.json").read_text())
    cells, anchor = data["cells"], data["rohf_energy"]
    spec = CAMPAIGN[atom]

    text = src.read_text()
    if text.count(BOUNDARY) != 1:
        sys.exit(f"{src}: no {BOUNDARY!r} boundary -- place the "
                 "section family first (the atom's pass)")
    old = re.search(r"^## Method ladder.*?(?=^## )", text, re.M | re.S)
    if old is None:
        sys.exit(f"{src}: no '## Method ladder' section to replace")

    sections = sections_of(spec, cells, anchor)
    body = "\n\n\n".join("\n".join(s) for s in sections)
    new = text[:old.start()] + body + "\n\n" + text[old.end():]

    # structural assertions (the TABLES protocol): no pending rows,
    # no starred gauge numbers, every table row well-formed
    assert "pending" not in body
    assert not re.search(r"\d\*", body)      # no starred gauge numbers
    for line in body.splitlines():
        if line.startswith("|"):
            assert line.count("|") == HEAD.count("|"), line

    tmp = src.with_suffix(".md.tmp")
    tmp.write_text(new)
    assert tmp.read_text() == new
    tmp.replace(src)
    n = sum(1 for sec in sections for line in sec
            if line.startswith(("| free | CASSCF", "| free | SA-CASSCF",
                                "| pinned | CASSCF",
                                "| pinned | SA-CASSCF")))
    print(f"{src}: method ladder rendered from cells.json "
          f"({n} optimizer rows, ladder shape)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Ti")
