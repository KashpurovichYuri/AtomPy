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

import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.cells import (grid_keys, lane_keys, route_keys,  # noqa: E402
                            window_grid)
from _tables_common import (DASH, EV, cas_label,  # noqa: E402
                            comm_cell, conv_cell, fold_col,
                            l2_text, minus)
BOUNDARY = "## Studies (frozen history; ticket 28)"
POINTER = ("Declared, not yet run -- the optimizer cells are void "
           "until their CASCI anchors are read (the anchor-first "
           "scenario). One regen away: `regen_cells.py {atom} "
           "--window {window}`.")

# the ladder grammar (one shape with the CAS ladder table); the
# family separator row renders as a horizontal rule mid-table
HEAD = ("| method | route | window | CAS | fold | dets | E, Eh | "
        "ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | "
        "‖[H,Lz]‖ | ‖[H,L²]‖ | note |")
SEP = "|" + "---|" * 16
RULE = "|" + " --- |" * 16

# identity -> display (the ONLY branching the renderer does)
METHOD_NAMES = {"casci": "CASCI", "casscf": "CASSCF",
                "sa-casscf": "SA-CASSCF"}
START_NAMES = {"ao": "AO", "mo": "MO", "no": "NO"}
MASK_NAMES = {None: "free", "radial": "radial", "frozen": "frozen",
              "virtmix": "virtual-mix"}


def _method(k) -> str:
    """The method column: the full treatment name (method, mask,
    start) -- every axis of the run's identity is on the row: route,
    window, CAS and fold are columns; the substrate IS the fold
    (Ar = the canonical core, Ar*pur the purified span)."""
    return (f"{METHOD_NAMES[k.method]}, {MASK_NAMES[k.mask]}, "
            f"{START_NAMES[k.start]}")


def _run(p: dict | None) -> dict | None:
    """The payload of a RUN cell (None for pending and void cells)."""
    if p is None or "void_reason" in p:
        return None
    return p


def _e0(p: dict) -> float:
    """A payload's as-run comparison energy (the collapse ruler):
    the lane payloads' TERM total or the free cells' run total."""
    return float(p["total"]) if "total" in p \
        else float(p["energies_pyscf"][0])


def _row(k, p, anchor, spec) -> str:
    """One optimizer run cell in the ladder columns -- every run cell
    a row, fully named (no collapse: the absorbed rows' conv and
    meters were being hidden)."""
    route = ("pinned" if k.fold_tag is not None and k.substrate == "pur"
             else "free")
    fold = fold_col(k, spec)
    if k.fold_tag is not None:
        fold += "·lane"     # the raw-shell Loewdin fold, not the prefix
    e = _e0(p)
    de = e - anchor
    ptp = p.get("ptp")
    if ptp is None and len(p.get("energies_pyscf", [])) > 1:
        ptp = float(np.ptp(p["energies_pyscf"]))   # the SA root set
    note = "cycle cap" if p.get("cycles") == 200 else ""
    cells = [
        _method(k), route, k.span_label, cas_label(k), fold,
        f"{p['n_det']:,}",
        f"{e:.4f}", f"{de * 1e3:+.1f}", f"{e * EV:.4f}",
        f"{de * EV:+.4f}", conv_cell(p), l2_text(p),
        f"{ptp:.0e}" if ptp is not None else DASH,
        comm_cell(p, "comm_norm_h_lz", "res_h_lz"),
        comm_cell(p, "comm_norm_h_l2", "res_h_l2"),
        note,
    ]
    # sign minus on the ENERGY cells only (names keep their hyphens)
    for i in (6, 7, 8, 9):
        cells[i] = minus(cells[i])
    return "| " + " | ".join(cells) + " |"


def _rows(keys, cells, anchor, spec) -> list[str]:
    """Every OPTIMIZER run cell of the key set, in declaration order
    (CASCI cells are the survey's and the CAS ladder's story)."""
    return [_row(k, p, anchor, spec)
            for k in keys
            if k.method != "casci"
            and (p := _run(cells.get(k.address))) is not None]


def _anchor_row(anchor: float) -> str:
    """The ROHF anchor row: ΔE's zero (the caption states it)."""
    cells = ["ROHF (anchor)", "—", "—", "—", "—", "—",
             f"{anchor:.4f}", "0", f"{anchor * EV:.4f}", "+0.0000",
             DASH, DASH, DASH, DASH, DASH, "forced filling"]
    for i in (6, 8):
        cells[i] = minus(cells[i])
    return "| " + " | ".join(cells) + " |"


def _caption(anchor: float) -> list[str]:
    return ["Optimizer rows, every CAS window · ΔE vs the ROHF "
            f"anchor ({anchor:.4f} Eh) · source `cells.json` "
            "(as-run)."]


def ladder_section(spec, cells, anchor) -> list[str]:
    """THE METHOD LADDER: one table for the whole atom -- the primary
    window's optimizer grid (the ROHF anchor row leading), then every
    lane window whose route twins ran, window families separated by
    the rule row.  Gated planned windows follow as pointer lines; a
    window with no optimizer run cells is not this table's story."""
    from atomterm.cases import scenarios_for
    from atomterm.cells import void_reason

    body = [_anchor_row(anchor), *_rows(grid_keys(spec), cells,
                                        anchor, spec)]
    routes = route_keys(spec)
    families = []
    for leg in [k for k in routes
                if k.substrate == "pur" and k.method == "sa-casscf"]:
        twins = [k for k in routes
                 if k.window == leg.window and k.method != "casci"]
        if any(_run(cells.get(k.address)) for k in twins):
            families.append(_rows(twins, cells, anchor, spec))
    for fam in families:
        body += [RULE, *fam]

    out = ["## Method ladder", "", *_caption(anchor), "",
           HEAD, SEP, *body]

    scens = scenarios_for(spec)
    for w in spec.planned_windows:
        keys = window_grid(spec, w)
        if not any(k.method != "casci" and _run(cells.get(k.address))
                   for k in keys):
            out += ["", POINTER.format(atom=spec.symbol, window=w)]
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
    first = re.search(r"^## Method ladder", text, re.M)
    if first is None:
        sys.exit(f"{src}: no '## Method ladder' section to replace")

    body = "\n".join(ladder_section(spec, cells, anchor))
    # the splice eats EVERYTHING from the first ladder heading to the
    # boundary -- historical section forms cannot survive a render
    new = text[:first.start()] + body + "\n\n" \
        + text[text.index(BOUNDARY):]
    assert new.count("\n## Method ladder\n") == 1

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
    n = sum(1 for line in body.splitlines()
            if line.startswith(("| CASSCF", "| SA-CASSCF")))
    print(f"{src}: method ladder rendered from cells.json "
          f"({n} optimizer rows, ladder shape)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Ti")
