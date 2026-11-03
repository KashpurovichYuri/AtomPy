"""Render the method-ladder section of an element's TABLES.md from its
cells.json payloads -- the key-addressed provenance (the render
contract: the table and its source cannot diverge).

Usage: python scripts/build_methods_tables.py [Ti]

One table per RUN declared window, the SUMMARY SHAPE (one treatment
column composes the identity -- the C/Si/Sc convention), under the ROW
POLICY (the compression): a cell without a stored payload (pending,
or void under a scenario) leaves the table -- a window with no run
cells renders a pointer line instead of a pending grid (the full grid
is one regen away); within a (method, mask, substrate, fold) group
the AO row is the head and MO/NO runs agreeing with it within EPS
collapse into its treatment name -- disagreement stays as its own row
(the basin and construction classes live there).  No verdict column:
the classes and their row-level discriminators are the reading key in
README -- the one home.

Notes
-----
Column conventions and the reading key live in README (Methods and
Results) -- the one home; this renderer formats, it does not
re-state them. Two operational rules it does own:

- the anchor line in each caption is COMPUTED from the store's
  ``rohf_energy`` (Eh + eV), never typed;
- rows branch on key IDENTITY only (method/start/mask/substrate/
  fold) and on the ROW POLICY; order and grouping live once in the
  declaration builder (cells.window_grid) -- this renderer receives
  ready-made key lists and formats.

Applicability: the Ti-style factor-grid tables (WINDOWS below; the
window set and every cell address come from the declaration --
``cases.CAMPAIGN`` via ``cells.grid_keys``/``cells.planned``; this
file owns presentation only: headings, notes, value formatting). The
study inserts are frozen history; the per-atom summary tables (C/Si/
Sc) stay hand-curated until their manifests migrate.

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
from atomterm.cases import CAMPAIGN
from atomterm.cells import grid_keys, window_grid
from _tables_common import DASH, EV, comm_text, minus
EPS = 1e-4          # Eh -- agreeing starts collapse into the AO head

# the declared windows this file renders, with their presentation
# (headings). The window list MUST equal the declaration's (primary +
# planned) -- asserted in main; every cell address derives from the
# declaration, not from here. Captions stay minimal (reading
# conventions live in README; the anchor line is COMPUTED from the
# store's rohf_energy).
WINDOWS = [
    ("## Method ladder, CAS(4,6)/{3d,4s}", "3d4s"),
    ("### CAS(4,10)/{3d,4d} — declared, gated", "3d4d"),
    ("### CAS(4,14)/{3d,4s,4p,4d} — declared, gated", "3d4s4p4d"),
]
BOUNDARY = "## Studies (frozen history; ticket 28)"
POINTER = ("Declared, not yet run -- the optimizer cells are void "
           "until their CASCI anchors are read (the anchor-first "
           "scenario). One regen away: `regen_cells.py Ti --window "
           "{window}`.")

# the summary shape (the C/Si/Sc convention): ONE treatment column
# composes the identity (method, mask, substrate, fold, starts); the
# meters and the span composition carry the reading
HEAD = ("| treatment | E, Eh | ΔE, mEh | E, eV | ΔE, eV | ⟨L²⟩ | "
        "‖[H,Lz]‖ | ‖[H,L²]‖ | Span (Lowdin) |")
SEP = "|" + "---|" * 9

# identity -> display (the ONLY branching the renderer does)
METHOD_NAMES = {"casci": "CASCI", "casscf": "CASSCF",
                "sa-casscf": "SA-CASSCF"}
START_NAMES = {"ao": "", "mo": ", MO", "no": ", NO"}
MASK_NAMES = {None: "free", "radial": "radial", "frozen": "frozen",
              "virtmix": "virtual-mix"}
FOLD_NAMES = {"raw": "raw", "avg": "averaged", "frac": "fractional"}


def _treatment(k, collapsed=()) -> str:
    """The treatment column: the cell's identity composed into one
    name -- CASCI names its fold (pur implies averaged, the canonical
    pur-fold rule); an optimizer names mask and the pur substrate
    (pol implied); a non-AO start names itself; a collapsed head
    names the starts it absorbed."""
    name = METHOD_NAMES[k.method]
    if k.method == "casci":
        return (f"{name}, pur" if k.substrate == "pur"
                else f"{name}, {FOLD_NAMES[k.fold]}")
    name += f", {MASK_NAMES[k.mask]}"
    if k.substrate == "pur":
        name += ", pur"
    name += START_NAMES[k.start]
    if collapsed:
        name += f" (+{','.join(collapsed)})"
    return name


def _cfg(w: dict) -> str:
    """The converged window's SPAN composition (S-metric Löwdin
    weights, the payload ``window``) -- fractions of the window's
    columns per l, summing to 1 (NOT the electron configuration).
    Compact: s0 p.17 d.83 f4e-4."""
    out = []
    for l in ("s", "p", "d", "f"):
        v = float(w.get(l, 0.0))
        if v < 5e-5:
            out.append(f"{l}0")
        elif v < 0.001:
            out.append(f"{l}{v:.0e}")
        else:
            out.append(f"{l}{v:.2f}")
    return " ".join(out)


def _values(p: dict, anchor: float) -> list[str]:
    """The value cells of a run cell: the as-run E in Eh/eV, ΔE vs the
    anchor, the meters, the span composition -- sign minus
    normalized.  Kernel-only payloads (the dense-meter budget class)
    carry NaN meters -- they render the em-dash; the window
    composition falls back to empty (the payload records no span)."""
    asrun = float(p["energies_pyscf"][0])
    de_eh = asrun - anchor
    g = float(p.get("ground_l2", float("nan")))
    out = [f"{asrun:.4f}", f"{de_eh * 1e3:+.1f}", f"{asrun * EV:.4f}",
           f"{de_eh * EV:+.4f}",
           f"{g:.4f}" if g == g else DASH,          # NaN check
           comm_text(p.get("comm_norm_h_lz", float("nan"))),
           comm_text(p.get("comm_norm_h_l2", float("nan"))),
           _cfg(p.get("window", {}))]
    return [minus(v) for v in out]   # sign minus only


def _run(p: dict | None) -> dict | None:
    """The payload of a RUN cell (None for pending and void cells)."""
    if p is None or "void_reason" in p:
        return None
    return p


def _row(cells, k, anchor, collapsed=()):
    """One table row from a run cell: the composed treatment name and
    the value cells."""
    p = _run(cells.get(k.address))
    assert p is not None, "the policy renders run cells only"
    vals = _values(p, anchor)
    return f"| {_treatment(k, collapsed)} | " + " | ".join(vals) + " |"


def _groups(keys, cells) -> dict[tuple, dict[str, tuple]]:
    """The window's RUN cells grouped by identity (method, mask,
    substrate, fold) -> {start: (key, payload)} -- the policy's one
    home (rendering and the structural assertions share it)."""
    groups: dict[tuple, dict[str, tuple]] = {}
    for k in keys:
        p = _run(cells.get(k.address))
        if p is None:
            continue
        groups.setdefault(
            (k.method, k.mask, k.substrate, k.fold), {})[k.start] = (k, p)
    return groups


def _collapsed(groups) -> dict[tuple, tuple[str, ...]]:
    """Per group: the non-AO starts agreeing with the AO head within
    EPS (they collapse into its start cell). Groups without a RUN AO
    head collapse nothing."""
    out = {}
    for g, members in groups.items():
        head = members.get("ao")
        if head is None:
            continue
        e0 = float(head[1]["energies_pyscf"][0])
        out[g] = tuple(
            st for st in ("mo", "no") if st in members
            and abs(float(members[st][1]["energies_pyscf"][0]) - e0)
            <= EPS)
    return out


def render_window(keys, cells, anchor) -> list[str]:
    """Format one window's ready-made key list under the ROW POLICY:
    run cells only; within a (method, mask, substrate, fold) group the
    AO head absorbs MO/NO runs agreeing within EPS (named in its start
    cell); disagreeing runs stay. Declaration order preserved."""
    groups = _groups(keys, cells)
    collapsed = _collapsed(groups)

    rows = []
    for k in keys:
        p = _run(cells.get(k.address))
        if p is None:
            continue
        g = (k.method, k.mask, k.substrate, k.fold)
        if k.start != "ao" and k.start in collapsed.get(g, ()):
            continue
        rows.append(_row(cells, k, anchor,
                         collapsed.get(g, ()) if k.start == "ao" else ()))
    return rows


def _anchor_note(anchor: float) -> str:
    """The computed anchor line (sign-minus normalized)."""
    return minus(f"ROHF anchor: {anchor:.4f} Eh ({anchor * EV:.4f} eV).")


def window_section(heading, window, keys, cells, anchor) -> str:
    """One window's generated block: the compressed table when any
    cell has run, else the pointer line (declared, gated/pending --
    no pending grid in the committed table)."""
    if not any(k for k in keys if _run(cells.get(k.address))):
        return "\n".join([heading, "", _anchor_note(anchor), "",
                          POINTER.format(window=window)])
    rows = render_window(keys, cells, anchor)
    return "\n".join([
        heading, "",
        "Route free · partition the fold column · source "
        "`cells.json` (the committed draw).",
        "", _anchor_note(anchor), "", HEAD, SEP, *rows])


def main(atom: str) -> None:
    """Rebuild ``<atom>``'s method-ladder section from cells.json under
    the render contract (WINDOWS cover exactly the declared windows --
    asserted; the TABLES edit protocol's structural checks; verified
    move)."""
    src = Path(f"elements/{atom}/TABLES.md")
    data = json.loads(Path(f"elements/{atom}/cells.json").read_text())
    cells, anchor = data["cells"], data["rohf_energy"]
    spec = CAMPAIGN[atom]

    # the declaration's per-window grids (primary + planned). Render
    # contract: this file's WINDOWS cover EXACTLY the declared windows,
    # primary first -- presentation here, content from the declaration.
    primary = "".join(spec.shell_patterns)
    per_window = {primary: grid_keys(spec)}
    per_window.update({w: window_grid(spec, w)
                       for w in spec.planned_windows})
    assert [w for _, w in WINDOWS] == [primary, *spec.planned_windows], \
        "the renderer's WINDOWS and cases.CAMPAIGN diverged"

    sections = [window_section(heading, window, per_window[window],
                               cells, anchor)
                for heading, window in WINDOWS]

    text = src.read_text()
    b0 = WINDOWS[0][0]
    assert text.count(b0) == 1 and text.count(BOUNDARY) == 1
    new = text[:text.index(b0)] + "\n\n\n".join(sections) + "\n\n" \
        + text[text.index(BOUNDARY):]

    # structural assertions (the TABLES protocol): one heading each,
    # the declared windows covered, and NO pending rows anywhere
    for heading, window in WINDOWS:
        assert new.count(heading) == 1, heading
    span = new[new.index(b0):new.index(BOUNDARY)]
    assert "pending" not in span
    assert not re.search(r"\d\*", span)          # no starred gauge numbers
    for heading, window in WINDOWS:
        h = new.index(heading)
        e = (new.index("\n### ", h) + 1 if "\n### " in new[h:]
             else new.index(BOUNDARY))
        rows = [r for r in new[h:e].splitlines() if r.startswith("|")]
        groups = _groups(per_window[window], cells)
        expect = len(render_window(per_window[window], cells, anchor))
        if expect:
            assert len(rows) - 2 == expect, (heading, len(rows), expect)
        else:
            assert not rows, heading        # the pointer line, no table

    tmp = src.with_suffix(".md.tmp")
    tmp.write_text(new)
    assert tmp.read_text() == new
    tmp.replace(src)
    n = sum(1 for _, w in WINDOWS for k in per_window[w]
            if _run(cells.get(k.address)))
    print(f"{src}: method ladder rendered from cells.json "
          f"({n} run cells, compressed)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Ti")
