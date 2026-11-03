"""Render the CASCI window-survey section of an element's TABLES.md
from its cells.json payloads (the build_*_tables family's survey
member; the render contract: the table and its source cannot diverge).

The block between the ``<!-- survey:begin/end -->`` markers is
GENERATED: one row per declared CI cell -- the lattice spans (AO/MO
starts), the primary window's CASCI fold variants (raw/frac/pur), and
the pinned ladder legs (their own folds) -- as-run; a cell without a
stored payload renders ``pending``.  The prose around the markers is
HAND physics and stays untouched.

Usage:
    python scripts/build_survey_tables.py [atoms]

Atoms (default: every declared atom with survey rows, CAMPAIGN order).
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.cells import CellKey, grid_keys, lane_keys, lattice_keys  # noqa: E402
from _tables_common import (EV, SHELL_COLS, cas_label, comm_text,  # noqa: E402
                            fold_label, l2_text, minus, payload_of,
                            replace_block)
from atomterm.windows import NOBLE_CORES  # noqa: E402

HEAD = ("| window | CAS | fold | start | route | dets | E, Eh | "
        "ΔE, mEh | E, eV | ΔE, eV | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ |")
SEP = "|" + "---|" * 13
BEGIN, END = "<!-- survey:begin -->", "<!-- survey:end -->"


def _noble(ncore: int) -> str:
    """The case prefix's noble name (the fold of every untagged cell):
    the noble whose closed electron count fills ``ncore`` orbitals."""
    for name, shells in NOBLE_CORES.items():
        if sum(SHELL_COLS[s[-1]] for s in shells) == ncore:
            return name.capitalize()
    return f"n={ncore}"


def survey_keys(spec) -> list[CellKey]:
    """The declared CI cells, ordered by window family (the lattice
    span order, then the ladder-only spans): per span -- the free
    starts, the fold variants, then the pinned leg.  Deduped by
    address (the primary window's pol-avg CI cell IS its lattice AO
    row)."""
    lat = lattice_keys(spec)
    seen = {k.address for k in lat}
    variants = [k for k in grid_keys(spec)
                if k.method == "casci" and k.start == "ao"
                and (k.fold in ("raw", "frac") or k.substrate == "pur")
                and k.address not in seen]
    legs = [k for k in lane_keys(spec) if k.address not in seen]
    from dataclasses import replace

    twins = [replace(k, substrate="pol", method="casci")
             for k in lane_keys(spec)]
    legs += [k for k in twins if k.address not in seen
             and k.address not in {x.address for x in legs}]
    by_span: dict[str, list[CellKey]] = {}
    for k in lat + variants + legs:
        by_span.setdefault(k.span_label, []).append(k)
    order = list(spec.lattice_windows) + [
        s for s in dict.fromkeys(k.span_label for k in legs)
        if s not in spec.lattice_windows]
    return [k for s in order for k in by_span.get(s, [])]


def _fold_col(k: CellKey, spec) -> str:
    """The fold column: the pinned legs' own tag; the free cells'
    prefix (the variants name their deviation)."""
    if k.fold_tag is not None:
        return fold_label(k.fold_tag)
    base = _noble(spec.ncore)
    if k.fold == "raw":
        return f"{base}·raw"
    if k.fold == "frac":
        return f"{base}·frac"
    if k.substrate == "pur":
        return f"{base}·pur"
    return base


def section(spec, cells: dict) -> str:
    """The generated block between the markers: one row per declared
    CI cell, as-run, in window-family order.  ``ΔE`` anchors on the
    pinned leg spanning the case's primary window (asserted present --
    the column's zero)."""
    keys = survey_keys(spec)
    primary = "".join(spec.shell_patterns)
    anchor = next((k for k in keys if k.fold_tag is not None
                   and k.substrate == "pur"
                   and k.span_label == primary), None)
    if anchor is None:
        sys.exit(f"{spec.symbol}: no pinned leg spans the primary window "
                 f"{primary!r} -- the ΔE column needs it")
    a0 = payload_of(cells, anchor)
    if a0 is None:
        sys.exit(f"{spec.symbol}: the anchor leg {anchor.address} has no "
                 "stored payload -- run the lane first "
                 "(regen_cells.py --group lane)")
    e0 = float(a0["total"])

    rows = []
    for k in keys:
        p = payload_of(cells, k)
        route = ("pinned" if k.substrate == "pur" and k.fold_tag
                 is not None else "free")
        head = [k.span_label, cas_label(k), _fold_col(k, spec),
                k.start.upper(), route]
        if p is None:
            rows.append("| " + " | ".join(head + ["pending", "", "", "",
                                                  "", "", "", ""]) + " |")
            continue
        e = float(p["total"] if route == "pinned"
                  else p["energies_pyscf"][0])
        de = (e - e0) * 1e3
        vals = [f"{p['n_det']:,}", f"{e:.4f}",
                "0" if k is anchor else f"{de:+.1f}",
                f"{e * EV:.4f}", f"{(e - e0) * EV:+.4f}", l2_text(p),
                comm_text(p["comm_norm_h_lz"]),
                comm_text(p["comm_norm_h_l2"])]
        rows.append("| " + " | ".join(
            head + [minus(v) for v in vals]) + " |")
    return "\n".join([BEGIN, HEAD, SEP, *rows, END])


def run(atom: str) -> None:
    """Rebuild ``elements/<Atom>/TABLES.md``'s survey block from its
    cells.json payloads under the TABLES edit protocol."""
    spec = CAMPAIGN[atom]
    keys = survey_keys(spec)
    if not keys:
        sys.exit(f"no declared CI cells for {atom}")
    src = Path(f"elements/{atom}/TABLES.md")
    cells = json.loads(
        Path(f"elements/{atom}/cells.json").read_text())["cells"]
    new_block = section(spec, cells)
    replace_block(src, BEGIN, END, new_block)

    # structural assertions (the TABLES protocol)
    body = [r for r in new_block.splitlines() if r.startswith("|")]
    assert len(body) - 2 == len(keys), body
    for r, k in zip(body[2:], keys):
        assert r.startswith(f"| {k.span_label} |"), r
    n = sum(1 for k in keys if payload_of(cells, k))
    print(f"{src}: CI survey rendered from cells.json "
          f"({len(keys)} rows, {n} run)")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("atoms", nargs="*",
                   help="declared atoms with survey rows (default: all)")
    ns = p.parse_args(argv)
    chosen = ns.atoms or list(CAMPAIGN)
    unknown = [a for a in chosen if a not in CAMPAIGN]
    if unknown:
        sys.exit(f"not declared: {', '.join(unknown)}")
    for a in chosen:
        if not survey_keys(CAMPAIGN[a]):
            print(f"{a}: no declared CI cells -- skipped")
            continue
        run(a)


if __name__ == "__main__":
    main()
