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
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.cells import CellKey, grid_keys, lane_keys, lattice_keys  # noqa: E402
from atomterm.windows import NOBLE_CORES, fold_shells  # noqa: E402

SHELL_COLS = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9, "h": 11}
EV = 27.2114
DASH = "—"
HEAD = ("| window | CAS | fold | start | route | dets | E, Eh | "
        "ΔE, mEh | E, eV | ΔE, eV | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ |")
SEP = "|" + "---|" * 13
BEGIN, END = "<!-- survey:begin -->", "<!-- survey:end -->"


def minus(text: str) -> str:
    """Sign-minus normalization (the tables' typographic convention)."""
    return re.sub(r"(?<!e)-", "−", text)


def _noble(ncore: int) -> str:
    """The case prefix's noble name (the fold of every untagged cell):
    the noble whose closed electron count fills ``ncore`` orbitals."""
    for name, shells in NOBLE_CORES.items():
        if sum(SHELL_COLS[s[-1]] for s in shells) == ncore:
            return name.capitalize()
    return f"n={ncore}"


def fold_label(tag: str) -> str:
    """The fold column's README notation: noble names capitalized,
    ``semi`` as itself, explicit lists braced."""
    if tag in NOBLE_CORES:
        return tag.capitalize()
    if tag == "semi":
        return tag
    return "{" + ",".join(fold_shells(tag)) + "}"


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
    by_span: dict[str, list[CellKey]] = {}
    for k in lat + variants + legs:
        by_span.setdefault(k.span_label, []).append(k)
    order = list(spec.lattice_windows) + [
        s for s in (k.span_label for k in legs)
        if s not in spec.lattice_windows]
    return [k for s in order for k in by_span.get(s, [])]


def _payload(cells: dict, k: CellKey) -> dict | None:
    """The run payload (None for pending and void cells)."""
    p = cells.get(k.address)
    if p is None or "void_reason" in p:
        return None
    return p


def _l2_text(p: dict) -> str:
    """⟨L²⟩: the free cells' ground root; the pinned legs' components
    -- one value when they agree, the min–max range when they
    scramble."""
    if "l2" in p:
        l2 = p["l2"]
        if max(l2) - min(l2) < 1e-6:
            return f"{sum(l2) / len(l2):.4f}"
        return f"{min(l2):.1f}–{max(l2):.1f}"
    return f"{p['ground_l2']:.4f}"


def _comm(v) -> str:
    """One commutator norm, or the em-dash when the sector sat above
    the dense-meter cap."""
    return f"{v:.1e}" if v == v else DASH          # NaN check: v != v


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


def _cas(k: CellKey) -> str:
    """The CAS(nelec, ncas) label: electrons from the key's sector,
    columns from the window span's summed m-set sizes."""
    cols = sum(SHELL_COLS[l] for l in
               re.findall(r"\d+([spdfgh])", k.span_label))
    return f"({sum(k.sector)},{cols})"


def section(spec, cells: dict) -> str:
    """The generated block between the markers: one row per declared
    CI cell, as-run, in window-family order.  ``ΔE`` anchors on the
    pinned leg spanning the case's primary window (asserted present --
    the column's zero)."""
    keys = survey_keys(spec)
    primary = "".join(spec.shell_patterns)
    anchor = next((k for k in keys if k.fold_tag is not None
                   and k.span_label == primary), None)
    if anchor is None:
        sys.exit(f"{spec.symbol}: no pinned leg spans the primary window "
                 f"{primary!r} -- the ΔE column needs it")
    a0 = _payload(cells, anchor)
    if a0 is None:
        sys.exit(f"{spec.symbol}: the anchor leg {anchor.address} has no "
                 "stored payload -- run the lane first "
                 "(regen_cells.py --group lane)")
    e0 = float(a0["total"])

    rows = []
    for k in keys:
        p = _payload(cells, k)
        route = "pinned" if k.fold_tag is not None else "free"
        head = [k.span_label, _cas(k), _fold_col(k, spec),
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
                f"{e * EV:.4f}", f"{(e - e0) * EV:+.4f}", _l2_text(p),
                _comm(p["comm_norm_h_lz"]), _comm(p["comm_norm_h_l2"])]
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

    text = src.read_text()
    if text.count(BEGIN) != 1 or text.count(END) != 1:
        sys.exit(f"{src}: expected exactly one {BEGIN}/{END} marker pair"
                 " (place them around the survey table on the first pass)")
    head, rest = text.split(BEGIN)
    _, tail = rest.split(END)
    new = head + new_block + tail

    # structural assertions (the TABLES protocol)
    assert new.count(BEGIN) == 1 and new.count(END) == 1
    body = [r for r in new[new.index(BEGIN):new.index(END)].splitlines()
            if r.startswith("|")]
    assert len(body) - 2 == len(keys), body
    for r, k in zip(body[2:], keys):
        assert r.startswith(f"| {k.span_label} |"), r

    tmp = src.with_suffix(".md.tmp")
    tmp.write_text(new)
    assert tmp.read_text() == new
    tmp.replace(src)
    n = sum(1 for k in keys if _payload(cells, k))
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
