"""Render an element's CELL-SOURCED window-lattice section of
TABLES.md from its cells.json payloads (the build_*_tables family's
lattice member; the render contract: the table and its provenance
cannot diverge).

Usage: python scripts/build_lattice_tables.py [atoms]

The whole section -- heading, tag line, table -- is GENERATED from
the declaration (``CaseSpec.lattice_windows`` at ``lattice_sector``
via ``cells.lattice_keys``), as-run: a cell without a stored payload
renders ``pending``.  Column conventions and the reading live in
README -- the one home; this file owns presentation only.

TABLES edit protocol lives here (as in build_methods_tables): the
section is rebuilt between its heading and the next heading,
structural assertions checked, temp file + verified move.  First pass
accepts a legacy window-lattice heading and rewrites it.
"""

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.cells import lattice_keys  # noqa: E402
from _tables_common import SHELL_COLS, minus  # noqa: E402

HEAD = "| window | CAS | dets | AO, Eh | ⟨L²⟩ | MO, Eh | ⟨L²⟩ |"
SEP = "|" + "---|" * 7
HEADING = "## Window lattice"
# first pass: accept any legacy lattice heading (the pre-grammar spellings)
_LEGACY = re.compile(r"^## (?:Production — )?[Ww]indow lattice.*$", re.M)


def _cas(cols: int, sector: tuple[int, int]) -> str:
    """The CAS(nelec, ncas) label: electrons from the declared sector."""
    return f"CAS({sum(sector)},{cols})"


def _cell(cells: dict, key) -> dict | None:
    """The key's payload, or None when the cell has not run."""
    p = cells.get(key.address)
    if p is None or "void_reason" in p:
        return None
    return p


def _e(p: dict | None) -> str:
    """The as-run ground total, or ``pending``."""
    return minus(f"{float(p['energies_pyscf'][0]):.4f}") if p else "pending"


def _l2(p: dict | None) -> str:
    """The ground state's ⟨L²⟩ at table precision (blank when pending)."""
    return f"{p['ground_l2']:.3f}" if p else ""


def section(spec, cells: dict) -> str:
    """The generated section (heading, tag line, table): one row per
    declared span, as-run, in declaration order; the two
    constructions (AO / MO) are its column pairs."""
    keys = lattice_keys(spec)
    sector = keys[0].sector
    by_span: dict[str, dict[str, object]] = {}
    for k in keys:
        by_span.setdefault(k.span_label, {})[k.start] = _cell(cells, k)

    rows = []
    for w in spec.lattice_windows:
        parts = by_span[w]
        ao = parts.get("ao")
        dets = f"{ao['n_det']:,}" if ao else "pending"
        cols = sum(SHELL_COLS[l] for l in
                   re.findall(r"\d+([spdfgh])", w))
        row = [w, _cas(cols, sector), dets]
        for start in ("ao", "mo"):
            p = parts.get(start)
            row += [_e(p), _l2(p)]
        rows.append("| " + " | ".join(row) + " |")

    return "\n".join([
        HEADING,
        "",
        f"Route free · partition averaged fold · sector {sector} · "
        "source `cells.json` (as-run).",
        "",
        HEAD,
        SEP,
        *rows,
    ])


def run(atom: str) -> None:
    """Rebuild ``<atom>``'s lattice section from cells.json under the
    TABLES edit protocol.  The spans, sector, and every cell address
    derive from the declaration; this file owns presentation only."""
    spec = CAMPAIGN[atom]
    keys = lattice_keys(spec)
    if not keys:
        sys.exit(f"no declared lattice for {atom} "
                 "(CaseSpec.lattice_windows is empty)")
    src = Path(f"elements/{atom}/TABLES.md")
    cells = json.loads(
        Path(f"elements/{atom}/cells.json").read_text())["cells"]
    text = src.read_text()
    if HEADING not in text and "<!-- survey:begin -->" in text:
        # the survey supersedes this section (the lattice spans are
        # its first window family)
        print(f"{src}: the lattice rows render in the CI survey -- skipped")
        return
    new_section = section(spec, cells)

    if HEADING in text:
        start = text.index(HEADING)
    else:
        m = _LEGACY.search(text)
        if m is None:
            sys.exit(f"{src}: no window-lattice heading found -- place "
                     f"'{HEADING}' on the first pass")
        start = m.start()
    end = text.index("\n## ", start) + 1 if "\n## " in text[start:] \
        else len(text)
    new = text[:start] + new_section + "\n\n" + text[end:]

    # structural assertions (the TABLES protocol), scoped to the
    # section (heading to the next heading -- later sections carry
    # their own tables)
    assert new.count(HEADING) == 1
    sect = new[new.index(HEADING):
               new.index("\n## ", new.index(HEADING)) + 1
               if "\n## " in new[new.index(HEADING):] else len(new)]
    body = [r for r in sect.splitlines() if r.startswith("|")]
    assert len(body) - 2 == len(spec.lattice_windows), body
    for r, w in zip(body[2:], spec.lattice_windows):
        assert r.startswith(f"| {w} |"), r

    tmp = src.with_suffix(".md.tmp")
    tmp.write_text(new)
    assert tmp.read_text() == new
    tmp.replace(src)
    n = sum(1 for k in keys if _cell(cells, k))
    print(f"{src}: lattice rendered from cells.json "
          f"({len(spec.lattice_windows)} spans, {n}/{len(keys)} cells run)")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("atoms", nargs="*",
                   help="declared atoms with lattices (default: all)")
    ns = p.parse_args(argv)
    chosen = ns.atoms or [a for a in CAMPAIGN if lattice_keys(CAMPAIGN[a])]
    unknown = [a for a in chosen if a not in CAMPAIGN]
    if unknown:
        sys.exit(f"not declared: {', '.join(unknown)}")
    for a in chosen:
        if not lattice_keys(CAMPAIGN[a]):
            print(f"{a}: no declared lattice -- skipped")
            continue
        run(a)


if __name__ == "__main__":
    main()
