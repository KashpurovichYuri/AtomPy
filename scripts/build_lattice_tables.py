"""Render an element's CELL-SOURCED window-lattice sections of
TABLES.md from its cells.json payloads (the cell wave's render
contract: the table and its provenance cannot diverge).

This renderer owns the sections whose every value is a cell payload;
log-sourced sections (Ti's 3d²4s² series, from run_tm_lattice.py --
the 0.9 mEh lattice-vs-campaign bookkeeping trace lives there) stay
hand-curated until their own cells run.

Value conventions (the lattice contract):

- E = payload ``e_ground`` -- the AVERAGED-fold reconstruction (the
  ruler), matching the historic lattice numbers' semantics in
  run_tm_lattice.run_tier. (The Methods tables instead show the run's
  own total -- finding 2; the two differ by the raw↔avg fold, ~2.8
  mEh on Ti CASCI cells.)
- Δ reading vs the ROHF anchor (cells.json ``rohf_energy``), 4-decimal
  Eh; ⟨L²⟩ 3 decimals; dets = payload ``n_det``.
- constructions per row: AOs (the m-complete shell span), PySCF MO
  (the canonical slice -- m-delocalized control; non-integer ⟨L²⟩ =
  the span does not close), NO (the full-span occupation-group prefix
  of the same source's term-averaged density, fact 23 -- equals AOs
  by construction; its closure certificate is the point), custom
  (declared empty this wave).

TABLES edit protocol lives here (as in build_methods_tables):
sections rebuilt from constructed blocks, structural assertions, temp
file + verified move.

Usage: python scripts/build_lattice_tables.py Ti
"""

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, "src")
from atomterm.cells import CellKey

SECTOR = (2, 0)          # d² Ms=1: the 4s² pair closed in the core
D2_SERIES = [            # (window label, CAS label)
    ("3d", "CAS(2,5)"),
    ("3d4p", "CAS(2,8)"),
    ("3d4p4d", "CAS(2,13)"),
]
HEAD = ("| window | CAS | dets | AOs, Eh | ⟨L²⟩ | PySCF MO, Eh | ⟨L²⟩ | "
        "NO, Eh | ⟨L²⟩ | custom, Eh | ⟨L²⟩ |")
SEP = "|---|---|---|---|---|---|---|---|---|---|---|"


def _e(p: dict) -> str:
    if not p:
        return "pending"
    return re.sub(r"(?<!e)-", "−", f"{p['e_ground']:.4f}")  # sign minus only


def _l2(p: dict) -> str:
    return f"{p['ground_l2']:.3f}" if p else ""


def main(atom: str) -> None:
    src = Path(f"elements/{atom}/TABLES.md")
    data = json.loads(Path(f"elements/{atom}/cells.json").read_text())
    cells = data["cells"]
    anchor = data["rohf_energy"]

    def payload(window, start):
        k = CellKey(atom, window, SECTOR, method="casci", start=start,
                    fold="avg")
        return cells.get(k.address)

    rows = []
    for w, cas in D2_SERIES:
        ao, mo = (payload(w, s) for s in ("ao", "mo"))
        no_label = {"3d": "3dno1", "3d4p": "3d4pno2", "3d4p4d": "3d4p4dno3"}[w]
        no = payload(no_label, "no")
        dets = f"{ao['n_det']:,}" if ao else "pending"
        rows.append(f"| {w} | {cas} | {dets} | {_e(ao)} | {_l2(ao)} | "
                    f"{_e(mo)} | {_l2(mo)} | {_e(no)} | {_l2(no)} | "
                    f"— | — |")

    section = [
        "## Window lattice — d² series",
        "",
        "CASCI, averaged; ROHF reference "
        + re.sub(r"(?<!e)-", "−", f"{anchor:.4f}")
        + " Eh. Constructions per row: the shell AOs"
        " (m-complete),",
        "the PySCF MO slice (non-integer ⟨L²⟩ = the span does not"
        " close), and",
        "the full-span NO group prefix of the same source (fact 23:"
        " equals the",
        "AOs by construction -- the certificate is the point)."
        " Source: `cells.json`",
        "(`*_2-0_*` addresses).",
        "",
        HEAD,
        SEP,
        *rows,
    ]
    new_section = "\n".join(section)

    text = src.read_text()
    b0, b1 = "## Window lattice — d² series", "## Window lattice — 3d²4s² series"
    assert text.count(b0) == 1 and text.count(b1) == 1, "section boundaries"
    new = text[:text.index(b0)] + new_section + "\n\n" + text[text.index(b1):]

    # structural assertions (the TABLES protocol)
    assert new.count(b0) == 1 and new.count(b1) == 1
    span = new[new.index(b0):new.index(b1)]
    body = [r for r in span.splitlines() if r.startswith("|")]
    assert len(body) - 2 == len(D2_SERIES), body
    for r, (w, _) in zip(body[2:], D2_SERIES):
        assert r.startswith(f"| {w} |"), r
        assert r.count("pending") <= 1, r      # a cell, never the row

    tmp = src.with_suffix(".md.tmp")
    tmp.write_text(new)
    assert tmp.read_text() == new
    tmp.replace(src)
    print(f"{src}: d² lattice rendered from cells.json "
          f"({len(rows)} rows, anchor {anchor:.4f})")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Ti")
