"""Render an element's CELL-SOURCED window-lattice sections of
TABLES.md from its cells.json payloads (the cell wave's render
contract: the table and its provenance cannot diverge).

Usage: python scripts/build_lattice_tables.py Ti

Notes
-----
This renderer owns the sections whose every value is a cell payload;
log-sourced sections (Ti's 3d²4s² series, from the retired
run_tm_lattice.py -- git history; the 0.9 mEh lattice-vs-campaign
bookkeeping trace lives there) stay hand-curated until their own
cells run.

Value conventions (the lattice contract):

- E = ``energies_pyscf[0]`` -- the run's OWN total (as-run), directly
  comparable with the ROHF anchor (the variational reading: the MO
  rows sit exactly AT it; the AO rows are unbounded, their span does
  not contain the reference determinant). One convention per family:
  tables show as-run; the AVERAGED-fold reconstruction stays the
  analysis ruler the meters (commutators, ⟨L²⟩) are computed on --
  it sits ~2.8 mEh above as-run on Ti CASCI cells and would bend
  every anchor comparison by exactly that.
- ⟨L²⟩ 3 decimals; dets = payload ``n_det``; no ΔE column -- the
  anchor lives in the caption.
- constructions per row: AO (the m-complete shell span), MO (the
  canonical slice -- m-delocalized control; non-integer ⟨L²⟩ = the
  span does not close), custom (declared empty this wave). The NO
  construction is NOT a column here: the full-span group prefix of
  the same source equals the AO span by construction (fact 23) --
  the machinery certificate, not a methodology; the compression
  ladder from a bigger source lives in the 3d²4s² section.

TABLES edit protocol lives here (as in build_methods_tables):
sections rebuilt from constructed blocks, structural assertions, temp
file + verified move.
"""

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.cells import CellKey, lattice_keys  # noqa: E402

SHELL_COLS = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9, "h": 11}
HEAD = ("| window | CAS | dets | AO, Eh | ⟨L²⟩ | MO, Eh | ⟨L²⟩ | "
        "custom, Eh | ⟨L²⟩ |")
SEP = "|---|---|---|---|---|---|---|---|---|"
HEADING = "## Window lattice — 3d² series"
BOUNDARY = "## Window lattice — 3d²4s² series"


def cas_label(window: str) -> str:
    """The free-ion CAS label of a shell-span window: the electron
    count is the lattice's (2), the column count the summed m-set
    sizes of its shells."""
    cols = sum(SHELL_COLS[l] for l in re.findall(r"\d+([spdfgh])", window))
    return f"CAS(2,{cols})"


def _e(p: dict) -> str:
    """The as-run ground total (sign-minus normalized), or ``pending``
    when the cell has not run."""
    if not p:
        return "pending"
    return re.sub(r"(?<!e)-", "−",
                  f"{float(p['energies_pyscf'][0]):.4f}")


def _l2(p: dict) -> str:
    """The ground state's ⟨L²⟩ at table precision (blank when
    pending)."""
    return f"{p['ground_l2']:.3f}" if p else ""


def main(atom: str) -> None:
    """Rebuild ``<atom>``'s 3d² lattice section from cells.json under
    the TABLES edit protocol (structural assertions, verified move).
    The spans and sector derive from the declaration
    (``CaseSpec.lattice_windows`` via ``cells.lattice_keys``); this
    file owns presentation only."""
    spec = CAMPAIGN[atom]
    keys = lattice_keys(spec)
    if not keys:
        sys.exit(f"no declared lattice for {atom} "
                 f"(CaseSpec.lattice_windows is empty)")
    sector = keys[0].sector
    spans = spec.lattice_windows
    src = Path(f"elements/{atom}/TABLES.md")
    data = json.loads(Path(f"elements/{atom}/cells.json").read_text())
    cells = data["cells"]
    anchor = data["rohf_energy"]

    def payload(window, start):
        k = CellKey(atom, window, sector, method="casci", start=start,
                    fold="avg")
        return cells.get(k.address)

    rows = []
    for w in spans:
        ao, mo = (payload(w, s) for s in ("ao", "mo"))
        dets = f"{ao['n_det']:,}" if ao else "pending"
        rows.append(f"| {w} | {cas_label(w)} | {dets} | {_e(ao)} | {_l2(ao)} | "
                    f"{_e(mo)} | {_l2(mo)} | — | — |")

    section = [
        HEADING,
        "",
        "CASCI; E = the run's own total (as-run). ROHF reference "
        + re.sub(r"(?<!e)-", "−", f"{anchor:.4f}") + " Eh.",
        "Reading — README · Results (window lattice). Source:"
        " `cells.json` (`*_2-0_*`).",
        "",
        HEAD,
        SEP,
        *rows,
    ]
    new_section = "\n".join(section)

    text = src.read_text()
    # accept the pre-rename spelling on the first pass
    m = re.search(r"^## Window lattice — (?:3d²|d²) series$", text, re.M)
    assert m and text.count(BOUNDARY) == 1, "section boundaries"
    new = text[:m.start()] + new_section + "\n\n" + text[text.index(BOUNDARY):]

    # structural assertions (the TABLES protocol)
    assert new.count(HEADING) == 1 and new.count(BOUNDARY) == 1
    span = new[new.index(HEADING):new.index(BOUNDARY)]
    body = [r for r in span.splitlines() if r.startswith("|")]
    assert len(body) - 2 == len(spans), body
    for r, w in zip(body[2:], spans):
        assert r.startswith(f"| {w} |"), r
        assert r.count("pending") <= 1, r      # a cell, never the row

    tmp = src.with_suffix(".md.tmp")
    tmp.write_text(new)
    assert tmp.read_text() == new
    tmp.replace(src)
    print(f"{src}: 3d² lattice rendered from cells.json "
          f"({len(rows)} rows, as-run, anchor {anchor:.4f})")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Ti")
