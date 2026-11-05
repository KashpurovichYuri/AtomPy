"""Render the method-ladder section of an element's TABLE.md from its
cells.json payloads -- the key-addressed provenance (the render
contract: the table and its source cannot diverge).

Usage: python scripts/build_methods_tables.py [atom]

THE ONE LADDER TABLE (the element file IS the table): every run cell
of the atom a row, ROHF anchor first, the UHF anchor beside it when
the store carries the record; families in declaration order,
separated by the rule row -- the lattice spans (CI only), the
primary window (each fold group's CASCI BESIDE its optimizer
variants -- CI interleaved, the row-collapse policy is DELETED),
then every lane window whose route twins ran (the pinned and free
CI, the free CASSCF/SA, the pinned SA).  NO-start cells are not
ladder rows.  The PLANNED windows follow as pointers when gated
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
  fold) and on the ROW POLICY; order and grouping live once in the
  declaration builder (cells.window_grid / cells.route_keys) -- this
  renderer receives ready-made key lists and formats.

Applicability: any declared atom (``cases.CAMPAIGN``; windows and
cell addresses derive from the declaration -- this file owns
presentation only: headings, notes, value formatting).

ladder edit protocol lives here: sections are rebuilt from constructed
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
from atomterm.cells import (grid_keys, lattice_keys,  # noqa: E402
                            route_keys, window_grid)
from _tables_common import (DASH, EV, cas_label,  # noqa: E402
                            comm_cell, conv_cell, fold_col,
                            l2_text, minus, payload_of, replace_block)
BOUNDARY = "## Studies (frozen history; ticket 28)"
POINTER = ("Declared, not yet run -- the optimizer cells are void "
           "until their CASCI anchors are read (the anchor-first "
           "scenario). One regen away: `regen_cells.py {atom} "
           "--window {window}`.")

# the ladder grammar (one shape with the CAS ladder table); the
# family separator row renders as a horizontal rule mid-table
HEAD = ("| method | route | window | CAS | fold | dets | E, Eh | "
        "ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | "
        "‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |")
SEP = "|" + "---|" * 18
RULE = "|" + " --- |" * 18
BEGIN, END = "<!-- ladder:begin -->", "<!-- ladder:end -->"

# identity -> display (the ONLY branching the renderer does)
METHOD_NAMES = {"casci": "CASCI", "casscf": "CASSCF",
                "sa-casscf": "SA-CASSCF"}
START_NAMES = {"ao": "AO", "mo": "MO", "no": "NO"}
MASK_NAMES = {None: "free", "radial": "radial", "frozen": "frozen",
              "virtmix": "virtual-mix"}


CSV_HEADER = ["method", "route", "window", "CAS", "fold", "dets", "E_Eh",
              "dE_mEh", "E_eV", "dE_eV", "conv", "L2", "ptp_Eh", "PLz",
              "PLperp", "HLz", "HL2", "note"]


def csv_of(body: str) -> str:
    """The CSV mirror of the ladder table: the same data rows (the
    rule rows and the markdown furniture dropped), ASCII header
    names, csv-quoted (the braced folds carry commas)."""
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
    start rides in the window cell, the substrate in the fold
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


def _row(k, p, anchor, spec) -> str:
    """One run cell in the ladder columns -- every run cell a row,
    fully identified (method, route, window/start, CAS, fold)."""
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
    cyc, cap = p.get("cycles"), p.get("max_cycle")
    # the cap-hitting flag by IDENTITY on payloads carrying the cap
    # (the store's max_cycle field -- the run's actual mc cap);
    # payloads predating the field keep the legacy == 200 key (the
    # cap the measurement runs declared; PySCF's default is 50)
    note = "cycle cap" if cyc and cyc == (cap if cap is not None else 200) \
        else ""
    cells = [
        _method(k), route, f"{k.span_label} / {START_NAMES[k.start]}",
        cas_label(k), fold, f"{p['n_det']:,}",
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


def _rows(keys, cells, anchor, spec) -> list[str]:
    """Every RUN cell of the key set except NO starts (the
    occupation-group windows -- diagnostics, not ladder rows), in
    the given order."""
    return [_row(k, p, anchor, spec)
            for k in keys
            if k.start != "no"
            and (p := payload_of(cells, k)) is not None]


def _anchor_row(anchor: float) -> str:
    """The ROHF anchor row: ΔE's zero (the caption states it)."""
    cells = ["ROHF (anchor)", "—", "—", "—", "—", "—",
             f"{anchor:.4f}", "0", f"{anchor * EV:.4f}", "+0.0000",
             DASH, DASH, DASH, DASH, DASH, DASH, DASH, "forced filling"]
    for i in (6, 8):
        cells[i] = minus(cells[i])
    return "| " + " | ".join(cells) + " |"


def _uhf_anchor_row(anchor: float, uhf: dict) -> str:
    """The UHF anchor row: the spin-broken reference beside the ROHF
    one -- ΔE against the same zero, the meter reading (⟨S²⟩, the
    weakest corresponding-orbital pair) in the note."""
    e = uhf["e_tot"]
    de = (e - anchor) * 1000.0
    cells = ["UHF (anchor)", "—", "—", "—", "—", "—",
             f"{e:.4f}", f"{de:+.1f}", f"{e * EV:.4f}",
             f"{de * EV / 1000.0:+.4f}",
             DASH, DASH, DASH, DASH, DASH, DASH, DASH,
             f"⟨S²⟩={uhf['s2']:.4f}, σ_min={uhf['sigma_min']:.4f}"]
    # sign minus on the ENERGY cells only (names keep their hyphens)
    for i in (6, 7, 8, 9):
        cells[i] = minus(cells[i])
    return "| " + " | ".join(cells) + " |"


def ladder_section(spec, cells, anchor, uhf: dict | None = None) -> list[str]:
    """THE LADDER: one table for the whole atom.  Families in
    declaration order, separated by the rule row -- the lattice
    spans (CI only), the primary window (each fold group's CASCI
    BESIDE its optimizer variants), the planned windows whose CASCI
    anchors are read (the anchor-first protocol's readable half: the
    tier families render their anchor rows while the optimizer cells
    stay void at the key, ``--full`` materializing them into the same
    families), then every lane window whose route twins ran (the
    pinned and free CI, the free CASSCF/SA, the pinned SA).  The ROHF
    anchor row leads, the UHF anchor row follows when the store
    carries the record; NO-start cells are not ladder rows; a lane
    family without run optimizer twins is not this table's story
    yet."""
    body: list[str] = []
    seen: set[str] = set()

    def _fam(keys) -> None:
        rows = _rows([k for k in keys
                      if k.address not in seen and k.start != "no"],
                     cells, anchor, spec)
        seen.update(k.address for k in keys)
        if rows:
            body.extend([RULE] if body else [])
            body.extend(rows)

    for span in dict.fromkeys(k.window for k in lattice_keys(spec)):
        _fam([k for k in lattice_keys(spec) if k.window == span])

    # the primary family: fold groups in declaration order, each
    # group's CASCI beside its optimizer variants
    grid = [k for k in grid_keys(spec) if k.start != "no"]
    groups: dict[str, list] = {}
    for k in grid:
        groups.setdefault(fold_col(k, spec), []).append(k)
    for keys in groups.values():
        _fam(keys)

    routes = route_keys(spec)
    for leg in [k for k in routes
                if k.substrate == "pur" and k.method == "sa-casscf"]:
        twins = [k for k in routes if k.window == leg.window]
        if any(k.method != "casci" and k.start != "no"
               and payload_of(cells, k) for k in twins):
            _fam([k for k in twins if k.start != "no"])

    # the planned windows whose CASCI anchors are read: the tier
    # families (the whole window's non-NO cells at once -- the optimizer
    # cells are void at the key under the anchor-first gate and join
    # the family when --full materializes them)
    for w in spec.planned_windows:
        _fam([k for k in window_grid(spec, w) if k.start != "no"])

    out = [BEGIN, HEAD, SEP, _anchor_row(anchor)]
    if uhf is not None:
        out.append(_uhf_anchor_row(anchor, uhf))
    out.extend([*body, END])

    for w in spec.planned_windows:
        if not any(payload_of(cells, k) for k in window_grid(spec, w)):
            out[-1:-1] = ["", POINTER.format(atom=spec.key, window=w)]
    return out


def main(atom: str) -> None:
    """Rebuild ``<atom>``'s ladder block (the one table) from
    cells.json under the ladder edit protocol (the marker pair;
    structural checks; verified move) and write the CSV mirror."""
    src = Path(f"elements/{atom}/TABLE.md")
    if not src.exists():
        sys.exit(f"{src}: not a one-table atom (no TABLE.md) -- the "
                 "atom's pass migrates it")
    data = json.loads(Path(f"elements/{atom}/cells.json").read_text())
    cells, anchor = data["cells"], data["rohf_energy"]
    spec = campaign_spec(atom)

    body = "\n".join(ladder_section(spec, cells, anchor,
                                    uhf=data.get("uhf")))

    # structural assertions (the ladder protocol): no pending rows,
    # no starred gauge numbers, every table row well-formed
    assert "pending" not in body
    assert not re.search(r"\d\*", body)      # no starred gauge numbers
    for line in body.splitlines():
        if line.startswith("|"):
            assert line.count("|") == HEAD.count("|"), line

    replace_block(src, BEGIN, END, body)
    src.with_suffix(".csv").write_text(csv_of(body))
    n = sum(1 for line in body.splitlines()
            if line.startswith(("| CASCI", "| CASSCF", "| SA-CASSCF")))
    print(f"{src}: ladder rendered ({n} run rows) + TABLE.csv mirror")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Ti")
