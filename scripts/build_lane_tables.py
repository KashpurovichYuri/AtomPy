"""Render the ladder section of an element's TABLES.md from its
cells.json payloads (the build_*_tables family's ladder member; the
TABLES edit protocol verbatim -- structural assertions, the verified
.tmp move).

The block between the ``<!-- lane:begin/end -->`` markers is
GENERATED: one row per declared ladder leg
(``CaseSpec.lane_windows`` via ``cells.lane_keys``), as-run -- a leg
without a stored payload renders ``pending``.  Route/method columns
carry the treatment identity (``pinned``/``CI``).  An atom whose
ladder rows render in its CI SURVEY (the survey markers present) is
skipped -- the survey supersedes this section.  The prose around the
markers (the intro, the readings) is HAND physics and stays
untouched.

Usage:
    python scripts/build_lane_tables.py [atoms]

Atoms (default: every declared atom with ladder legs, CAMPAIGN order).
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.cells import lane_keys  # noqa: E402
from _tables_common import (EV, cas_label, comm_cell,  # noqa: E402
                            fold_label, l2_text, minus, payload_of,
                            replace_block)

HEAD = ("| route | method | window | CAS | fold | dets | E, Eh | "
        "ΔE, mEh | E, eV | ΔE, eV | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | "
        "‖[H,L²]‖ |")
SEP = "|" + "---|" * 14
BEGIN, END = "<!-- lane:begin -->", "<!-- lane:end -->"
SURVEY_BEGIN = "<!-- survey:begin -->"


def section(spec, payloads: dict) -> str:
    """The generated block between the markers: one row per declared
    leg, as-run, in declaration order.  ``ΔE`` anchors on the leg
    whose span is the case's primary window (asserted present -- the
    column's zero)."""
    keys = lane_keys(spec)
    primary = "".join(spec.shell_patterns)
    anchor = next((k for k in keys if k.span_label == primary), None)
    if anchor is None:
        sys.exit(f"{spec.symbol}: no ladder leg spans the primary window "
                 f"{primary!r} -- the ΔE column needs it")
    a0 = payloads.get(anchor.address)
    if a0 is None:
        sys.exit(f"{spec.symbol}: the anchor leg {anchor.address} has no "
                 "stored payload -- run the lane first "
                 "(regen_cells.py --group lane)")
    e0 = float(a0["total"])

    rows = []
    for k in keys:
        p = payload_of(payloads, k)
        if p is None:
            rows.append(f"| pinned | CI | {k.span_label} | {cas_label(k)} "
                        f"| {fold_label(k.fold_tag)} | pending | | | | "
                        f"| | | | |")
            continue
        e = float(p["total"])
        de = (e - e0) * 1e3
        de_txt = "0" if k is anchor else f"{de:+.1f}"
        row = [f"{e:.4f}", de_txt, f"{e * EV:.4f}",
               f"{(e - e0) * EV:+.4f}", l2_text(p),
               f"{p['ptp']:.0e}",
               comm_cell(p, "comm_norm_h_lz", "res_h_lz"),
               comm_cell(p, "comm_norm_h_l2", "res_h_l2")]
        rows.append("| pinned | CI | "
                    + " | ".join([k.span_label, cas_label(k),
                                  fold_label(k.fold_tag),
                                  f"{p['n_det']:,}"]
                                 + [minus(v) for v in row])
                    + " |")
    return "\n".join([BEGIN, HEAD, SEP, *rows, END])


def run(atom: str) -> None:
    """Rebuild ``elements/<Atom>/TABLES.md``'s ladder block from its
    cells.json payloads under the TABLES edit protocol."""
    spec = CAMPAIGN[atom]
    keys = lane_keys(spec)
    if not keys:
        sys.exit(f"no declared ladder legs for {atom}")
    src = Path(f"elements/{atom}/TABLES.md")
    payloads = json.loads(
        Path(f"elements/{atom}/cells.json").read_text())["cells"]
    if SURVEY_BEGIN in src.read_text():
        print(f"{src}: the ladder rows render in the CI survey -- skipped")
        return
    new_block = section(spec, payloads)
    replace_block(src, BEGIN, END, new_block)

    # structural assertions (the TABLES protocol)
    body = [r for r in new_block.splitlines() if r.startswith("|")]
    assert len(body) - 2 == len(keys), body
    for r, k in zip(body[2:], keys):
        assert r.startswith(f"| pinned | CI | {k.span_label} "), r
    print(f"{src}: ladder section rendered from cells.json "
          f"({len(keys)} legs, as-run)")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("atoms", nargs="*",
                   help="declared atoms with ladder legs (default: all)")
    ns = p.parse_args(argv)
    chosen = ns.atoms or [a for a in CAMPAIGN if lane_keys(CAMPAIGN[a])]
    unknown = [a for a in chosen if a not in CAMPAIGN]
    if unknown:
        sys.exit(f"not declared: {', '.join(unknown)}")
    for a in chosen:
        if not lane_keys(CAMPAIGN[a]):
            print(f"{a}: no declared ladder legs -- skipped")
            continue
        run(a)


if __name__ == "__main__":
    main()
