"""Render the ladder section of an element's TABLES.md from its
cells.json payloads (the build_*_tables family's ladder member; the
TABLES edit protocol verbatim -- structural assertions, the verified
.tmp move).

The block between the ``<!-- lane:begin/end -->`` markers is
GENERATED: one row per declared ladder leg
(``CaseSpec.lane_windows`` via ``cells.lane_keys``), as-run -- a leg
without a stored payload renders ``pending``.  Route/method columns
carry the treatment identity (``pinned``/``CI`` today; the free-route
rows arrive with the route axis).  The prose around the markers (the
intro, the readings) is HAND physics and stays untouched.

Usage:
    python scripts/build_lane_tables.py [atoms]

Atoms (default: every declared atom with ladder legs, CAMPAIGN order).
"""

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.cells import CellKey, lane_keys  # noqa: E402
from atomterm.windows import NOBLE_CORES, fold_shells  # noqa: E402

SHELL_COLS = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9, "h": 11}
EV = 27.2114
DASH = "—"
HEAD = ("| route | method | window | CAS | fold | dets | E, Eh | "
        "ΔE, mEh | E, eV | ΔE, eV | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | "
        "‖[H,L²]‖ |")
SEP = "|" + "---|" * 14
BEGIN, END = "<!-- lane:begin -->", "<!-- lane:end -->"


def minus(text: str) -> str:
    """Sign-minus normalization (the tables' typographic convention)."""
    return re.sub(r"(?<!e)-", "−", text)


def fold_label(tag: str) -> str:
    """The fold column's README notation: noble names capitalized,
    ``semi`` as itself, explicit lists braced."""
    if tag in NOBLE_CORES:
        return tag.capitalize()
    if tag == "semi":
        return tag
    return "{" + ",".join(fold_shells(tag)) + "}"


def cas_label(key: CellKey) -> str:
    """The CAS(nelec, ncas) label: electrons from the key's sector,
    columns from the window span's summed m-set sizes."""
    cols = sum(SHELL_COLS[l] for l in
               re.findall(r"\d+([spdfgh])", key.span_label))
    return f"({sum(key.sector)},{cols})"


def _l2_text(payload: dict) -> str:
    """⟨L²⟩ across the branch components: one value when they agree
    (the term reading), the min–max range when they do not (the
    scramble -- the ptp column carries its energy side)."""
    l2 = payload["l2"]
    if max(l2) - min(l2) < 1e-6:
        return f"{sum(l2) / len(l2):.4f}"
    return f"{min(l2):.1f}–{max(l2):.1f}"


def _comm(v) -> str:
    """One commutator norm, or the em-dash when the sector sat above
    the dense-meter cap (structural for the certified fold class)."""
    return f"{v:.1e}" if v == v else DASH          # NaN check: v != v


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
        p = payloads.get(k.address)
        if p is None or "void_reason" in p:
            rows.append(f"| pinned | CI | {k.span_label} | {cas_label(k)} "
                        f"| {fold_label(k.fold_tag)} | pending | | | | "
                        f"| | | | |")
            continue
        e = float(p["total"])
        de = (e - e0) * 1e3
        de_txt = "0" if k is anchor else f"{de:+.1f}"
        row = [f"{e:.4f}", de_txt, f"{e * EV:.4f}",
               f"{(e - e0) * EV:+.4f}", _l2_text(p),
               f"{p['ptp']:.0e}", _comm(p["comm_norm_h_lz"]),
               _comm(p["comm_norm_h_l2"])]
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
    new_block = section(spec, payloads)

    text = src.read_text()
    if text.count(BEGIN) != 1 or text.count(END) != 1:
        sys.exit(f"{src}: expected exactly one {BEGIN}/{END} marker pair"
                 " (place them around the ladder table on the first pass)")
    head, rest = text.split(BEGIN)
    _, tail = rest.split(END)
    new = head + new_block + tail

    # structural assertions (the TABLES protocol)
    assert new.count(BEGIN) == 1 and new.count(END) == 1
    body = [r for r in new[new.index(BEGIN):new.index(END)].splitlines()
            if r.startswith("|")]
    assert len(body) - 2 == len(keys), body
    for r, k in zip(body[2:], keys):
        assert r.startswith(f"| pinned | CI | {k.span_label} "), r

    tmp = src.with_suffix(".md.tmp")
    tmp.write_text(new)
    assert tmp.read_text() == new
    tmp.replace(src)
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
