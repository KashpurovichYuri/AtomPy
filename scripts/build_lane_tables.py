"""Render the pinned-lane section of an element's TABLES.md from its
cells.json payloads (the build_*_tables family's lane member; the
TABLES edit protocol verbatim -- structural assertions, the verified
.tmp move).

The block between the ``<!-- lane:begin/end -->`` markers is
GENERATED: one row per declared lane leg (``CaseSpec.lane_windows``
via ``cells.lane_keys``), as-run -- a leg without a stored payload
renders ``pending``.  The prose around the markers (the intro, the
readings, the mechanism record) is HAND physics and stays untouched.

Usage:
    python scripts/build_lane_tables.py [atoms]

Atoms (default: every declared atom with lane legs, CAMPAIGN order).
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
L_LETTER = "SPDFGHIKLMNOQRTUVWXYZ"
HEAD = "| window | CAS | fold | dets | E, Eh | vs declared | term |"
SEP = "|---|---|---|---|---|---|---|"
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


def term_verdict(payload: dict, spec) -> str:
    """The generated term verdict from the meters: exact when every
    branch component reads L(L+1) and S(S+1) and the components sit
    degenerate; otherwise NOT a term window (the scramble's nuance
    lives in the hand readings)."""
    L = spec.term[0]
    name = term_name(spec)
    if (all(abs(v - L * (L + 1)) < 1e-6 for v in payload["l2"])
            and all(abs(v - S_of(spec)) < 1e-6 for v in payload["s2"])
            and payload["ptp"] < 1e-8):
        return f"{name} exact"
    return f"not a {name} window"


_SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹"


def term_name(spec) -> str:
    """The term symbol of the declared ground term (³F for 2S=2, L=3)."""
    return _SUP[spec.spin + 1] + L_LETTER[spec.term[0]]


def S_of(spec) -> float:
    """S(S+1) of the declared spin (2S = spin)."""
    s = spec.spin / 2
    return s * (s + 1)


def section(spec, payloads: dict) -> str:
    """The generated block between the markers: one row per declared
    leg, as-run, in declaration order.  ``vs declared`` anchors on the
    leg whose span is the case's primary window (asserted present --
    the column's zero)."""
    keys = lane_keys(spec)
    primary = "".join(spec.shell_patterns)
    anchor = next((k for k in keys if k.span_label == primary), None)
    if anchor is None:
        sys.exit(f"{spec.symbol}: no lane leg spans the primary window "
                 f"{primary!r} -- the 'vs declared' column needs it")
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
            rows.append(f"| {k.span_label} | {cas_label(k)} | "
                        f"{fold_label(k.fold_tag)} | pending | | | |")
            continue
        e_txt = minus(f"{float(p['total']):.4f}")
        de = (float(p["total"]) - e0) * 1e3
        vs = "0" if k is anchor else f"{de:+.1f} mEh"
        rows.append(
            f"| {k.span_label} | {cas_label(k)} | {fold_label(k.fold_tag)}"
            f" | {p['n_det']:,} | {e_txt} | {minus(vs)} |"
            f" {term_verdict(p, spec)} |")
    return "\n".join([BEGIN, HEAD, SEP, *rows, END])


def run(atom: str) -> None:
    """Rebuild ``elements/<Atom>/TABLES.md``'s lane block from its
    cells.json payloads under the TABLES edit protocol."""
    spec = CAMPAIGN[atom]
    keys = lane_keys(spec)
    if not keys:
        sys.exit(f"no declared lane legs for {atom}")
    src = Path(f"elements/{atom}/TABLES.md")
    payloads = json.loads(
        Path(f"elements/{atom}/cells.json").read_text())["cells"]
    new_block = section(spec, payloads)

    text = src.read_text()
    if text.count(BEGIN) != 1 or text.count(END) != 1:
        sys.exit(f"{src}: expected exactly one {BEGIN}/{END} marker pair"
                 " (place them around the lane table on the first pass)")
    head, rest = text.split(BEGIN)
    _, tail = rest.split(END)
    new = head + new_block + tail

    # structural assertions (the TABLES protocol)
    assert new.count(BEGIN) == 1 and new.count(END) == 1
    body = [r for r in new[new.index(BEGIN):new.index(END)].splitlines()
            if r.startswith("|")]
    assert len(body) - 2 == len(keys), body
    for r, k in zip(body[2:], keys):
        assert r.startswith(f"| {k.span_label} "), r

    tmp = src.with_suffix(".md.tmp")
    tmp.write_text(new)
    assert tmp.read_text() == new
    tmp.replace(src)
    print(f"{src}: lane section rendered from cells.json "
          f"({len(keys)} legs, as-run)")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("atoms", nargs="*",
                   help="declared atoms with lane legs (default: all)")
    ns = p.parse_args(argv)
    chosen = ns.atoms or [a for a in CAMPAIGN if lane_keys(CAMPAIGN[a])]
    unknown = [a for a in chosen if a not in CAMPAIGN]
    if unknown:
        sys.exit(f"not declared: {', '.join(unknown)}")
    for a in chosen:
        if not lane_keys(CAMPAIGN[a]):
            print(f"{a}: no declared lane legs -- skipped")
            continue
        run(a)


if __name__ == "__main__":
    main()
