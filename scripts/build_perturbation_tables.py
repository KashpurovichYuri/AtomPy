"""Render PERTURBATION_TABLE.md (+ the CSV mirror) from
``elements/<Atom>/perturbation.json`` -- the render contract: the
table renders ONLY from the payload (byte-determinism: no render-time
stamps; the provenance stamp is the payload's digest, so a re-render
of the same payload is byte-identical), one table between markers,
and a SCHEMA refusal for stale payloads (the stale-render precedent:
an old-schema payload refuses with the collector hint -- the payload
is regenerated, never hand-patched).

Usage: python scripts/build_perturbation_tables.py <atom> [<atom> ...]

Groups branch on key IDENTITY only (the methods-table rule): the span
ladder (casci, no/no1 absent, no frozen mask), the NO twins (casci +
no start), the reference axis (casscf / sa-casscf), the frozen pair
(frozen mask).  Every stored root enters the table through its cell's
row: root 0's values, the block spread beside them (the declared
roots policy lives in the payload; a degenerate block renders one
line with spread ~0, a split block carries its spread).
"""

import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _tables_common import atomic_write_text, replace_block  # noqa: E402
from run_pt2_ladder import _declared, _declared_lane_cells  # noqa: E402

BEGIN = "<!-- pt2 begin -->"
END = "<!-- pt2 end -->"

HEAD = ("| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) "
        "| E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |")
SEP = "|---" * 9 + "|"


def _required_schema(payload: dict) -> None:
    """Refuse a stale/old-schema payload (the stale-render precedent)."""
    prov = payload.get("provenance", {})
    need_prov = ("engine", "pyscf_version", "collected", "date",
                 "declared_roots_policy", "e_ref_partition",
                 "wobble_note", "stored_roots_bands", "meters")
    missing = [k for k in need_prov if k not in prov]
    cells = payload.get("cells")
    if missing or not isinstance(cells, list) or not cells:
        sys.exit(f"stale payload for {payload.get('atom')!r}: missing "
                 f"provenance keys {missing or 'cells'} -- regenerate "
                 "with scripts/run_pt2_ladder.py")
    for c in cells:
        if c.get("disposition") not in ("rows", "refusal"):
            sys.exit(f"stale payload: cell {c.get('address')!r} has no "
                     "disposition -- regenerate")
        if c.get("lane") and "blocks" not in c:
            sys.exit(f"stale payload: lane cell {c.get('address')!r} "
                     "carries no blocks -- regenerate")


def _group(c: dict, atom: str) -> str:
    """The row group, by key identity only (the group vocabulary is
    the manifest's -- one home)."""
    addr = c["address"]
    if c.get("lane") or any(a == addr
                            for a, _ in _declared_lane_cells(atom)):
        return ("the pinned lane (the core-partition axis; the D2h "
                "components)")
    if "_frozen_" in addr:
        return "the frozen pair (the core-relaxation disclosure)"
    if "_casci_no_-" in addr:
        return "the NO twins (CASCI)"
    if "_casscf" in addr:
        return "the reference axis (single-root CASSCF vs SA)"
    size = _declared(atom).get("size", {})
    if any(addr.startswith(f"{w}_") for w in size):
        return "the CAS-size axis (the flagship R-test)"
    size_group = _declared(atom).get("size_group", ())
    if any(addr.startswith(f"{w}_") for w in size_group):
        return "the CAS-size axis (the flagship R-test)"
    dist = _declared(atom).get("dist", ())
    if any(addr.startswith(f"{w}_") for w in dist):
        return "the distilled prior family (CASCI ao)"
    return "the span ladder (CASCI: the unoptimized windows)"


def _rows(payload: dict) -> list[str]:
    atom = payload["atom"]
    # GROUP-COLLECTED: one section per group, the cells in payload
    # order within it -- a window-interleaved declared set must not
    # scatter one axis across repeating headers (the U14e flagship
    # read the trio had to be reassembled by hand)
    sections: dict[str, list[str]] = {}
    order: list[str] = []
    for c in payload["cells"]:
        g = _group(c, atom)
        if g not in sections:
            sections[g] = []
            order.append(g)
        sections[g].append(_cell_row(c))
    out = []
    for g in order:
        out += ["", f"**{g}**", "", HEAD, SEP]
        out += sections[g]
    return out


def _cell_row(c: dict) -> str:
    """One cell's row(s): the refusal line, the lane rows (one per
    block root), or the ground-block row."""
    if c["disposition"] == "refusal":
        reason = c["refusal_reason"].split("--")[0].strip()
        return (f"| `{c['address']}` | — | — | — | REFUSED: "
                f"{reason} | — | — | — | — |")
    if c.get("lane"):
        # the lane rows: one per block root; the spread column reads
        # the TERM's across-component E(2) ptp
        spread = c.get("term_e2_ptp", float("nan")) * 1e3
        out = []
        for blk in c["blocks"]:
            r0 = blk["roots"][0]
            out.append(
                f"| `{c['address']}` ·{blk['irrep']} "
                f"| {blk['mult']} "
                f"| {r0['l2']:.6f} "
                f"| {r0['e_ref']:.8f} "
                f"| {r0['e_corr'] * 1e3:+.4f} "
                f"| {r0['e_corr_eV']:+.4f} "
                f"| {r0['e_total']:.8f} "
                f"| {r0['e_total_eV']:+.4f} "
                f"| {spread:.1e} |")
        return "\n".join(out)
    roots = c["roots"]
    r0 = roots[0]
    # the spread reads the GROUND BLOCK only: the roots degenerate
    # with root 0 (the gauge block); broken spans store excited
    # roots beside the ground one and those are different states,
    # not gauges of it
    block = [r for r in roots
             if abs(r["e_ref"] - r0["e_ref"]) <= 1e-6]
    e2 = [r["e_corr"] for r in block]
    spread = (max(e2) - min(e2)) * 1e3
    return (
        f"| `{c['address']}` | {c['nroots']} "
        f"| {r0['l2']:.6f} "
        f"| {r0['e_ref']:.8f} "
        f"| {r0['e_corr'] * 1e3:+.4f} "
        f"| {r0['e_corr_eV']:+.4f} "
        f"| {r0['e_total']:.8f} "
        f"| {r0['e_total_eV']:+.4f} "
        f"| {spread:.1e} |")


def _csv(payload: dict) -> str:
    atom = payload["atom"]
    lines = ["address,nroots,l2,e_ref_eh,e2_meh,e2_ev,e_pt2_eh,"
             "e_pt2_ev,spread_meh,disposition"]
    for c in payload["cells"]:
        if c["disposition"] == "refusal":
            lines.append(f"{c['address']},,,,,,,,"
                         f",refused")
            continue
        if c.get("lane"):
            for blk in c["blocks"]:
                r0 = blk["roots"][0]
                lines.append(
                    f"{c['address']}·{blk['irrep']},{blk['mult']},"
                    f"{r0['l2']:.6f},{r0['e_ref']:.8f},"
                    f"{r0['e_corr'] * 1e3:+.4f},{r0['e_corr_eV']:+.4f},"
                    f"{r0['e_total']:.8f},{r0['e_total_eV']:+.4f},"
                    f"{c.get('term_e2_ptp', 0.0) * 1e3:.1e},rows")
            continue
        r0 = c["roots"][0]
        block = [r for r in c["roots"]
                 if abs(r["e_ref"] - r0["e_ref"]) <= 1e-6]
        e2 = [r["e_corr"] for r in block]
        lines.append(
            f"{c['address']},{c['nroots']},{r0['l2']:.6f},"
            f"{r0['e_ref']:.8f},{r0['e_corr'] * 1e3:+.4f},"
            f"{r0['e_corr_eV']:+.4f},{r0['e_total']:.8f},"
            f"{r0['e_total_eV']:+.4f},"
            f"{(max(e2) - min(e2)) * 1e3:.1e},rows")
    return "\n".join(lines) + "\n"


def render_atom(atom: str, workdir: str = "elements") -> None:
    payload_path = Path(workdir) / atom / "perturbation.json"
    if not payload_path.exists():
        sys.exit(f"no {payload_path} -- collect first "
                 "(scripts/run_pt2_ladder.py)")
    payload = json.loads(payload_path.read_text())
    _required_schema(payload)
    digest = hashlib.sha256(payload_path.read_bytes()).hexdigest()
    prov = payload["provenance"]
    block = "\n".join([
        BEGIN,
        f"PT2 analyzation: {atom} -- {prov['engine']}, "
        f"pyscf {prov['pyscf_version']}; payload collected "
        f"{prov['collected']}.",
        "",
        "Reading key: README · Analysis (the reference-class axis, "
        "the gauge-spread disclosure, the basis-thinness caveat). "
        "E_ref within its own partition (kernel totals only); "
        "E_PT2 = E_ref + E(2).  The row reads the GROUND root; the "
        "block spread runs over the roots degenerate with it (the "
        "gauge block) -- broken spans store excited roots beside it "
        "(different states, not gauges; full per-root records live "
        "in the payload).",
        "",
    ] + _rows(payload) + [
        "",
        f"<!-- payload sha256 {digest[:16]}... "
        f"bands: {prov['stored_roots_bands']} -->",
        END,
    ]) + "\n"
    table = Path(workdir) / atom / "PERTURBATION_TABLE.md"
    if not table.exists():
        atomic_write_text(table, f"# {atom} -- PT2\n\n" + block)
    else:
        replace_block(table, BEGIN, END, block)
    atomic_write_text(Path(workdir) / atom / "PERTURBATION_TABLE.csv",
                      _csv(payload))
    print(f"rendered {table} + PERTURBATION_TABLE.csv "
          f"(payload sha256 {digest[:16]}...)")


def main() -> None:
    atoms = sys.argv[1:]
    if not atoms:
        sys.exit("usage: build_perturbation_tables.py <atom> [...]")
    for atom in atoms:
        render_atom(atom)


if __name__ == "__main__":
    main()
