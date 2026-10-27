"""Render the Methods-compared factor tables of an element's TABLES.md
from its cells.json payloads -- the key-addressed provenance (the cell
wave's render contract: the table and its source cannot diverge).

One table per declared window: ONE grid with an explicit method column
(CASCI / CASSCF / SA-CASSCF), an em dash for the CASCI non-axis (the
mask; CASCI has no optimizer), blank factor cells inheriting from the
row above, and ``pending`` for declared-but-not-run cells.

Value conventions (the manifest contract, README · Methods):

- E = the run's OWN total (``energies_pyscf[0]``), Eh + eV (×27.2114)
  -- the partition-changed class never shows our reconstruction
  (finding 2); ΔE vs the ROHF anchor (``rohf_energy``); the anchor
  is also stated in each caption;
- ⟨L²⟩, ⟨L̂z²⟩, the commutator norms: our analysis of the run (payload);
- Span (Lowdin) = the converged window's S-metric shell weights
  (fractions of the window's columns per l -- the span, not the
  electron configuration; window identity is what the verdicts name);
- Elec (Lowdin) = the ground root's electron count per l
  (analysis.shell_occupations) -- partitions the active electron
  count; blank in payloads that predate the field;
- the fold-shift bookkeeping (|E(analysis fold) − E(as run)|) is NOT a
  column: its reading lives in README · Methods (Partition
  bookkeeping);
- verdict: the HUMAN step (AUDIT below) -- basin verdicts branch on
  window weights (facts 15g/20), never on noise. A cell renders filled
  only when BOTH payload and verdict exist; otherwise pending.

Applicability: the Ti-style factor-grid tables (WINDOWS below). The
study inserts stay with atomterm.study.render; the per-atom summary
tables (C/Si/Sc) stay hand-curated until their manifests migrate.

TABLES edit protocol lives here: sections are rebuilt from constructed
blocks, structural assertions checked (each heading once, row counts,
no starred numbers), temp file + verified move.

Usage: python scripts/build_methods_tables.py [Ti]
"""

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, "src")
from atomterm.cells import CellKey

EV = 27.2114
SECTOR = (3, 1)          # Ti CAS(4,x), Ms=1
DASH = "—"

# the declared windows: (heading, window label, NO label, note lines).
# The NO label is the full-span occupation-group prefix of the window's
# own source density (k = its group count: {3d,4s} -> [s,d] = 2,
# {3d,4s,4p,4d} -> [s,d,d',p] = 4) -- the NOs start at matched span.
# Captions stay minimal (reading conventions live in README · Methods;
# the anchor is the separately-specified ΔE reference).
WINDOWS = [
    ("## CAS(4,6) / {3d,4s}", "3d4s", "3d4sno2", [
        "Reading — README · Methods (⟨L̂z²⟩ is a gauge diagnostic, not a"
        " label). ROHF anchor: −848.4065 Eh (−23086.3167 eV)."]),
    ("### CAS(4,10) / {3d,4d} — planned (constrained)", "3d4d", "3d4dno2", [
        "No 4s (constrained partition). Pending = declared, not yet run; "
        "the CASCI anchor lives in the 3d²4s² lattice above."]),
    ("### CAS(4,14) / {3d,4s,4p,4d} — planned", "3d4s4p4d",
     "3d4s4p4dno4", [
        "The complete valence window. Pending = declared, not yet run; "
        "the CASCI anchor lives in the 3d²4s² lattice above."]),
]

# the audit verdicts (the human step): cell address -> verdict. Basins
# branch on window weights: this draw's free/virtmix CASSCF escaped
# (s 1/6, p 1/2, d 1/3, f-bearing), the free SAs swapped 4s<->3p
# (s -> 4e-6, p -> 1/6) -- facts 15g/20/21c. "MO start was not closed" = the
# construction's span never rotation-closes (the canonical MO slice,
# fact 22d) -- not an optimizer outcome.
AUDIT = {
    "ti_3d4s_3-1_casci_ao_-_pol_raw": "baseline (detector)",
    "ti_3d4s_3-1_casci_ao_-_pol_avg": "exact",
    "ti_3d4s_3-1_casci_ao_-_pur_avg": "exact (pur)",
    "ti_3d4s_3-1_casci_ao_-_pol_frac": "different functional",
    "ti_3d4s_3-1_casscf_ao_-_pol_avg": "escape",
    "ti_3d4s_3-1_casscf_mo_-_pol_avg": "escape (same minimum)",
    "ti_3d4s_3-1_casscf_ao_virtmix_pol_avg": "escape (diagnostic)",
    "ti_3d4s_3-1_casscf_ao_radial_pol_avg": "gauge",
    "ti_3d4s_3-1_casscf_ao_radial_pur_avg": "gauge",
    "ti_3d4s_3-1_casscf_ao_frozen_pol_avg": "baseline",
    "ti_3d4s_3-1_casscf_ao_frozen_pur_avg": "exact (pur)",
    "ti_3d4s_3-1_sa_ao_-_pol_avg": "swap",
    "ti_3d4s_3-1_sa_ao_-_pur_avg": "swap (pur)",
    "ti_3d4s_3-1_sa_mo_-_pol_avg": "escape",
    "ti_3d4s_3-1_sa_ao_radial_pol_avg": "floor",
    "ti_3d4s_3-1_sa_ao_radial_pur_avg": "floor",
    "ti_3d4s_3-1_sa_ao_frozen_pur_avg": "exact",
    # ---- the grid completion (this draw) ----
    "ti_3d4s_3-1_casscf_ao_-_pur_avg": "escape",
    "ti_3d4s_3-1_casscf_mo_-_pur_avg": "escape (same minimum)",
    "ti_3d4sno2_3-1_casscf_no_-_pol_avg": "escape (same minimum)",
    "ti_3d4sno2_3-1_casscf_no_-_pur_avg": "escape (same minimum)",
    "ti_3d4s_3-1_casscf_mo_frozen_pol_avg": "MO start was not closed",
    "ti_3d4s_3-1_casscf_mo_frozen_pur_avg": "MO start was not closed",
    "ti_3d4sno2_3-1_casscf_no_frozen_pol_avg": "baseline",
    "ti_3d4sno2_3-1_casscf_no_frozen_pur_avg": "exact (pur)",
    "ti_3d4s_3-1_sa_mo_-_pur_avg": "escape",
    "ti_3d4sno2_3-1_sa_no_-_pol_avg": "swap",
    "ti_3d4sno2_3-1_sa_no_-_pur_avg": "swap (pur)",
    "ti_3d4s_3-1_sa_ao_frozen_pol_avg": "baseline",
    "ti_3d4s_3-1_sa_mo_frozen_pol_avg": "MO start was not closed",
    "ti_3d4s_3-1_sa_mo_frozen_pur_avg": "MO start was not closed",
    "ti_3d4sno2_3-1_sa_no_frozen_pol_avg": "baseline",
    "ti_3d4sno2_3-1_sa_no_frozen_pur_avg": "exact",
}

HEAD = ("| method | start | mask | substrate | fold | E, Eh | ΔE (mEh) | "
        "E, eV | ΔE (eV) | ⟨L²⟩ | ⟨L̂z²⟩ | ‖[Ĥ,L̂z]‖ | ‖[Ĥ,L²]‖ | "
        "Span (Lowdin) | Elec (Lowdin) | verdict |")
SEP = "|" + "---|" * 16

# CASCI fold variants: (substrate axis, fold axis, display fold)
CASCI_FOLDS = [("pol", "raw", "raw"), ("pol", "avg", "averaged"),
               ("pur", "avg", "averaged"), ("pol", "frac", "fractional")]


def key(atom, window, method, start, mask, substrate, fold):
    return CellKey(atom, window, SECTOR, method=method, start=start,
                   mask=mask, substrate=substrate, fold=fold)


def _cfg(w: dict) -> str:
    """The converged window's SPAN composition (S-metric Löwdin
    weights, the payload ``window``) -- fractions of the window's
    columns per l, summing to 1 (NOT the electron configuration):
    the span is what diagnoses window identity (the swap reads
    s -> 0). Compact: s0 p.17 d.83 f4e-4."""
    out = []
    for l in ("s", "p", "d", "f"):
        v = float(w.get(l, 0.0))
        if v < 5e-5:
            out.append(f"{l}0")
        elif v < 0.001:
            out.append(f"{l}{v:.0e}")
        else:
            out.append(f"{l}{v:.2f}")
    return " ".join(out)


def _elec(n: dict | None) -> str:
    """The ground root's electron count per l (payload ``electrons``,
    analysis.shell_occupations) -- partitions the active electron
    count. Blank when the payload predates the field. Compact:
    s1.2 p0 d2.8 f0."""
    if n is None:
        return ""
    out = []
    for l in ("s", "p", "d", "f"):
        v = float(n.get(l, 0.0))
        out.append(f"{l}{v:.1f}" if v >= 5e-3 else f"{l}0")
    return " ".join(out)


def _values(p: dict, anchor: float) -> list[str]:
    asrun = float(p["energies_pyscf"][0])
    de_eh = asrun - anchor
    out = [f"{asrun:.4f}", f"{de_eh * 1e3:+.1f}", f"{asrun * EV:.4f}",
           f"{de_eh * EV:+.4f}", f"{p['ground_l2']:.4f}",
           f"{p['ground_lz2']:.4f}", f"{p['comm_norm_h_lz']:.1e}",
           f"{p['comm_norm_h_l2']:.1e}", _cfg(p["window"]),
           _elec(p.get("electrons"))]
    return [re.sub(r"(?<!e)-", "−", v) for v in out]   # sign minus only


def _row(cells, audit, addr, method, start, mask, substrate, fold, anchor):
    verdict = audit.get(addr)
    p = cells.get(addr)
    head = [method, start, mask, substrate, fold]
    if p is None or verdict is None:
        return "| " + " | ".join(head + ["pending"] + [""] * 9
                                 + ["pending"]) + " |"
    return "| " + " | ".join(head + _values(p, anchor) + [verdict]) + " |"


def render_window(atom, window, no_label, cells, audit, anchor) -> list[str]:
    out = []
    # CASCI block: the fold IS the row's treatment -- substrate and fold
    # explicit on EVERY row (blank inheritance would misread row 2 as
    # raw); the mask dash repeats (the CASCI non-axis)
    for i, (sub, fold, fold_name) in enumerate(CASCI_FOLDS):
        addr = key(atom, window, "casci", "ao", None, sub, fold).address
        out.append(_row(cells, audit, addr,
                        "**CASCI**" if i == 0 else "",
                        f"AO {{{window}}}" if i == 0 else "",
                        DASH,
                        "polarized" if sub == "pol" else "purified",
                        fold_name, anchor))
    # optimizer blocks: method x mask x substrate x start. Factor cells
    # are explicit on each group head (mask x substrate, the AOs row)
    # and blank-inherited after. RADIAL rows exist only for start=ao
    # (the mask rotates member-labeled AO slots; MO/NO columns carry
    # none) and NOs rows appear everywhere else -- each declared cell
    # is executable through the seam.
    for method, m in (("CASSCF", "casscf"), ("SA-CASSCF", "sa")):
        first = True
        for mask, mask_name in ((None, "free"), ("radial", "radial"),
                                ("frozen", "frozen window")):
            for sub in ("pol", "pur"):
                sub_name = "polarized" if sub == "pol" else "purified"
                starts = ("ao", "mo", "no") if mask != "radial" else ("ao",)
                for st in starts:
                    st_name = {"ao": f"AO {{{window}}}",
                               "mo": "PySCF MO",
                               "no": "NO"}[st]
                    # the NOs row addresses the full-span group prefix
                    # of the source (k = the window's group count)
                    w = no_label if st == "no" else window
                    addr = key(atom, w, m, st, mask, sub, "avg").address
                    head_row = st == "ao"
                    out.append(_row(
                        cells, audit, addr,
                        f"**{method}**" if first else "",
                        st_name,
                        mask_name if head_row else "",
                        sub_name if head_row else "",
                        "averaged" if head_row else "", anchor))
                    first = False
                if method == "CASSCF" and mask is None and sub == "pol":
                    addr = key(atom, window, "casscf", "ao", "virtmix",
                               sub, "avg").address
                    out.append(_row(cells, audit, addr, "",
                                    f"AO {{{window}}}", "virtual-mix",
                                    "", "", anchor))
    return out


def main(atom: str) -> None:
    src = Path(f"elements/{atom}/TABLES.md")
    data = json.loads(Path(f"elements/{atom}/cells.json").read_text())
    cells, audit, anchor = data["cells"], AUDIT, data["rohf_energy"]

    sections = []
    for heading, window, no_label, notes in WINDOWS:
        sec = [heading, ""] + notes + ["", HEAD, SEP]
        sec += render_window(atom, window, no_label, cells, audit, anchor)
        if heading.startswith("## "):
            sec += ["", "Source: `cells.json` (this draw; verdicts — the "
                    "AUDIT table of `scripts/build_methods_tables.py`)."]
        sections.append("\n".join(sec))

    text = src.read_text()
    b0, b1 = "## CAS(4,6)", "## Open questions"
    assert text.count(b0) == 1 and text.count(b1) == 1
    new = text[:text.index(b0)] + "\n\n\n".join(sections) + "\n\n" \
        + text[text.index(b1):]

    # structural assertions (the TABLES protocol)
    assert new.count("## CAS(4,6)") == 1
    assert new.count("### CAS(4,10)") == 1 and new.count("### CAS(4,14)") == 1
    span = new[new.index("## CAS(4,6)"):new.index("## Open questions")]
    assert "C E, Eh" not in span and "S E, Eh" not in span
    assert not re.search(r"\d\*", span)          # no starred gauge numbers
    # row count per window: 4 CASCI folds + per method (free 2 substrates
    # x 3 starts + radial 2 x 1 (ao only) + frozen 2 x 3) + the virtmix
    # diagnostic once (CASSCF free pol)
    per_method = 2 * 3 + 2 * 1 + 2 * 3
    expect = len(CASCI_FOLDS) + 2 * per_method + 1
    for chunk in span.split("\n### ")[1:]:
        chunk = chunk.split("\n#")[0]
        rows = [r for r in chunk.splitlines() if r.startswith("|")]
        assert len(rows) - 2 == expect, len(rows)
    head = span.split("\n### ")[0]
    rows = [r for r in head.splitlines() if r.startswith("|")]
    assert len(rows) - 2 == expect, len(rows)
    tmp = src.with_suffix(".md.tmp")
    tmp.write_text(new)
    assert tmp.read_text() == new
    tmp.replace(src)
    print(f"{src}: methods tables rendered from cells.json "
          f"({sum(1 for a in audit if a in cells)}/{len(audit)} "
          "manifest cells audited)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Ti")
