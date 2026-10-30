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

Applicability: the Ti-style factor-grid tables (WINDOWS below; the
window set and every cell address come from the declaration --
``cases.CAMPAIGN`` via ``cells.grid_keys``/``cells.planned``; this file
owns presentation only: headings, notes, value formatting). The study
inserts are frozen history (the study module retired 2026-10-29 --
its generated sections stay as committed data); the per-atom summary tables
(C/Si/Sc) stay hand-curated until their manifests migrate.

THE BRANCH RULE (ticket 29 T1b): rows branch on key IDENTITY only
(method/start/mask/substrate/fold); order and grouping live once in the
declaration builder (cells.window_grid) -- this renderer receives
ready-made key lists and formats.

TABLES edit protocol lives here: sections are rebuilt from constructed
blocks, structural assertions checked (each heading once, row counts,
no starred numbers), temp file + verified move.

Usage: python scripts/build_methods_tables.py [Ti]
"""

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from atomterm.cases import CAMPAIGN
from atomterm.cells import grid_keys, window_grid

EV = 27.2114
DASH = "—"

# the declared windows this file renders, with their presentation
# (headings + note lines). The window list MUST equal the declaration's
# (primary + planned) -- asserted in main; the NO full-span prefix and
# every cell address derive from the declaration, not from here.
# Captions stay minimal (reading conventions live in README · Methods;
# the anchor is the separately-specified ΔE reference).
WINDOWS = [
    ("## CAS(4,6) / {3d,4s}", "3d4s", [
        "Reading — README · Methods (⟨L̂z²⟩ is a gauge diagnostic, not a"
        " label). ROHF anchor: −848.4065 Eh (−23086.3167 eV)."]),
    ("### CAS(4,10) / {3d,4d} — planned (constrained)", "3d4d", [
        "No 4s (constrained partition). Pending = declared, not yet run; "
        "the CASCI anchor lives in the 3d²4s² lattice above."]),
    ("### CAS(4,14) / {3d,4s,4p,4d} — planned", "3d4s4p4d", [
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
    "ti_3d4s_3-1_sa-casscf_ao_-_pol_avg": "swap",
    "ti_3d4s_3-1_sa-casscf_ao_-_pur_avg": "swap (pur)",
    "ti_3d4s_3-1_sa-casscf_mo_-_pol_avg": "escape",
    "ti_3d4s_3-1_sa-casscf_ao_radial_pol_avg": "floor",
    "ti_3d4s_3-1_sa-casscf_ao_radial_pur_avg": "floor",
    "ti_3d4s_3-1_sa-casscf_ao_frozen_pur_avg": "exact",
    # ---- the grid completion (this draw) ----
    "ti_3d4s_3-1_casscf_ao_-_pur_avg": "escape",
    "ti_3d4sno2_3-1_casscf_no_-_pol_avg": "escape (same minimum)",
    "ti_3d4sno2_3-1_casscf_no_-_pur_avg": "escape (same minimum)",
    "ti_3d4s_3-1_casscf_mo_frozen_pol_avg": "MO start was not closed",
    "ti_3d4sno2_3-1_casscf_no_frozen_pol_avg": "baseline",
    "ti_3d4sno2_3-1_casscf_no_frozen_pur_avg": "exact (pur)",
    "ti_3d4sno2_3-1_sa-casscf_no_-_pol_avg": "swap",
    "ti_3d4sno2_3-1_sa-casscf_no_-_pur_avg": "swap (pur)",
    "ti_3d4s_3-1_sa-casscf_ao_frozen_pol_avg": "baseline",
    "ti_3d4s_3-1_sa-casscf_mo_frozen_pol_avg": "MO start was not closed",
    "ti_3d4sno2_3-1_sa-casscf_no_frozen_pol_avg": "baseline",
    "ti_3d4sno2_3-1_sa-casscf_no_frozen_pur_avg": "exact",
}

HEAD = ("| method | start | mask | substrate | fold | E, Eh | ΔE (mEh) | "
        "E, eV | ΔE (eV) | ⟨L²⟩ | ⟨L̂z²⟩ | ‖[Ĥ,L̂z]‖ | ‖[Ĥ,L²]‖ | "
        "Span (Lowdin) | Elec (Lowdin) | verdict |")
SEP = "|" + "---|" * 16

# identity -> display (the ONLY branching the renderer does)
METHOD_NAMES = {"casci": "**CASCI**", "casscf": "**CASSCF**",
                "sa-casscf": "**SA-CASSCF**"}
START_NAMES = {"mo": "PySCF MO", "no": "NO"}
MASK_NAMES = {None: "free", "radial": "radial", "frozen": "frozen window",
              "virtmix": "virtual-mix"}
SUBSTRATE_NAMES = {"pol": "polarized", "pur": "purified"}
FOLD_NAMES = {"raw": "raw", "avg": "averaged", "frac": "fractional"}


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


def _row(cells, audit, k, prev, anchor):
    """One table row from a declared key; every presentation choice
    (bold method heads, blank factor inheritance, the CASCI dash)
    derives from the key's identity and its predecessor's."""
    verdict = audit.get(k.address)
    p = cells.get(k.address)
    m_new = prev is None or prev.method != k.method
    # the head-row rule: factor cells fill on AO-start rows (and every
    # CASCI row names its substrate/fold -- blank inheritance would
    # misread row 2 as raw); the virtual-mix diagnostic names only its
    # mask; CASCI repeats the mask dash (no optimizer)
    head = [
        METHOD_NAMES[k.method] if m_new else "",
        "" if (k.method == "casci" and not m_new) else
        (f"AO {{{''.join(k.shell_patterns)}}}" if k.start == "ao"
         else START_NAMES[k.start]),
        DASH if k.method == "casci" else
        (MASK_NAMES[k.mask] if k.start == "ao" else ""),
        SUBSTRATE_NAMES[k.substrate]
        if k.method == "casci" or (k.start == "ao" and k.mask != "virtmix")
        else "",
        FOLD_NAMES[k.fold]
        if k.method == "casci" or (k.start == "ao" and k.mask != "virtmix")
        else "",
    ]
    if p is None or verdict is None:
        return "| " + " | ".join(head + ["pending"] + [""] * 9
                                 + ["pending"]) + " |"
    return "| " + " | ".join(head + _values(p, anchor) + [verdict]) + " |"


def render_window(keys, cells, audit, anchor) -> list[str]:
    """Format one window's ready-made key list (the declaration's order
    -- cells.window_grid); the renderer never constructs addresses."""
    return [_row(cells, audit, k, prev, anchor)
            for prev, k in zip([None] + keys[:-1], keys)]


def main(atom: str) -> None:
    src = Path(f"elements/{atom}/TABLES.md")
    data = json.loads(Path(f"elements/{atom}/cells.json").read_text())
    cells, audit, anchor = data["cells"], AUDIT, data["rohf_energy"]
    spec = CAMPAIGN[atom]

    # the declaration's per-window grids (primary + planned). Render
    # contract: this file's WINDOWS cover EXACTLY the declared windows,
    # primary first -- presentation here, content from the declaration.
    primary = "".join(spec.shell_patterns)
    per_window = {primary: grid_keys(spec)}
    per_window.update({w: window_grid(spec, w)
                       for w in spec.planned_windows})
    assert [w for _, w, _ in WINDOWS] == [primary, *spec.planned_windows], \
        "the renderer's WINDOWS and cases.CAMPAIGN diverged"

    sections = []
    for heading, window, notes in WINDOWS:
        sec = [heading, ""] + notes + ["", HEAD, SEP]
        sec += render_window(per_window[window], cells, audit, anchor)
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
    # row count per window: 4 CASCI folds + per method (free: 3 starts x
    # 2 substrates; radial: ao only x 2 substrates; frozen: 3 starts x 2
    # substrates -- mo x pur IS in the grid since the projection
    # builder) + the virtmix diagnostic once (CASSCF free pol)
    per_method = 3 + 3 + 2 + 6
    expect = 4 + 2 * per_method + 1
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
