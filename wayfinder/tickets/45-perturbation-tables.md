# Ticket 45 — The PT2 runs + PERTURBATION_TABLE.md

## Goal

The step-by-step PT2 application on the stored C/Si solutions, then
the render. Steps 3–4 of the spec's ladder: (3) one C SA cell —
per-root rows + the two probes (gauge redraw, component identity; the
symmetry answer measured per README · Analysis); (3b) the
reference-class axis — NEVPT2 on stored single-root CASSCF/CASCI
cells vs SA, read against symmetry amplitude; (4) the C then Si
window-ladder read + the render.

## Status: 3b MEASURED 2026-11-19 (the reference-class axis, C)

| reference (compact span unless noted) | ⟨L²⟩ | E(2) | PT2 total | probes |
|---|---|---|---|---|
| SA (primary) | 2.000000 | −58.1340 mEh (−1.5819 eV) | −37.74052528 | clean |
| CASCI pure-AO span | 2.000000 | **−265.7585 mEh (−7.2317 eV)** | −37.69736624 | clean — 2nd L=1 block, different basis |
| single-root @r1 | 2.000000 | −58.0119 mEh (−1.5786 eV) | −37.74042974 | clean |
| single-root mo-start | 1.917610 | −58.0119 mEh (−1.5786 eV) | −37.74042974 | E(2) gauge-clean; ⟨L²⟩ control fails = the amplitude |
| single-root compact ao | 1.683618 | **REFUSED** (stored-roots guard) | — | declared void |

READINGS: (1) PT2 washes out the rotational defect at this scale —
the broken mo-start reference prices the SAME correction as its kept
twin (E(2) is orbital-sensitive, rotation-insensitive). (2) The
window axis dwarfs the reference axis: 58.0 → 265.8 mEh across
windows, ≤ 0.12 mEh across reference classes on one span. (3) PT2
does NOT repair the thin window: its corrected total stays 43.2 mEh
above the SA-corrected total — PT2 closes 207.6 of the 250.8 mEh
poverty-class gap, the remaining 43.2 mEh is genuinely static (the
radial correlation belongs in the window). (4) The refusal: the
16-cycle near-converged draw's stored CI/orbital pair is not
self-consistent on today's substrate (max_de_vs_pyscf 5.3e-5) —
recovery = a regen (user-gated), the row stays void until then.
Note: single-root cells carry a two-ruler 53 µEh split (kernel
stored E = −37.68241788 vs analysis e_ground −37.68236463) — rows
read e_ref from KERNEL totals only (the standing discipline).

## Status: BLOCKED BY 44 (superseded by the 3b record above; step 4 next)

## Inputs

- `elements/<Atom>/perturbation.json` (the payload the adapter
  writes; rows per cell × root × variant, per-class pairs, meters,
  provenance);
- `scripts/build_perturbation_tables.py` (the renderer: ONE table +
  CSV mirror between markers; refuses stale provenance — the
  stale-render precedent);
- `scripts/validate_pt2.py` (the gate: anchors, identities, probes;
  green output committed to `logs/`).

## Verification

- The table renders ONLY from the payload; byte-determinism pins;
  declared-roots policy (no silent draws).
- The reading key (README · Analysis + the element files): R over the
  ladder, the gauge-spread disclosure, the reference-class axis, the
  basis-thinness caveat, the per-L expectations.
- The gate's green log committed with the render (provenance).
- Fast tier green; counts lockstep.
