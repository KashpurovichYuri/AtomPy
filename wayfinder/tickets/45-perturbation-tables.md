# Ticket 45 — The PT2 runs + PERTURBATION_TABLE.md

## Goal

The step-by-step PT2 application on the stored C/Si solutions, then
the render. Steps 3–4 of the spec's ladder: (3) one C SA cell —
per-root rows + the two probes (gauge redraw, component identity; the
symmetry answer measured per README · Analysis); (3b) the
reference-class axis — NEVPT2 on stored single-root CASSCF/CASCI
cells vs SA, read against symmetry amplitude; (4) the C then Si
window-ladder read + the render.

## Status: BLOCKED BY 44

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
