# Ticket 46 — The TABLE.md physics pass: the significance policy

## Goal

A reviewed, pruned TABLE.md render — RENDER-LEVEL ONLY (the
`cells.json` stores keep every payload; pruning is the row policy,
never the record). DESIGN-FIRST: the design session proposes the
significance criteria (candidates: the m-incomplete MO/MOC controls'
render fate, the duplicate raw/avg ruler rows, the no2 variants,
per-family row multiplicity), produces the per-atom exclusion census
UNDER each criterion, and the user adjudicates the exclusions in
discussion before any re-render. Nothing disappears from the record —
the store keeps it; the table keeps the physics.

## Status: BLOCKED BY nothing (design-first, user-gated per atom;
default slot: after ticket 45 — one combined tables discussion)

## Inputs

- `scripts/build_methods_tables.py` (the row policy's one home);
- the per-atom `TABLE.md`/`cells.json` pairs (the census substrate).

## Verification

The re-render loses no payload; the exclusion census is committed
with the re-render (the record of WHAT was excluded and WHY, per
criterion); the fast tier stays green.
