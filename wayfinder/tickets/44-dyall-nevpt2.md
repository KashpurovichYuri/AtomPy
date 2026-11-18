# Ticket 44 — Dyall H₀ + NEVPT2 strongly contracted (rungs 2–3)

Spec: `wayfinder/specs/pt2-perturbation.md`.

## Goal

The parameter-free PT2: the Dyall diagonal (the three-way
inactive/active/virtual range separation; canonical orbital energies
in all ranges; the empty-active limit is exact MP2), the active RDM
stack to the 4-cumulant, the 7 contracted classes
(`[1V][2V][3V][4V][V][W][X]`) as small per-class solves with
denominator-floor meters — intruder-freeness VERIFIED (per-class
floors ≥ 0 on every case), never assumed. The validation script
grows the new layers' certificates.

## Status: BLOCKED BY 42, 43

## Inputs

- rung 0–1 machinery (the FOIS enumeration, the solve grammar, the
  meters);
- `natural.py`'s lazy-driver pattern for the RDM builders' imports
  (the analysis-row import discipline);
- the density-cumulant hierarchy (the connected parts of the RDMs —
  γ, Δ, up to the 4-cumulant).

## Verification

(c) the empty-active → MP2 identity (the textbook known answer —
any deviation is a bug, not physics); (a) cumulant reconstruction
identities (the 4-cumulant reproduces the lower RDMs elementwise);
per-class floor ≥ 0 signs on every case; class dimensions vs the
FOIS class counts; the carbon fixture pins the first NEVPT2 number
(observed fresh). CASPT2 (rung 4) stays specified-not-built — its
level shift masks exactly the intruder signal the floors read.
