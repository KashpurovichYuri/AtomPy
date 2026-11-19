# Ticket 44 — Dyall H₀ + NEVPT2 strongly contracted (rungs 2–3)

Spec: `wayfinder/specs/pt2-perturbation.md`.

## Status: IN PROGRESS (44a Dyall DONE 2026-11-18; 44b cumulant stack
+ 44c NEVPT2-SC remain)

44a (DONE): `dyall_h0`/`dyall_e2`/`dyall_e2_dense` — the partition
with the folded-core convention stated honestly (the inactive range
IS the folded core, absent from the columns; the textbook
empty-active→MP2 anchor is UNREACHABLE in this convention —
substituted by the closed-form canonical-Fock anchor; the map
records the deviation).  THREE caught bugs, all machine-certified
fixed: (1) the Fock build needs PHYSICAL integrals — the chemist
8-fold (the true detbasis pair swap is transpose(1,0,3,2), NOT
(2,3,0,1) — guarded, refusal test); (2) THE ACTIVE-FIRST RELABEL IS
NOT THE IDENTITY for the beta block (beta q → n_val + q, never
n_act + q) — `_cas_one_rdm` builds γ valence-embedded, `DyallH0.act/
virt` carry the split sets; (3) the closed-form anchor's own
negative-indexing garbage — the test fixed.  The script grows the
Dyall certificates (e0 anchor + dense twin per root).  Fast tier
271; validation C 239 / Si 307 PASS 0 FAIL.  A READING-KEY FINDING
for 45: on `c_2p_..._frozen_pol_avg` root0 the Dyall E(2) reads
−199.4 mEh at w_ref 0.787 (the suspect class; EN read −7.8 mEh at
0.996 on the same root) — the Dyall partition's denominators expose
a near-degenerate external direction the EN exact-diagonal
denominators masked.  Flag, never celebrate.

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
