# Si — silicon (p²)

cc-pVDZ; reference **ROHF −288.846437 Eh**. Standard sectors Ms = 1 /
Ms = 0 — "as carbon" throughout, with its own numbers. Ticket 09's
patient: the multi-solution well was a construction bug (non-orthonormal
build_full_mo), not physics.

Provenance: `results/*_ladder.md`, `tests/test_ladder.py` (Si block),
`tests/test_divergence_guard_fires_on_broken_orthonormality` (the
ticket-09 guard's negative control). Cross-atom verdicts: README
§7.3–7.5.

## Term content

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,3) / {3p} | Ms=1 / Ms=0 | 3 / 9 | as carbon (³P, ¹D, ¹S) |

## Treatment ladder (ΔE vs ROHF, mEh)

| rung | ΔE | note |
|---|---|---|
| CASCI, averaged fold | +221.5 | exact — the confinement cost |
| CASCI, raw core fold | +221.5 | comm 5.6e-4 |
| CASCI, frac reference | +217.8 | |
| CASSCF, free | +0.1 | lands symmetric |
| CASSCF, window frozen | +217.5 | |
| CASSCF, radial restriction | +0.1 | exact: the whole gain is the θ\* rotation |
| SA-CASSCF, free | +0.04 | lands symmetric |
| SA-CASSCF, frozen, purified core | +217.7 | |
| SA-CASSCF, radial, purified core | +0.04 | exact |

## Scan anatomy (§7.4)

One m-symmetric radial rotation θ\* = 0.85 over 3p → 4p (−126.7 mEh) at
exact symmetry; single-member tilts exactly flat. Ticket-16a: the
ENTIRE free-CASSCF gain is radially reachable, radial rung == free rung
to ≤3e-7 Eh, with exact shell purity — which gauge a single-root
optimizer lands in is chaotic (Si drew 1.2e-1 member asymmetry where C
drew 5e-11).

## Window-definition note

max|dE| for CASCI rungs ≈ 9.5e-2 Eh vs the raw-window definition is the
raw-vs-orthogonalized window offset (ticket 14/09), not fold anisotropy.
The fact-17 baseline: the CASSCF input-window orthogonalization defect
~1e-3 on Si (the pur substrate removes it).
