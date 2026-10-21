# Ticket 22: window methodologies — from AOs to NOs

Status: OPEN (write-only; execute after ticket 21 part 2, before/with 23)

Origin: user direction 2026-10-21. Frame: the standard field approach
seeds active spaces from SCF solutions; our U campaign used AO
primitives; BOTH are legal under our closure theorem, and the choice
should be an experiment axis, not an accident of history.

## The theorem that licenses all of it

A window is exactly rotation-closed iff per l it is (a radial subspace
E_l) (x) (all 2l+1 members) — for ANY E_l (fact 20). So the hard
mechanical constraints are ONLY: (i) orthonormality of the window
columns, (ii) m-complete replication over ALL members (the fact-19
AVAS lesson: per-member thresholds produce rotation-defective spans —
0.9855 overlap yet O(0.5) closure defects, measured). Everything else
is a CHOICE, to be tried and evaluated.

## The methodologies (compared on one axis)

- **A. AO primitive bands/ranks** — the U scan protocol: per-primitive
  open-MO weights measured from the gated reference; bands are
  contiguous weight-sets (never hand-picked exponent intervals);
  conclusions require neighboring-band stability (one more / one fewer
  primitive per side -> same ordering).
- **B. Replicated SCF radials** — the ROHF valence radial(s) AND
  correlating virtuals: extract member-resolved radial blocks from the
  MOs (the purified_core_orbitals construction generalized from core
  to valence/virtuals), replicate each radial over all members.
  Occupied + virtual both admissible; including doubly-occupied
  orbitals (7s active) or semicore (6p) varies the CORE/WINDOW split —
  a parallel choice axis, there is no single correct answer.
- **C. NO-occupation selection** — natural orbitals with significant
  fractional occupations (the 0 < occ < 2 heuristic): NEW small helper
  (1-RDM from our CI vectors in the det basis -> diagonalize per spin
  block -> NOs); window = top-occupation NO radials replicated.
- **D. Hybrid** — scan in primitives (unbiased convergence), then
  compress the converged window to its NO basis (same span quality at
  a fraction of the determinants: 4f-quality answers at 2f-sized cost).

## Evaluation criteria (identical for every methodology)

Exact <L^2>/degeneracy (machine class); low-spectrum ORDERING
stability (within and across neighboring windows); occupancy content
(<N> per block, parity labels); reported window identity (span overlap
vs candidates). Amendment tools applicable to any candidate: Loewdin,
replication, pur substrate, SA restoration.

## Scope / acceptance

- Ticket 24 runs FIRST the matched-cost A-vs-B comparison and the
  README §7.x lattice section (user frame, 2026-10-21) — its results
  are this ticket's methodology-B column; build on them, don't redo.
- The selector registry: one signature per methodology; candidates
  combinable with any fold/method treatment (the axes product).
- U as the testbed: run B/C/D against the A-series results (the
  1f/2f/3f/4f tables); does the NO window reproduce the 5L° ground
  and the f3d1/f2d2 ordering at smaller determinant counts?
- Tests: the closure certificate (window_closure machine class) for
  every methodology's output; NO helper pinned against CI occupancies.
- Record: the comparison table (methodology x dets x ordering x
  exactness) in results/u_14e; §3 entries for any new friction.

## Out of scope

Sparse O.v mode (7s-active 6e windows stay deferred); the ergonomics
runner (ticket 20 consumes this registry); NIST comparison.
