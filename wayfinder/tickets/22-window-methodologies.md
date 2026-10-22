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

## Findings (2026-10-23, first execution round)

1. NO HELPER LANDED: `src/atomterm/natural.py` (ci_rdm1 /
   natural_orbitals / occupation_groups / select_no_window) + the
   SELECTORS registry (shell / primitive / scf_mo / no, one signature,
   Window(label, ao_rows | c_win, m_complete)). Tests:
   tests/test_natural.py (brute-force 1-RDM agreement, projector
   property, root-linearity, the m-set claim, the fact-1 primitive
   guard, real-AO certificates).
2. THE CLAIM, REFINED AND MEASURED (tests + U pilot): a single
   |L, M_L> root has an m-DIAGONAL but m-UNEQUAL density (p^2 3P
   M_L=0: occ(m=0) = 2, occ(+-1) = 0) -- occupation selection on one
   root is the AVAS mistake. A term-AVERAGED root set (rotation-
   invariant) has an m-SCALAR density: NOs arrive in complete
   l-pure m-sets. Root averaging is what makes NO selection
   symmetry-safe.
3. U PILOT (scripts/run_u_no.py, results/u_14e/u_run_no.log): the
   2f+2d term-averaged density gives exactly 4 groups (d, f, d, f =
   5/7/5/7 columns, occ 0.399/0.286/0.0005/0.0002). Every group
   window: certificate 1e-13/1e-12, <L^2> = 72.0000 exact, <Nf> = 2
   exact, the 17-fold 5L ground intact. Compression: 12 cols (495
   dets) +4.56 mEh vs the full window; 17 cols (2,380) +3.43; AO
   lattice at the same 495 dets: 431 mEh higher. The D-methodology
   promise (full-window quality at a fraction of the cost) CONFIRMED.
   The missing 3.4-4.6 mEh is the second f radial (occ 2e-4).
4. CORRECTION to the campaign labels: the f2d2 ground blocks are
   5L (L = 8, <L^2> = 72) -- the 56.000 labels in the conversation
   tables (and briefly in README 7.5 / the ticket-24 appendix) were
   wrong; campaign logs are authoritative (README fixed this round).
5. E_core bookkeeping: the closed-shell pur-core formula reproduces
   the campaign total to -0.000 mEh (validation built into the pilot).

## Remaining (next round)

- Sc/Ti NO columns of the matched tables (same pilot pattern; the Sc
  1,274-det CASCI is seconds).
- NO windows vs the B (raw-MO) side at matched det counts on Sc/Ti
  (the NO column the tables reserve).
- Methodology B proper (replicated SCF radials via the
  purified_core_orbitals construction generalized to valence) --
  partially covered by ticket 24's MO column; decide after the NO
  columns land whether a separate B-variant row adds anything.
- Tests: the NO helper is pinned; add a NO-window CASCI regression
  (small Sc case) if the NO columns enter the standard tables.
