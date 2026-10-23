# Ticket 23: U SA wave — orbital relaxation vs window growth

Status: EXECUTING (2026-10-24; `scripts/run_u_sa.py`, staged:
prep ✅ NO cache [groups 5/7/5/7, anchors reproduce the pilot log] →
casscf ladder [NO nk=1..3, shell-12 vs NO-12 matched pair, shell-24 vs
NO-24 same-span pair] → sa nk=2,3 → relax nk=1..3 → six13)

Origin: the settled U plan (2026-10-21). The campaign's open physics
question, and one of ticket 22's comparisons in disguise.

## The question

The window-convergence series (1f: f1d3 -> 2f: f2d2 -> 3f: 5L° odd
ground) established that RADIAL CONTENT decides the configuration
ordering. SA-CASSCF optimizes window orbitals freely — can relaxation
recover the missing radial character INSIDE a smaller window and reach
the converged ordering without window growth? If yes, the expensive
40k-detector sectors become unnecessary for identification; if no,
window growth is irreducible and the scan protocol is mandatory.

## The runs (maximal Ms=2 sector, (4,0))

0. **Free single-root CASSCF FIRST** (user direction 2026-10-21, before
   any SA): from the shell start in the 3f/4f windows. This is where the
   "many different minima characterized by a specific term" phenomenology
   lives (the fact-8/15/16 class: symmetry-broken basins, member-gauge
   chaos, escape/swap windows). Measure: which term the optimizer lands
   in, ⟨L²⟩ breaking, window weights, ortho_err; branch on basins.
1. **17-root SA over the 5L° M_L set** (L=8: M_L = -8..+8; weights
   1/17 each) in the 2f and 3f windows. Read: does the ensemble window
   borrow the missing f-radial character (occupancies/weights per
   rung); does the per-root table (our CASCI pass in the SA window)
   reproduce the 3f ordering from the 2f window?
2. **relax / pur-relax rungs** (frozen window, core relaxes; pur
   substrate): TM precedent says machine-exact closure expected —
   these are the honest symmetric rungs for U.
3. Free CASSCF from OTHER starts (canonical, NO-like) = the basin map,
   only after 0-2 established the reference phenomenology.

## Acceptance

- The relaxation-vs-growth verdict, measured (ordering comparison
  across windows x {frozen CI, SA, relaxed}).
- Per-rung fingerprints in the ticket + results/u_14e; chaos branches
  documented; §3 entries for friction.
- Tests: only deterministic objects (frozen-window rung closure on
  pur; the SA-window CASCI cross-checks) — optimizer reruns are not
  pinned (fact 10 discipline).

## §3 Friction

- 2026-10-24: `report()` derived `ncas = mo.shape[1] - ncore` — the
  fact-7 trap on the OUTPUT side: a free CASSCF returns the FULL
  (172-col) matrix, so the "active" slice grabbed the whole virtual
  tail (C(167,4) = 35.2M dets -> 17.6 PiB alloc -> crash). Same bug
  was waiting in the SA and six13 stages. Rule: `ncas` comes from
  `mc.ncas`, slices are `[ncore : ncore+ncas]` — on outputs as well as
  inputs. (The crashing run itself had converged and persisted; only
  forensics died. Third hiding place found later the same day: the
  PRE-SA `ground_block` call — every consumer of a constructed or
  returned MO matrix needs the bounded slice, not just the obvious
  ones.)
- Script stdout to a file is block-buffered: a "silent" stage may be
  20 cores deep in a macro iteration. Run `python -u` for staged
  drivers.
- PySCF single-root CASSCF returns `mc.ci` as a LIST wrapping the
  vector (`CASResult.ci` becomes doubly wrapped); normalize with
  `np.asarray(ci).ravel()` before any contraction.

## Findings (running)

- 2026-10-24, FREE RUNG IS DEAD ON U (u_sa_casscf.log, u_sa_diag.log):
  free single-root CASSCF from EVERY start (NO nk=1/2/3, shell 12/24,
  NO-24) dives ~540 mEh below the frozen curve — macro step 1 alone
  takes −493 mEh (window identity destroyed in one rotation step) —
  and lands in a symmetry-annihilated flat region: <L^2> = 0.6..31
  vs the 72 anchor, span defects |[P,Lperp]| = 6.5..7.8, ||[H,Lz]||
  up to 270, <S^2> = 6 exact throughout. The diag (100 cycles,
  verbose) shows the endpoint is NOT a certified stationary point:
  |grad_o| stuck at ~1e-2, energy crawling −18 uEh/cycle (cycle 50 ->
  100 gained only 0.44 mEh) — basin-like, never converged. The NO
  start does NOT tame the escape; at the same span (24 cols) shell vs
  NO land 0.14 mEh apart with chaotic gauges (<L^2> 30.9 vs 23.3).
  5S° ATTRACTOR (ticket-24 connection): the nk=1 free run lands AT
  the ROHF reference energy (−50.378690 vs −50.381852, +3.2 mEh) with
  <L^2> = 0.64; the diag's cycle 1 passes through −50.38279 en
  route. The 5S° reference is an attractor of the U correlation
  landscape, not an SCF accident (AOC-smearing mechanism still open).
  Honest rungs: SA + relax (running).
- 2026-10-24, SA DOES NOT HOLD 5L ON U EITHER (u_sa12.log, u_sa17.log):
  SA over the measured 17-root 5L block (weights 1/17) CONVERGES
  (certified, unlike free) but to broken windows. sa/no-12: ensemble
  −50.423810, span cert 5e-4/8e-4 (internally nearly closed) yet
  per-root <L^2> = 24.7..30.7 (not 72) and the window identity traded
  ~3.5 columns away (f 7.000 -> 4.706, d 5.000 -> 3.764): the whole
  block collapsed into ONE broken well (roots within 10 uEh). sa/no-17:
  span cert 4.3/5.4 (open, free-class), <L^2> = 26.6..31.2, root spread
  2.6 mEh. VERDICT CLASS: on U the escape attraction (534 mEh, toward
  the 5S°-like region) beats state averaging — qualitatively new vs the
  TM precedent (fact 12/15: SA restored <L^2> exactly). LESSON: the
  span certificate measures closure, not identity — SA rebuilt the
  window from other (m-complete) material at 1e-4 closure while the
  fold went anisotropic and the CASCI roots in it are broken. The
  <L^2> = 27-class with cert 1e-4 is consistent: closed span, broken
  fold, broken eigenstates.
- 2026-10-24, REVISION PENDING: the SA/relax <L^2> labels came from
  PySCF ci vectors -- the Davidson straggle in degenerate clusters (the
  run_u.py disease) can fake broken <L^2> inside a perfectly symmetric
  fold (relax/no-12: span cert 1e-13, ||[H,Lz]|| = 9e-14, yet
  <L^2>(mc.ci) = 27.4 -- impossible for a true eigenstate). Honest
  labels = exact eigh in the converged window; stages now print both.
  Solver-independent so far: the relax ENERGIES (frozen span + relaxed
  core): no-5 -49.9423 (still descending, converged=False, <L^2>
  exact 6 by symmetry of the 5-det d4 space), no-12 -49.8997 (-10.00
  mEh), no-17 -49.9004 (-9.51 mEh); the SA window-identity drift
  (f 7.000 -> 4.706). six13 first pass was BUGGY (7s column placed
  FIRST in cols while the occupation indices assumed LAST; negative
  <N7s> = -0.110 printed -- impossible for a normalized expectation;
  never ship a negative occupation). Rerunning with assertions +
  exact-eigh labels.


- 2026-10-24, RELAX VERDICT FINAL (u_sa_relax.log, exact-eigh labels):
  the Davidson straggle CONFIRMED as the <L^2> liar -- relax/no-12
  mc.ci said <L^2> = 27.4 while the EXACT eigh ground in the same
  window/fold is <L^2> = 72.0000. Honest relax rungs on the NO
  substrate:
      relax/no-12: -10.00 mEh vs frozen-12, <L^2> = 72.0000,
                   span cert 1e-13, ||[H,Lz]|| = 1.3e-13, converged
      relax/no-17: -9.51 mEh, <L^2> = 72.0000, cert 1e-12, converged
      relax/no-5 : EXCLUDED (converged=False; E moved -49.942 ->
                   -49.800 across reruns -- truncated endpoints are
                   chaotic, converged ones reproduced exactly:
                   fact-10 discipline re-confirmed empirically)
  RELAXATION VS GROWTH (the ticket question, honest rungs only):
  core relaxation at frozen window = 10.0 mEh at exact symmetry;
  the ENTIRE frozen growth axis 12->24 cols = 4.56 mEh. Relaxation
  outbuys growth ~2:1 WITHIN the symmetric budget; the remaining
  ~520 mEh exists ONLY behind symmetry breaking (the 5S-deg
  attractor's price -- free/SA territory). Growth is NOT irreducible
  for CORRELATION (relax buys more at exact symmetry), but the frozen
  windows remain the only space where the terms are certified; no
  optimizer rung (free or SA) keeps 5L on U while taking the deep
  money. The SA "broken <L^2>" verdict PENDS the exact-eigh re-read
  (queued): only the window-identity drift stands so far.

- 2026-10-24, six13 "negative occupations" RESOLVED -- user error, not
  physics: the expectation of a diagonal operator is SUM |c_i|^2 n_i,
  not the signed SUM c_i n_i. The |vec|^2 diagnostic gave <N7s> =
  1.9885 (sane: the 7s^2 pair intact in every low root) where the
  signed dot gave -0.41. The unit-vector x non-negative-counts
  "paradox" was the wrong contraction. Rule: |c_i|^2-weighting for
  diagonal pops; sandwich (ci @ M @ ci) for L^2/S^2; the
  <N_total> = nele assert now guards the six13 table permanently.
  Bonus physics already visible: <N7s> ~ 1.99 in ALL low roots -- the
  7s pair is passive (not correlation-active) in the NO substrate.

- 2026-10-24, SA VERDICT REVERSED BY THE EXACT RE-READ (u_sa12.log):
  the SA-12 endpoint's exact eigh ground is <L^2> = 72.0000 EXACT
  with ||[H,Lz]|| (SA fold) = 1.3e-4 and E_tot = -50.4238 -- the
  earlier "<L^2> ~ 27 broken" labels were Davidson straggle on the
  17-fold cluster, not physics. SA-17 stays degraded (fold broken,
  [H,Lz] = 10.2, <L^2> = 71.88). REVISED VERDICT TABLE (exact labels
  only):
      frozen-12      -49.8897  <L2>=72.0000  exact (anchor)
      frozen-24      -49.8943  <L2>=72.0000  exact (growth: +4.6 mEh)
      relax@12       -49.8997  <L2>=72.0000  exact (core relax: +10.0)
      SA@12          -50.4238  <L2>=72.0000  fold 1.3e-4, cert 8e-4
      SA@17          ~-50.428  <L2>=71.88    fold BROKEN (10.2) -- out
      free (any)     ~-50.43   broken        unconverged -- out
  HEADLINE: SA over the 5L block from the NO-12 start = -50.4238 at
  EXACTLY quantized 5L -- 42 mEh BELOW the ROHF 5S-deg mean field
  (-50.3819) and 470 mEh below the best frozen CASCI. The first U
  object to beat the 5S-deg reference at certified symmetry (the
  ticket-24 quest). Relaxation vs growth resolves STRONGLY YES on U --
  via the ENSEMBLE anchor: the deep money (~530 mEh) is real and
  symmetry-accessible, but only SA holds the gauge; free single-root
  is destroyed by the same attraction it chases. The SA window's
  identity drifted (f 4.7/d 3.8) -- l-content readout added to
  stage_sa to identify what entered (suspect: m-complete core sets =
  core-pair correlation).
- 2026-10-24, six13 IDENTIFICATION BUG: pur core columns are
  [s, s, p, p, p]; col 4 is a 6p MEMBER, not the 7s. The first six13
  run promoted a 6p pair to active (mo_energy[:ncore] -- canonical
  order meaningless on this basis, the standing warning!). Its
  numbers are relabeled: the "7s" occupation 1.989 = the promoted
  6p member's; span defect 2.0 = the lone m-incomplete p member;
  <L2> = 72.02 near-quantized because that pair stayed ~closed.
  REAL 7s2 probe rerunning with the l-sector identification
  (s-columns {0,1}; 7s = the one on the highest-energy OCCUPIED
  canonical). Rule: identify orbitals by l-sector content first,
  then by overlap on the OCCUPIED canonical block.

- 2026-10-24, FINAL BLOCK (u_six13.log, u_sa12.log):
  six13 (correct 7s identification): E = -49.894866, <L2> = 72.0000
  EXACT, span cert 2e-13/4e-13 machine-exact, <N7s> = 1.972, <Nf> =
  2.000, <Nd> = 2.028 -- 7s2-as-valence moves NOTHING (-0.58 mEh);
  the ground stays f2d2 5L; all low roots = the degenerate 5L block.
  The buggy run's -551 mEh was the accidental 6p-member promotion:
  INNER core pairs carry huge correlation, the outermost s pair none
  (next tier candidate: core-pair valence done properly, m-complete).
  SA-12 MECHANISM IDENTIFIED: l-content of the converged window is
  exactly l2 + l3 = 12.0 (l0 = l1 = 0.000) -- the window NEVER left
  the f/d shells; the "identity drift" (f-band weight 7.000 -> 4.706)
  is RADIAL REDISTRIBUTION across exponent bands. The missing-radial
  axis of the growth lattice, recovered by relaxation: the ticket-23
  question answered YES in its purest form, via the ensemble anchor.
  Reproducibility: SA-12 ensemble E reproduced EXACTLY (-50.423810)
  across independent reruns (converged SA = deterministic; truncated
  runs chaotic -- fact 10 re-measured).
  Status of the certification: SA-12 is certified by the state
  (<L2> = 72.0000, fold [H,Lz] = 1.3e-4), not by machine-exact span
  closure (8e-4). TODO before closing: sa 3 l-content (degraded?
  why), the SA-12 per-band radial weights (which radials entered --
  forensics for the next ticket), full suite, commit, README done.

- 2026-10-24, RADIAL MECHANISM REFINED (l-content readouts, both SA logs):
  SA-12: total l-content UNCHANGED (l2=5.0, l3=7.0, l0=l1=0.000) but
  band-row weight 8.47/12 -> 3.53 columns of radial weight moved to
  f/d primitives OUTSIDE the selected exponent bands = the neighbor
  radials (ticket-24's growth direction) reached by relaxation, with
  m-completeness kept (cert 8e-4). SA-17: pulled 5 MORE f members in
  (l3 7->12, still shell-pure) but the span is OPEN (4.3/5.4) and the
  fold broken ([H,Lz]=10.2, <L2>=71.88): fact-8's refined statement
  live -- l-purity is necessary, NOT sufficient for closure; 12 f
  members over the f radial chain need not be m-complete. Verdict:
  radial relaxation is the deep axis; the ensemble holds the gauge at
  12 columns and loses it at 17.
