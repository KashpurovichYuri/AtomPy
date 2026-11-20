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
| single-root compact ao | 1.683618 | −58.0119 mEh (−1.5786 eV) | −37.74042975 | clean (recovered — see below) |

READINGS: (1) PT2 washes out the rotational defect AT EVERY MEASURED
AMPLITUDE — all three single-root references (⟨L²⟩ = 2.000 / 1.918 /
1.684, e_ref = ROHF exactly in all) price the SAME correction
−58.0119 mEh to print precision; the SA row differs by 0.12 mEh (its
orbitals differ, e_ref 27 uEh below ROHF).  E(2) is orbital-
sensitive, rotation-insensitive.  (2) The window axis dwarfs the
reference axis: 58.0 → 265.8 mEh across windows, ≤ 0.12 mEh across
reference classes on one span.  (3) PT2 does NOT repair the thin
window: its corrected total stays 43.2 mEh above the SA-corrected
total — PT2 closes 207.6 of the 250.8 mEh poverty-class gap, the
remaining 43.2 mEh is genuinely static (the radial correlation
belongs in the window).  (4) THE RECOVERY (user-approved regen +
the seam's fourth clause): the first refusal was a FALSE positive —
the readout re-solve reproduces the stored state with |overlap| =
1.000000, e_cas equal to 1e-14, E(2) equal to 2.8e-8 mEh, carrying
only a ~1e-4 unconverged-davidson mixing of a quasi-degenerate
partner (the 16-cycle cap-hitting kernel stores its last davidson
state; the re-solve lands a slightly different mixture).  The
band structure is now three clauses with honest powers: e_tol the
energy band, o_tol the re-solve ANGLE band (|sin|, linear in the
mixing, from the orthogonal residual — sqrt(1-cos^2) cancels
catastrophically at small angles), c_tol the exact-amplitude band;
the joint projector stays the multi-root block's invariant.  The
regen (deterministic draw on this box — e_ground reproduced to 15
digits) refreshed the store; the row stands recovered.  Note:
single-root cells carry a two-ruler 53 µEh split (kernel
stored E = −37.68241788 vs analysis e_ground −37.68236463) — rows
read e_ref from KERNEL totals only (the standing discipline).

## Status: step 4 (C) DONE 2026-11-19 -- 26 cells / 66 roots, gate ALL GREEN

The chain: `run_pt2_ladder.py` (the declared set: span ladder
casci x {ao,mo,moc} x {compact,@r1,@s1}; NO twins; the reference axis
{casscf,sa-casscf} x {compact,@r1,@s1} x {ao,mo}; the frozen pair)
-> `elements/C/perturbation.json` -> `build_perturbation_tables.py`
(ONE table + CSV, byte-deterministic, the stamp = the payload digest,
the schema refusal) -> `validate_pt2.py` ALL GREEN
(logs/validate_pt2_C.log): the anchors at 1e-6, the class sums 0e+00
over 66 roots / 528 classes, totals 3.5e-15, eV columns exact, the
meter identity 2.8e-17, 26/26 rows-or-refused, the probe clean.

THE C READING (the ladder read): (1) every OPTIMIZED cell (SA/CASSCF,
all windows and starts) prices -58.11..-58.14 mEh and corrects to
-37.7405 -- the ladder SATURATES at the optimized level; the
optimizer repairs the poverty class, and the windows are equivalent
gauges of one orbital optimum.  (2) The poverty class lives in the
UNOPTIMIZED thin spans: compact-ao CASCI and 2pno1 price -265.76 mEh
(-7.23 eV) -- the fuel reading -- and still land 43.2 mEh short of
the optimized corrected total.  (3) E(2) tracks the orbital DRAW, not
the span's rotational quality: the moc-start rows carry <L2> = 0.082
(yet E(2) = -58.1341), the broken single roots carry 1.68..1.92 (yet
E(2) = -58.0119 across all).  (4) The frozen-core disclosure: the
compact frozen row prices -243.33 mEh -- 4.2x its free twin (the
unrelaxed core's Sir class explodes); the @r1 frozen row reads the
free value.  (5) Broken spans store excited roots beside the ground
one (mo-start roots 1-2: +28.8 mEh, E(2) -75.97) -- the table's
spread column reads the GROUND BLOCK only; full per-root records in
the payload.

DECLINED (the lane axis): the lane payload is its own store class
(per-irrep-block pinned solutions, term energies -- no ci_stack/
core_energy); the stored-roots seam needs a lane-native readout twin
before it can price lane rows.  The irrep route-invariance gate row
waits with it.  Recorded as the follow-up.

## Status: the LANE SEAM + step 5 (Sc, Ti) DONE 2026-11-20; the uncs DONE; U14e in flight

THE LANE SEAM (1edb838): cells.readout_lane -- the runner's own
re-solve recipe asserted against the record (energies strict, the
meters within the re-solve noise band), per-block kernel-ready mc
objects out through perturb.nevpt2_rows_on; casscf.irrep_solvers
promoted (the pin's one construction, runner and twin share it);
cases.lane_substrate moved in from scripts (the substrate family's
second member).  MEASURED (C): the CI lane's ground block prices
-682.3334 mEh = the free twin's root EXACTLY (the irrep pinning is
pricing-inert); the SA lane's B1g prices -58.1340 = the certified
grid SA value.  The lane blocks hold the D2h det-class states (the
meters carry them) -- the lane rows are the core-partition
disclosure axis, as declared.  Route invariance: C delta 0,
Si 1.04e-5 mEh, Sc 5.4e-7 Eh -- pricing-inert where the substrates
agree; Ti DISCLOSES 1.76 mEh (the lane's D2h-canonical draw and the
free rung land the SAME basin and energy with DIFFERENT orbital
gauges -- E(2) moves at the mEh scale: the wave's orbital-draw
finding read at the route axis).

THE SEAM'S BANDS reached their per-atom scale (each measured, none
an occasion-tolerance):
- readout_mc(ci0=): the reconstruction starts FROM the stored stack
  -- a quasi-degenerate manifold's plain re-solve RE-PARTITIONS the
  roots (Sc compact: the same nine energies to 1e-11, a different
  span at 1.7e-1).
- assert_stored_roots(e_noise=1e-6): the energy miss consults the
  vector certificate (the Sc single-root CASSCF class: angle 9e-6,
  energy 2.7e-11..9.2e-9).
- p_tol 1e-4: the tail sweep measured -- Si 2.2e-7, Ti 3.9e-6,
  Ti-unc 5.0e-5; an exactly degenerate PAIR may rotate inside
  itself (Ti's SA roots 0-1: angle 0.17, the pair span invariant).
- the meters went MATRIX-FREE (detbasis.operator_apply -- bitwise
  equal to the dense values on C, 0.0 over 69 roots): the dense
  meter build on the M-shell pol twins exploded to ~1.35 TB (the
  door's price in logs/collect_sc.tmp.log; the 3d3p4d@semi pol twin
  PRICED through it).
- the probe instrument went term-aware: the controls read L(L+1)
  and scale with the BLOCK's own non-degeneracy (Sc's SA residual
  2.3e-7; Ti's 3F block: 3.08 mEh split -> controls 3.7e-2); the
  gate's probe predicate carries the per-L rule (L >= 2 admits
  LARGE as the declared spread disclosure).

Sc: 52 declared -> 43 rows / 9 refused (4 Lowdin construction
drifts on the semi legs -- the store predates today's strict
Lowdin; 5 genuine substrate drifts -- the regen loop, gated).
Ti: 46 -> 41 rows / 5 refused (four lane legs OUT: the same-l pairs
3d4d@ar / 3d4s4p4d@ar -- the distilled window's S-orthonormality,
the declared future seam, their stores schema-stale -- and the semi
legs, the Lowdin needle on fresh construction; the CI lanes
regenerated fresh).  Sc-unc 25/1, Ti-unc 24/2.  The uncs' LANE legs
declared out: the lane substrate's basin steer cannot build today
(the plain D2h draw lands a wrong basin, the C1-seeded retry
diverges).  regen_cells.py: the run_queue/_materializes/keep/
parse_args block restored (lost in the lane_substrate move; the
unit tier never imports scripts, so it stayed silent).
Gates ALL GREEN: logs/validate_pt2_Sc.log, _Ti.log, _unc.log.
READINGS: the lane rows carry the core-partition axis (the CI lanes'
det-class disclosure: ptp 262/152 mEh on C/Si; the SA lanes' ptp
17/5.5 mEh); the frozen rows explode as declared (Ti-unc's @r1
frozen: -710.27 mEh); the route axis reads the orbital-draw
sensitivity at mEh scale.

## Status: step 4 (Si) DONE 2026-11-19 -- 34 cells / 104 roots; C+Si gate ALL GREEN

Si through the same chain: payload + table + CSV + the combined gate
(logs/validate_pt2_C_Si.log -- 170 roots certified).  ONE refusal
during collection (si_3p@s2 sa-casscf ao) drove a SEAM discovery:
the stored SA vectors are an INTRA-block rotation (full-rank cross
overlap, energies exact) carrying solver tails that move the block
projector at 2.2e-7 -- above the amplitude band the projector check
was borrowing.  The projector got its OWN band (p_tol 1e-6, the
tail scale; pinned: the tail refused at default/admitted widened,
o_tol forced off so the clause is exercised; the pin needs a
4-orbital window -- at the fixture's 3-orbital size the block spans
the whole sector).  The regen loop ran again (scoped --group planned
--window 3p@s2) and the row stands: -54.5525 mEh, <L2> 2.000000.

THE SI READING (mirrors C): the optimized cells saturate (all
windows price -55.79..-55.93 mEh; @s2: -54.4..-54.6 -- the bigger
window absorbed part of the dynamic correction); the thin span
fuels (compact-ao/3pno1: -236.09 mEh); E(2) tracks the orbital draw
(moc <L2> = 0.040 yet E(2) = -55.9263 = the clean spans' value);
the frozen-core row explodes (-226.08, 4x its free twin).  The
C/Si correspondence is cell-for-cell -- the axis structure, not the
atom, decides the reading.  (Si compact fuel -236.1 mEh vs the SA
-55.9: PT2 closes 180.2 of the 221.5 mEh poverty gap; the static
remainder ~41 mEh.)

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
