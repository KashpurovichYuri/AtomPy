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
