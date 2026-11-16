# Ticket 38 — the U14e declaration close-out: the tiers retire, the batch completes

**Status**: ACTIVE (spec synthesized 2026-11-14 per the to-spec /
to-tickets / codebase-design passes; the plan draft lives at
`~/.opencode/plan/u14e-closeout-spec-and-tickets.md` until promoted).
**Depends on**: ticket 36 (batch A's phases; the pacing probe ran).

## Problem Statement

The U14e case declared a run set its own construction could not
deliver: the band-lattice tier windows (`_U_TIER` — 3f2d, 3f, 4f,
2f2d, 4f2d; the case's PRIMARY was 3f2d = `9f10f11f12d13d`) carry
same-l multi-band token sets, and since the distilled declaration every
fresh cell builds through `windows.window_distilled`, whose assembled
S-orthonormality certificate refuses same-l band blocks (deviation
7.75e-01; the member-preserving cross-shell orthogonalization is the
declared future seam).  Measured 2026-11-14 (map Issues): the grid CI
pass crashed on the primary's mo cells; an unfiltered regen was
impossible; 36 of the wave's 111 priced optimizer legs were
phantom-runnable.  The user's direction (2026-11-14): the
CAS(4,20+)-class windows are not wanted — "what for?  we simply work
with basic discussed cases for now."

Companion: the no3 certificate cells were half-refused — the 7s pair
at the declared `no_cap=12000` (the 16731-determinant prefix = the
whole window re-gauged at (5,1)), the 14p pair never attempted (the
crashed queue skipped it; it is UNDER the cap, C(15,4)=1365).

## Solution

The declaration EQUALS the intent: the declared key-set is the
promoted windows only, with the primary re-declared at `10f13d`
(render_primary flips to True — the primary IS the results window now;
the suppression belonged to the retired tier primary).  The tiers
retire from the run set exactly as they retired from the render
(2026-11-06) — the store and git history keep their records; no guard
predicate, no new seam (Refactor 5: delete the fake degree of
freedom).  The no3 ladder closes honestly: the under-cap 14p pair
runs; the over-cap 7s pair refuses PER CELL (`natural.NoCapExceeded`,
the batch-D lesson — a construction refusal must not kill the chain)
and stays pending as the declared class.  Batch A then completes at
the MEASURED pricing: the pacing probe put the max-cost legs at
17:09 (single-root CASSCF) / 26:53 (17-root SA) wall at OMP=128 —
the optimizer batch is ≈ one day, not a week.

## Implementation Decisions

- `cases.py`, the U14e row: `shell_patterns=["10f", "13d"]`;
  `_U_TIER` deleted (the promoted windows are the whole declared set);
  `planned_windows=_U14E_PROMOTED[1:]` (the primary owns its grid);
  `window_names` loses the tier entry; `render_spans` untouched (the
  `("5f6d", "5f6d7s")` entries match the LANE spans, verified);
  `nroots_primary=48` stays (clamped to the sector).
- `shell_radials` keys stay AS IS — the lane construction reads the
  chemical keys; trimming risks KeyError (inert extras are safer).
- The `NoCapExceeded` refusal: raised at the owning seam
  (`natural.no_group_windows`'s cap check), caught per cell in the
  regen queue (`[refused]` print, nothing stored, the pending row
  holds the place).
- Store mechanics (DECLARED): after re-homing, an unfiltered regen
  rebuilds the primary grid ONLY — the planned addresses prune under
  the replace semantics; the canonical rebuild = grid pass +
  `--group planned`.
- Out of scope: the cross-shell seam (recorded future step), any
  CAS(4,20+)-class run, raising `no_cap`.

## Testing Decisions

- The derived pin: U14e declares NO same-l window (a new declaration
  row that repeats an l fails the pin, not a chain).
- The void pin recounts (184/36 → 138/27); the planned-sectors pin
  re-homes (the primary not a planned window); the chemistry pin
  flips with `render_primary`.
- The grid smoke pass + the planned CI pass = the tracer bullet's
  demo (green, the refusals per-cell).

## Progress

- [x] **T1** — the primary re-declaration + the tiers retire + the
      `NoCapExceeded` per-cell skip + the pins + the smoke pass (the
      grid green; the 14pno3 pair ran; the 7sno3 pair refused
      per-cell; render restored at 103 rows via the
      `render_primary` flip — the regression the smoke pass caught).
- [x] **T3a** — the (4,12) section's optimizer legs (the `10f13d` +
      `10f13dno2` classes) — LANDED 2026-11-16: 74 min wall at
      OMP=128; 24 computed, 4 the DECLARED radial-mask construction
      refusals (the case is distilled — the radial mask has no
      member-labeled rows), nothing lost; every optimizer below its
      own window's CI anchor, all energies finite (the audit's
      first-pass flags were the no-form windows having no `casci_ao`
      anchor — the correct `casci_no` anchors clear them all).
      Count correction (the record's "~73" vindicated): 108 optimizer
      cells declared, 18 stored pre-T3, 90 fresh = 74 computable +
      12 radial refusals + 4 no3 cap refusals; the computed total is
      74, matching the recorded estimate (my plan's 86 over-counted
      the radial class).
- [x] **T3b** — the (6,13) section (the `10f13d7s` family) — LANDED
      2026-11-16: 80 min wall at OMP=128 (the unloaded box — the
      probe's 17–27 min/leg priced a loaded machine; the payloads
      prove the work: full analysis fields, n_det 16731, all 22
      optimizer legs converged, every one variational below its CI
      anchor); 14 [refused] = the 10 declared no3-cap refusals (the
      frozen-mask no3 variants refuse too — the plan's arithmetic
      undercounted them) + the 4 declared radial ones; nothing lost.
      Store 85 → 107 data payloads; pendings 26.
- [ ] **T3c** — the (4,15) section (the `10f13d14p` family + the
      no3-form optimizer classes) — the user's enable; the batch-A
      close.
- [ ] **T4** — the review round over the wave (both axes + the
      design survey).

Dependencies: T5 (Stage 1c, ticket 37) LANDED 2026-11-14 — T3 is
UNBLOCKED (the three sections run on the user's enable); T6/T7
follow ticket 37.
