# Ticket 39 — the U14e prior axis: the @r/@s windows over the full-band families, the poverty question

**Status**: ACTIVE (approved 2026-11-17; execution per the Stages below)
**Depends on**: ticket 37 (complete — the prior-axis machinery), the queued CASCI review (executed 2026-11-17 — the poverty class classified)
**Origin**: the 2026-11-17 plan session — the user's picks: U14e first, CASCI only now (CASSCF planned as the follow-on), @s2 included despite the (4,24) analysis-ceiling class.

## Problem Statement

U14e's incumbent 5f6d D window is band-limited — the declared bands hold the
compact primitives (10f+13d), and its CASCI rows sit at +923 mEh with EXACT
closure certificates: radial poverty, not a symmetry defect (the review's
classification). The prior family (`@r<k>`/`@s<k>`) reads the density block over
the FULL chemical-shell families (5f → 9f-11f, 6d → 12d-14d, per the band
certificate) — the same axis that certified @s2 mean-field-sound on the twins.
The wave measures whether that axis cures the poverty, or the band-identity
poverty is structural.

## Solution

Declare the U14e configuration series (head f³d¹ + the Hund-quintet members f⁴
and f²d², all spin 4, discovery-then-declare), declare `prior_bands` = the
certified band families, run the CASCI batch 5f6d@{r1,s1,s2} × {raw,avg,pur}
(9 legs, (4,12)/495 dets; @s2 = (4,24), the kernel-only class), re-render, and
read the +923 question off the table.

## Stages (gated: fast tier + targeted files; one regen at a time)

- **Stage 0** — the TABLE.md prose notes removed (the user's directive: the
  columns carry the analysis; README's key owns the reading rules). `232c3bb`.
- **Stage 1a** — the `run_series` seam fix: member configs thread the declaring
  atom WHOLESALE (`dataclasses.replace`); the constructor call dropped
  `ecp`/`charge`/`tag` — an ECP atom's members would have built all-electron.
  Micro-pin + the head test on the real `AtomConfig`. `e9711c8`.
- **Stage 1b** — the series census (discovery): plain arms (accelerator × the
  gate's guesses) + `dm0=`-seeded arms for f²d² (the hook cannot prefer it — it
  shares the head's spans); fingerprints (per-l weight sums over the open block)
  discriminate f³d¹/f²d²/f⁴ where the character post-check cannot. Measured
  `e_ref`s + recipes become the declaration. Log: `elements/U14e/logs/
  series_discovery.log`.
- **Stage 1c** — the declaration: `series`, `prior_bands`, `planned_windows` +=
  5f6d@{r1,s1,s2}, the U14e prior-anchor-first scenario (the optimizer cells
  void — the comparison is CASCI), the certification runs writing the trusted
  series chks, the declaration pin.
- **Stage 2** — the CASCI batch (9 legs, threads 16–32, well under an hour) +
  `build_methods_tables.py U14e`.
- **Stage 3** — the read: does the prior cure the +923?  YES → the CASSCF
  follow-on gets its own go/no-go (entangled with the adoption decision).  NO →
  the band-identity poverty is structural; the record says so.  Map updated,
  ticket retired either way.

## Risks declared (2026-11-17)

- ECP-member SCFs are untested territory for the series machinery — the seam
  fix (Stage 1a) removes the known hazard; the discovery arms are the test.
- Member energy orderings quasi-degenerate-class — never pinned, draw-documented.
- @s2's d-shell k=2 needs the f²d² member (the head's d block is rank-1); a
  refused member refuses @s2's construction — the honest declared refusal (the
  C class), reported, never masked.
