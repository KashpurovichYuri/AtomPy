# Ticket 35 — the basis-preparation layer (the universal construction)

The wave: ONE construction everywhere — the contraction recipe as a
declared layer (stage 0b) between the reference and the treatments.
Opened 2026-11-10 from the handoff (xxii) skeleton + the review round
(fce521a).

## The design (user-settled, do not re-open)

- The DIVISION: basis PREPARATION (atom+primitives -> the contraction,
  its own certificates) vs TREATMENT comparison (the existing
  machinery unchanged over the FROZEN basis).  Within a table one
  basis; across tables the basis is a declared axis value.  The
  per-target radial optimization is the radial-CASSCF TREATMENT row,
  never a per-cell silent contraction.
- The source axis: the @r/@s window tags (the density-block family,
  `basis.series_window`) beside the incumbent untagged windows
  (`window_distilled`, the reference-source span-signal family).  Two
  declared families, never mixed in one comparison row.  k from the
  tag; the kept-signal tail certificate refuses a beyond-signal k
  (single-config blocks are rank-1 -- @r2 does not exist).
- The ladder: basin-gated forced-occupation references, the HEAD = the
  case's own trusted reference (never re-run).  The k-frozen rows
  stand alone (symmetry-exact, suboptimal); radial-CASSCF applies to
  k=1 and k=2 alike (the k=1 legs serve capped sectors).
- NO BASIS.md artifact (user): the prior comparison renders where the
  cells live; the k evidence/captures/discarded-fraction disclosures
  ride the payloads + table preambles; the exporter round's basis file
  becomes the record.

## Done (92bedbf + fce521a + 56bc03c + a4e61a8 + 2b4a291 + this wave)

- The review round (fce521a): the corrected record, the probe's dead
  code gone, the anchor circularity documented.
- `basis.py` + the grammar + the Sc/Ti declarations + the anchor-first
  gates (92bedbf).  THE RADIAL FLIP REPAIR: the probe's metric
  transport solved (block, s_r) -- Euclidean-pairing, distorted
  radials; the honest capture-maximizing problem is (block, s_r^-1)
  (the planted-radial pin caught it).
- The census (56bc03c): rohf_basin_census + the occ_ok certificate.
  THE NEWTON DISEASE: the second-order path freezes occupations from
  its early cycles -- converged arms can VIOLATE the hook (a
  3d-forced quartet landing the d2s1 solution); run_ladder
  post-checks.  The plain-DIIS arms diverge (LadderSpec.accel).  The
  true quartet basins: Sc d2s1 -759.69899090, Sc d3 -759.57064777,
  Ti d3s1 -848.38650790, Ti d4 -848.24887522 (the probe's draws
  155-177 mEh shallow garbage).
- THE FIRST HONEST PRIOR TABLE (Sc; Ti's 3d4s rows match the pattern):
  @r1 +102.9/-7.0 (3d/3d4s), @s1 -64.7/+5.5, @s2 -119.5/-16.4 -- the
  k=1 hierarchy the probe reported is INVERTED (the ladder mean is
  the strong k=1, not the worst); @s2 reproduces the probe EXACTLY
  (span-invariance -- the internal consistency certificate).
- Structural voids: the tagged windows carry no frac fold; Ti's
  same-l planned windows voided on the cross-shell seam (latent since
  the distillation wave); the 3d-only prior windows at the LATTICE's
  d2 sector (planned_sectors -- the case default (3,1) has no
  incumbent twin).

## Open (in order)

1. THE RADIAL-CASSCF LEGS (the per-state <= on the fixed basis):
   k=1 AND k=2 prior windows alike (the k=1 legs serve the capped
   sectors -- the relaxation where k=2 cannot go; the k=2 legs the
   maximal relaxation at the normal tier).  The optimizer cells of
   prior windows are void under the anchor-first gate -- materialize
   with `--group planned --window <label> --full` (the no-forms are
   their own labels).  Their rows join the dieted same-CAS blocks.
   Per-state <= check against the incumbent's optimizer rows.
2. SA-over-series (the transferable k=1 contraction); the twins +
   C/Si ports; the adoption gate (user); the exporter + convergence.
3. Still queued: mo x raw cells; the twins' UHF class repair (the
   class diverges on the uncontracted bases, measured 2026-11-11 --
   map Issues).

## Executed 2026-11-11, the review-stage rounds (the user's
## structure, delivered iteratively)

- (9274b04..f8e8a58) THE APPROVED RUNS, fully executed: every render
  span on every atom reads ZERO runnable-pending (Sc/Sc-unc/Ti/Ti-unc/
  C/Si; U14e deferred).  Sc-unc: the 3d lattice (30), the 3d@ar lane,
  the @r1 NO certificates.  Ti-unc: the @r1 priors at both sectors
  (18), the 3d lattice (32), the lane, the UHF anchor.  Sc: the 3d
  lattice (32) + the lane.  Ti: the kind blocks (12), the 3d lattice
  (30), the lane.  C/Si: the mo/moc legs (12 each).  THREE code facts
  the runs forced: (a) the @s guard went LABEL-level (the sector-
  bearing prefix leaked 18 Ti-unc cells past the planned_sectors
  mirror); (b) the radial construction refusal got ONE home
  (cells.radial_refused -- regen's queue died at the twins' lattice
  radial cells); (c) the same-l guard DERIVES from the declared
  windows (Sc's lattice ladder carried five unguarded same-l spans;
  the derivation also caught Ti's 3d4p4d).  VERDICT REVERSAL: both
  twins' UHF anchors CONVERGE spin-pure (the 'lands diverged'
  expectation was wrong); the diverged-note machinery stays for
  whoever needs it.

- (f2e5d04..328ee86) THE STRUCTURE ROUND: pending placeholders, the
  sectioned skeleton (=== anchors/span transitions, --- families,
  lane closes its section), the UHF place row, the NO class retired
  to store-only, render_spans review scopes, lattice_keys = the
  full grammar, 3d@ar lane windows on the pairs, Sc's preamble
  emptied.  (d9cc576) THE CAMPAIGN REPLAN: C/Si/Ti/Ti-unc on the
  skeleton (41/41/150/146); Ti-unc = the Sc-unc mirror (head-only
  series, guarded @s, mirrored prior sectors); the parity pins on
  BOTH pairs; U14e deferred.  The C/Si UHF records re-warmed (lost
  to the pre-fix merge bug).
- (3646175) README: the source axis explained accurately (the signal
  definitions, the two controlled pairs, the MO-informed/AO-form
  contracts).
- (d6375e6) THE SC LEGS: free + SA over all six prior windows at the
  ruler cell (free start-blind EXACTLY; SA @r1 = the incumbent basin
  in 34 vs 129 cycles; the SA stack anomaly surfaced) + the conv
  flag.  The radial subcase: undeclared on priors BY CONSTRUCTION.
- (cff3e96) THE SC-UNC @r1 PORT (the degeneracy certificate: D ≡ @r1
  on the twin) + the @s block recorded (the fresh-SCF disease).
- (c159c07) the 19-row comparison block (MOC's optimizer legs run;
  the SA-MOC stall found; THE STORE-MERGE BUG fixed -- regen had
  been deleting the UHF records).
- (4bfba78) the review-stage table: CaseSpec.render_spans; the NO
  optimizer legs rendered (later retired); the preamble compacted.
- (09922ca) the symmetric core groups (priors in every group; 12
  pur optimizer legs), the explicit per-start order, THE ROUTES
  FILTER BUG (the lane/pinned family was wrongly suppressed).
- (466b6de) the brevity round (notes 1-3 words; Sc's preamble
  emptied), the LOCAL-MINIMUM derived note (the SA-floor variational
  rule; fires on Sc/Ti lane rows + U14e's MOC rows), the NO pol
  cache repaired (old-wave 2-det solves -> fresh 90-det windows).
- (this commit) the NO start fully store-only (its window IS its
  source's span re-gauged -- the CI degeneracy + the optimizer
  equality are the certificate, not a row); Sc-unc's preamble
  dieted to the atom-specific minimum.

## Executed 2026-11-11 (the user's table review round)

- (82d6611) THE RENDERING DIET: the methods table renders the
  comparison, the store keeps everything -- NO-start cells never
  render (energy-identical to D), the prior windows contribute
  their distilled CI row at the ruler core only, the lattice one
  row per span, the frac class does not render, the prior labels
  verbatim (3d4s@r1 / D); Sc 122->65 run rows, Ti 123->70, C/Si
  32->29, the twins 14->10.
- (010688a) THE LOCKSTEP RENAME: fold -> core (CellKey.core,
  core_tag/core_shells/core_col/core_label/core_chemistry,
  _lane_core; the storage schema strings frozen by design), ladder
  retired (methods_section, the methods:begin/end markers,
  SeriesSpec/run_series/CaseSpec.series, rungs/series); the L+-
  physics and run_folded/folded keep their names by design.
- (9dd43fb) THE UHF ANCHORS: Sc/Ti/U14e warmed bit-exact to their
  declared gates, the anchor rows render; the twins BLOCKED (the
  UHF class diverges on the uncontracted bases -- map Issues).

## Executed 2026-11-10 (the table feedback round + the docs)

- (a4e61a8..1a2446d) the basis layer + the honest prior tables (the
  wave record above).
- (c82fa6d) the two-axis review's fixes: occupation_postcheck ONE
  home (public in reference; run_ladder and the census both call
  it), the census ref required, the measured instance -> Notes, the
  prior-tag grammar through cells.prior_tag_of, _no_label public as
  no_label, the reference layout row current.
- (b0e8c53) the same-CAS comparison: _comparison_key (incumbent ->
  r -> s1 -> double-budget s2 -> NO -> controls -> optimizers), the
  prior windows JOIN their span's family (a lattice span equal to
  the primary's renders once; a prior span without an incumbent
  twin renders as its own family -- Ti's 3d at the d2 sector);
  StageFacts.window_report (the construction's disclosure rides the
  payloads -- 108 prior cells regenerated fresh, energies
  bit-stable); the prior preambles in Sc/Ti TABLE.md; pins: the
  comparison order, the fallback family, the report carriage.
- (70c8517) the README restructure: every axis defined ONCE in the
  variety-axes table in one fragment shape (fold/core-field/
  substrate rows; the source row defines the configuration series
  and the head; the baseline protocol's one home; the label grammar
  as one syntax block); the three column legends merged; --93 lines
  of duplication.
