# AtomTerm — the effort map

The navigator: what is live, what is open, where the worked examples
live. History is git's job — resolved tickets are deleted, survey dumps
are not kept (open findings move to Issues below before their section
dies).

## Tickets

- [29 — gated measurement batches](tickets/29-ti-first.md): the T2 Ti
  cells + the M store migration; runs only on the user's enable.
- [33 — the symmetry-breaking anatomy deliverable](tickets/33-symmetry-anatomy.md):
  design-first — the frozen-landscape tables that replace the retired
  scan pointers.
- [34 — the folded reference + the window characterization matrix](tickets/34-folded-reference.md):
  **STEPS 1-5 COMPLETE (2026-11-08/09)**: the Tiunc census + the
  promotion; run_folded (reference.py) + its tests; the matrix on the
  folded field; the lane re-runs -- THE GAPS CLOSED (Tiunc +124.3 →
  5.5e-13 Eh; Scunc +38.7 → 0.085 mEh); the Sc port.  The riding-along
  tail (the anisotropy column, the per-radial spectrum, the docstring
  review, the table preambles) + the folded-field campaign-schema
  decision remain; then T-B/T-C/T-D/T-E per the ticket's After.

Resolved: ticket 32 (the symmetry lane) — the campaign complete
(C/Si/Sc/Ti on the one table, the route axis, the UHF anchors), the
wave closed 2026-11-05 on the gates (the full-suite run retired by
user decision; the partial verdict 91/91-isolation-green).

## Issues (open dispositions + watchlist)

- The folded-reference round, MEASURED (2026-11-08/09, ticket 34
  steps 1-5 -- the logs: elements/{Tiunc,Scunc}/logs/{census,
  folded_probe,matrix_folded}.log): (1) the unc ROHF basins are ALL
  s/d-HYBRID-open (~0.4 d per open MO; pure-3d opens are a
  contracted-subspace property) and no fresh start reaches the
  deepest basin -- the trusted cache is the one mechanism; (2)
  run_folded's pick rule is MOM-PRIMARY within dominance gates
  (character alone loses to the complement's tight core-shell
  radials at w = 1.0; energy alone declares the ion field's
  d-BELOW-s inverted filling, eps -1.53 vs -1.19 on Ti); (3) the
  ion-field valence hybridizes (the folded closed pair ~0.69 s on
  both twins) -- the C1-gauge Scunc cut carries only 0.4628 4s-band
  weight, below the compression's 0.5 dominance gate (the D2h
  solver's gauge lands >=0.5; the neutral field's signal is 0.8988)
  -- a free-route folded cut needs a declared dominance (0.4
  measured honest: L2 = 6 exact); (4) THE GAPS CLOSED: Tiunc
  +124.3 -> 5.5e-13 Eh, Scunc +38.7 -> 0.085 mEh (the residual = the
  C1/D2h gauge of the window cut); (5) the MO-picked windows fail
  the meters on BOTH fields (the m-incomplete set) -- the AO class
  is the certified default; the NO arm degenerates with the source
  at k=1 (the groups [1,5] = the chemistry).  OPEN DECISION: the
  folded-field cells' campaign identity (the substrate tag in the
  key grammar, the free-route dominance declaration, the table
  re-anchoring) -- user-visible before the tables carry the folded
  rows.
- The projected-Gram alternative, DECLINED (2026-11-08, the design
  round): forming the compression Gram from S-PROJECTED member
  vectors (bounded spectra, self-interpretable shares) instead of
  the truncated ones.  The cons beyond the recompute: the idiom fork
  (the truncated Gram is shared with the pur construction -- changing
  both is a global re-certification), the marginal gain (the post-hoc
  per-radial spectrum delivers the interpretability without touching
  the construction), and the unmeasured ranking difference.  REVISIT
  CONDITION: the per-radial spectrum showing truncated-ranking
  distortions (the riding-along tail installs the spectrum).
- The D2h substrate layer, UNDERSTOOD (2026-11-07, the user-directed
  probe -- elements/Scunc/logs/d2h_probe.log is the record):
  (1) the D2h ROHF landscape is the SAME multi-basin zoo as C1 --
  same-filling stationary points at the 0.3-11 mEh scale, each
  reverse-seed-stable on C1 (the D2h solutions are genuine C1
  stationary points, not symmetry artifacts); (2) the seed's job is
  the OCCUPATION basin only -- the solver's irrep-blocked diagonals
  erase the gauge and land the D2h-canonical solution (Scunc: plain
  and seeded land IDENTICALLY); (3) NEVER symmetrize a degenerate-
  gauge seed: pyscf.symm.symmetrize_orb AMPUTATES the open MO
  (measured: the open d became a pz, +183 mEh) and symmetrize_space
  scrambles the occupation-column pairing (measured: the runaway
  class) -- the user's MO-selection-by-irreps hypothesis evaluated
  and declined ON THE SEED SIDE; the solver already irrep-adapts;
  (4) the Tiunc basin zoo, CENSUS-RESOLVED (2026-11-08,
  elements/Tiunc/logs/census.log): two stationary basins on the unc
  basis -- the forced-Newton class -848.40648875 (the old trusted
  record) and the seed-refinement DIIS attractor -848.40744343,
  1 mEh deeper -- both with s/d-HYBRID opens (~0.4 d per open MO:
  pure-3d opens are a contracted-subspace property; every fresh
  forced start DIVERGES, so no ReferenceSpec gate is declared -- the
  trusted cache is the one mechanism).  The user promoted the deeper
  basin; the free cells re-run on it.
  Scunc: the trusted record is already the deeper of its two.
- The lane substrate's D2h basin wander (2026-11-07, the compressed
  round): the from-scratch plain D2h draw is basin-chaotic (the
  contracted Ti's plain D2h lands the RUNAWAY -848.2834 while Sc's
  hits its reference exactly) -- lane_substrate VERIFIES the plain
  energy against the trusted rungs/rohf.chk and retries C1-seeded on
  disagreement past 1e-5.  The seeded retry slides to its own
  solution -- which the census resolved: the attractor it slides to
  IS the adopted deepest basin (-848.40744343, reverse-seed-stable,
  reached identically from every seed tested), so after the
  promotion the retry lands the reference basin EXACTLY -- the
  substrate and the free route run the same field again.
- The Ti multi-basin finding (2026-11-07, the consistency probe):
  the Ti 3d2 4s2 ROHF hosts SAME-FILLING basins 11 mEh apart on the
  CONTRACTED basis itself (the record's gated -848.40648875 vs an
  ungated -848.3953); the unc plain aufbau lands the SHALLOW class
  (-848.39639) -- no variationality violation (the primitive sets
  are identical, the subspace exact; SEEDED from the contracted
  record through the cross-overlap, the unc basis reproduces
  -848.40648875 EXACTLY).  2026-11-08: the Tiunc reference is now
  the CENSUSED deepest basin (logs/census.log -- the seeded record
  itself re-derives: DIIS from it rolls to -848.40744343); the
  Tiunc cells re-run on it.  Fact 21 gains a measured instance:
  even at the right filling, plain aufbau can land 11 mEh high; the
  seed (not the gate) is the cure when the gate's cycling diverges
  on the basis.
- The queued CASCI review (the user, 2026-11-06; PRECEDES the next
  physics wave): the CASCI calculations get their own review pass —
  the window-poverty readings (+923 mEh class), the reference-echo
  roots, the raw-fold anchor <L2> mixing — before any construction
  rides them.  (From ticket 31, kept alive at its deletion.)
- The 61c5529..642ac8e review (2026-11-07, both axes; the
  shell-layer wave).  FIXED IN ROUND: the moc start bypassed the
  shell layer (`select_valence_mo` on raw patterns -- the unc moc
  cells selected single-band character; expanded through
  `realized_patterns`, the four cells rerun); the correction branch
  of the dense adjudication unpinned (the store-anchored
  `test_unc_lanes_load_exact_from_the_store` re-verifies the
  committed Scunc/Tiunc lanes -- ptp, L2/S2, core_pairs, the
  recorded corrections); `_fold_span` misnamed after the redesign
  (renamed `_fold_core` with its consumers); the ticket/code
  adjudication-cap mismatch; the duplicate `sys.path` in
  band_composition.  DECLINED (with reasons): the `shell_bands`
  threading as a data clump -- the substrate signature IS the repo's
  documented convention, and a context object for one new field is
  speculative (the store-vs-layer divergence the reviewer notes is
  real but bounded: the layer is case-declaration data like the
  basis itself, and the lane payloads now carry `core_pairs` to
  cross-check); moving the band-name derivation out of
  band_composition -- it is display-layer grouping, not AO selection
  (fact 1 governs selection; the docstring now says so); bundling
  the scripts' `_ANG`-class constants into windows -- no consumer
  variance, a hypothetical seam.  DESIGN SURVEY (the wave's seams,
  the codebase-design vocabulary): the shell layer is a deliberately
  small module whose value is LOCALITY (the chemical->band grammar
  in one home; the deletion test scatters it across six call sites);
  `confined_core_orbitals` is DEEP (a Roothaan, the cusp physics,
  and the ion limit behind three parameters; the interface is the
  test surface -- the unit pins and the ion anchor cross it
  directly); `_fold_core` concentrates the refusal family; the
  adjudication INVERTED the anchor's direction (our exact H
  adjudicates the pinned solver, both numbers kept) -- an epistemic
  shift recorded in the docstring and below.
- The adjudication watch (2026-11-07): the lane's independent-solver
  verification became SELF-verification under the adjudication
  budget -- the pinned solver's roots are a proposal our exact H
  corrects (deterministic Davidson skips measured on both validation
  atoms, all four Ti blocks).  Both numbers are kept (`anchor_dev`,
  `anchor_corrected`); the risk is a shared-convention error our H
  and the solver make identically -- covered by the dense commutator
  certificates and the ion anchor, not by solver independence.
  Revisit if a block class appears where adjudication is
  unaffordable (> the budget) AND the solver skips.
- The f41bb51..HEAD review (2026-11-04, both axes; fixes in the
  round's commit, ticket 32's closing entry is the record).
  DESIGN SURVEY (the delta seams; the full src/ survey stays at the
  wave close): reference.py got DEEPER (one cache/gate skeleton +
  one acceleration dispatch behind three runners; the object
  construction `_uhf_object` is now a REAL seam -- two adapters:
  the production runner and the census arms; the declined cells-axis
  seam is documented at the module contract, where a seam refusal
  belongs); the UHF gate's declaration rides the row's
  ReferenceSpec (the declaration-owns-intent rule held).  The one
  real shallowness finding is item (1) below.  The scripts'
  vocabulary module tightened (comm_text deleted).
  RESOLVED 2026-11-05 (the wave-close refactoring round, 6395de7,
  after the full-scope design survey): (1) the kernel/lane residual
  block found its one home (`_apply_family`; plus `_key_core`, the
  STRONGER triple-duplication the survey caught) -- both runners
  kept (the deletion test: the fold/cache/store paths genuinely
  differ); (2) the V1 splits done (the apply-boundary pin three
  ways, the ladder pin six ways); (3) the V7 citations captured on
  the frozen twins (and the capture CAUGHT a graze: test_symm_tm's
  l2 band widened 1e-8 -> 1e-6 with the observed 8.2e-9 cited --
  1.2x headroom on an SA-leg solver-noise quantity); (4) max_cycle
  rides the payload beside cycles (read off the live mc, flags by
  identity; the builder keeps the legacy ==200 key for stores
  predating the field); (5) cas_label's NO branch deleted (dead --
  cells makes the span+start combination unrepresentable and the
  builder filters start=="no").  Also in the round: the verbose
  duplication RESOLVED (the run_cell kwarg died; the budgets key is
  the one channel); run_carbon reads the DECLARED row (NOTE: its
  artifacts now label "C" where the pre-declaration path labeled
  "C_ccpvdz_s2" -- alignment with run.py's row-driven path, no
  committed artifact affected); scripts one-homes (payload_of,
  atomic_write_text).
- Declinations of the round (2026-11-05, reasons recorded here):
  stage_facts' c_act absorption -- rides the ECP adapter retirement
  (ticket 31, its one caller-change moment); the ECP term-exactness
  and symm_lane L^2 pins' V7 citations -- textbook-integer class,
  docstrings carry the physics, not the residual class the item
  targeted; the unknown-atom CLI guard dedup -- four lines, two
  CLIs, a dedup that adds an import for a TODO's worth of code.
- The `test_symm_tm` SA legs mispriced as tests (2026-11-05): the
  fixture's SA-level pins re-run campaign-scale optimizer pairs
  (1-2 h class) fresh every invocation -- a run in a test costume.
  RE-SCOPE (the AGENTS Tests boundary, user-directed 2026-11-05):
  fresh CI legs stay; the SA verdicts become STORE-ANCHORED (load
  the committed route-batch stores, re-verify the meters -- no
  optimizer reruns).  Lands at the file's next touch (the U14e wave
  needs no bump from it; the verdicts already live in the stores).
- The U14e frac reference's MOM-anchored fix (2026-11-05, evidence in
  ticket 31): the fractional block's cycle-to-cycle gauge is
  unanchored on the ECP atom (DIIS diverges; Newton collapses to
  -12.74 vs the aufbau RHF -49.04; Newton-from-ROHF does not leave
  the seed) -- the runner docstring's prescribed remedy is MOM-style
  anchoring of the character-selected block.  Until then the
  `frac-reference-unconverged` scenario holds the frac folds void.
  RIDES THE SAME VISIT: the diverged-reference poison-chk hazard
  (a non-converged SCF leaves a chk the trusted cache would load;
  the regen head's guards unlink at refusal -- the deeper fix is
  writers persisting only converged chks).

- Character-selected folds (2026-11-05, the promoted-pair wave): the
  fold-tag grammar cannot express a multi-band fold -- the basis
  carries the physical 6s on bands 14s+15s, the 6p trio on
  13p/14p/15p, the 7s on 18s/19s/20s (labels are exponent-grid
  coordinates, not principal numbers), and `closed_electrons`'s
  per-shell arithmetic assumes one band per shell (a tag like
  `@14s15s14p15p18s19s` counts 20 closed of 14 valence electrons).
  So the m-complete 6p PROMOTION -- CAS(10,15), the record's own
  "inner core-pair valence, to be done m-complete" open question --
  has no expressible cell: the occupied-prefix fold cannot skip the
  6p trio (m-set violation at ncore=2), and the p trio enters the
  window only as virtual columns (the `10f13d14p` probe).  The fix
  (when wanted): the fold span selected by orbital character (the
  dominant band set of the named MOs), with the closed count read
  off the reference occupations -- a real seam extension, not a tag
  hack.

- The ladder-wave review round (2026-11-01, `701b31c..c02d0c6`, both
  axes): applied — the module-docstring lockstep, the shared
  `_substrate` builder, the V1 test split, the `_FOLD_*` renames, the
  `LadderLeg` NamedTuple.  Declined: splitting the ladder sessions
  out of ticket 32 (the ticket IS the provenance mechanism, per the
  Conventions); the map queue-decision placement nit (live prose
  until the wave closes).  Resolved 2026-11-04: the emitter class is
  DELETED (the one-table builder owns the render; the sector stays
  implicit in dets/CAS by design).  Recorded partials (spec axis):
  the ladder
  one leg per construction, not per leg (the uniform-construction
  argument — revisit if a leg's meters ever disagree).
- Pinned cells' window intake (the symmetry lane's cross-route
  question, 2026-10-31): the natural pinned-cell windows are
  D2h-BORN sources only (raw shells, D2h-run MOs/NOs) — free-run
  windows never cross into pinned cells; the grammar states the
  restriction at the declaration wave.  The trap it avoids:
  `label_orb_symm` labels m-mixed columns by dominance, so a
  free-run window fed to a pinned solver would appear to work and
  produce a silently wrong solve — any future intake needs an
  explicit irrep-purity check at the seam.  The laundering routes
  (NO-density α, radial-transplant β) stay diagnostic experiments,
  not grammar values.
- `run_casscf`/`run_sa_casscf` are pass-throughs carrying the stage
  name — retire at the ECP wave together with `run_u*` (their second
  consumer); the stage `name=` lives at the call site (2026-10-30;
  the 2026-11-05 survey confirmed the inventory: one caller + two,
  both the ticket-31 retirement set).
- `stage_facts`'s `c_act` parameter: every production caller passes the
  same derived slice — absorb at the ECP wave with the adapter
  retirement (2026-11-05 survey re-confirmed; the deepening lands
  with its one caller-change moment).
- Ti escape comm band (test_treatments, the sa_free band): quiet
  through recent suites; capture-not-widen stands (no observed value
  recorded).
- The D2h twin's determinism certification is procedural (double-run
  sha256 in the session record), encoded nowhere in-repo; a generator
  `--verify` mode is the candidate — declined until a second regen
  cycle (2026-10-31).
- The D2h mol construction has two homes by design (the generator's
  draw, the fixture's load); a shared builder would serve one consumer
  each — revisit at a third frozen reference (2026-10-31).

## Examples (worked records in `elements/`)

- `elements/C/`, `elements/Si/` — the p² campaign: the methods grid,
  scans, multiplet tables.
- `elements/Sc/`, `elements/Ti/` — the transition-metal campaign:
  rungs, the radial axis, the d² lattice, term tables.
- `elements/U/` — the ECP-atom record (`run_u*` until declared).

## Architecture (edge-verified)

    pure math      system · detbasis · secondquant · natural
    substrate      windows (the hub)
    1-particle     angular (uses windows' shell_slots)
    assembly       operators (both PySCF boundaries) · irreps (the
                   D2h dictionary; the symm boundary)
    engines        reference · core · casscf
    analysis       analysis · stages · scans
    declaration    cases -> cells

`natural`'s mathematics is pure at module level (numpy + secondquant
only — the windows→natural edge retired with the window-kind
registry, 2026-10-30); its one driver (`no_group_windows`) lazily
imports the assembly layers. No cycles, no upward edges. Division of
labor over windows: declaration NAMES (labels), windows BUILDS (the
kinds + substrate builders), cells BINDS (the start axis picks the
builder).
