# Ticket 31 — The ECP wave: declare the ECP atom

STATUS: queued (opens after ticket 32's close; the door plan SETTLED
with the user 2026-11-05).

## Goal

Declare the ECP atom in `cases.CAMPAIGN` — ECP basis + pseudopotential
through the existing `AtomConfig.ecp` pairing — and retire the `run_u*`
scripts into the modern path (the U tables become the declared
campaign's output; `run_u*` stays the frozen record until then).

## The U14e door plan (settled with the user, 2026-11-05)

1. THE ELEMENT DIR IS `elements/U14e/`: rename the CAMPAIGN key
   `U` -> `U14e` (regen derives the element dir from the key -- that
   is the whole trick); the `tag` stays `U_stuttgart78_dyall` (the
   basis/psp label on chk names).  `elements/U/` STAYS FROZEN (the
   user's call; archive-or-keep revisited later).
2. CACHE-WARM: the record's reference chk
   (`elements/U/cache/U_stuttgart78_dyall_rohf.chk`) copies ONE TIME
   to `elements/U14e/rungs/rohf.chk` (the trusted-cache convention --
   cached chks never re-gate), or a gated fresh SCF.
3. THE UHF CENSUS: five arms at Dyall size (~30 min budget) ->
   declare `occ_e_ref_uhf` -> `scripts/uhf_anchor.py` -> the anchor
   row.  ONE FLAG: the census arms run without the Newton accelerator
   U's ROHF class needed (fact 22) -- if an arm diverges, the census
   learns `ref.conv_accelerator` (extend only on evidence).
4. THE FIRST REGEN IS ANCHOR-FIRST (~20 cells, NOT the full grid):
   declare an `anchor-first` Scenario for U14e (Ti's pattern) --
   the primary window's CI cells + the raw detector run; the four
   tier windows run their CASCI anchors only (~2/tier); SA and the
   optimizer variants stay VOID at the key (the multi-solution
   class; `run_casscf=False`'s spirit made explicit), one `--full`
   away.  U has NO scenarios today -- an unfiltered regen would
   charge into the 33-cell grid + all 132 tier cells.
5. THE RETIREMENT (G5): `run_u*` + `run_casscf`/`run_sa_casscf` die
   with the first modern U tables; `stage_facts`' `c_act` absorbs
   here (the one caller-change moment).
6. Ticket 29's T2/M batches stay user-gated.  The .chk restart
   anatomy needs NO work (every rung's chkfile IS the restart
   artifact: orbitals + CI vectors persist per address).

## Design (pre-verified)

- The seam already exists: PySCF's `get_hcore()` includes `int1e_ecp`
  when `mol.ecp` is set; h1 is built in exactly two places
  (`operators.active_space_integrals`, `operators.full_space_integrals`).
- The ncore-from-Z trap is CLOSED: `system.ncore_of` takes
  `mol.nelectron` (core electrons are removed under ECP — never
  recompute from nuclear charge).
- Optional explicitness: a named `h1_ao(mol, mf)` helper in operators
  would make the seam visible — decide during the wave.
- The ECP machinery is validated (pairing + p² and TM campaigns:
  ccECP + ccecpccpvdz; `tests/unit/test_ecp_wiring.py` +
  `tests/integration/test_ecp.py`).

## Retirement riding this wave

- `run_u*` scripts — the tables regenerate through the modern path.
- `casscf.run_casscf` + `run_sa_casscf` — pass-throughs carrying the
  stage name; the `name=` becomes explicit at the cases call site.
  `run_casci` KEEPS (real content: the CASCI class, `nroots`,
  `fci_opts`).

## Open physics (runs gated on the user's enable)

- **The U reference-landscape survey**: the gated Newton reference pins
  the ⁵S° basin (a numerically pure L=0 determinant); a
  symmetry-BROKEN determinant may sit LOWER (the generic mean-field
  pattern — every other campaign atom's reference IS broken). Run
  UNGATED, from anisotropic / m-localized starting occupations, both
  Newton and DIIS paths; ask whether any converged determinant lands
  below −50.38185233, and diagnose with ⟨L̂²⟩ + the closure
  certificates. Either answer sharpens the record: a lower broken
  basin re-reads every "below ROHF" anchor in the U tables; none
  found means the symmetric eigenstate-determinant IS the mean-field
  minimum (a noteworthy exception).
- The bare-shell lemma (deterministic, pure numpy — pin it in the
  unit tier when this ticket touches the code): NO single determinant
  of f³d¹ is a rotational scalar — the minimum ⟨L̂²⟩ over all
  all-alpha 4-electron occupations of a bare 1s+1f+1d shell is 6
  (exact 0 only for the half-filled parallel controls p³/d⁵ — the
  top-exterior-power scalars). Hence the ⁵S° purity must live in the
  RELAXED l-mixed orbitals (the wedge of the occupied subspace), not
  in any m-occupation — which is why the mechanism is interesting.

## Progress

- Part 1 DONE (2026-10-31, code-only): U DECLARED in `CAMPAIGN` — the
  raw specs load from `data/` (parsed once; the explicit tag carries
  the identity), the reference spec is the gated-Newton one, the
  window tiers are declared as their exact shell subsets (this basis
  is its own exponent grid, so the record's bands ARE whole-shell
  sets: 3f = 9f10f11f13d, 3f2d = 9f10f11f12d13d = the primary, the
  other four tiers planned), sa_weights = the 17-member ground-term
  M_L set, nroots_primary = 48, run_casscf off (the multi-solution
  class). Code-only verification: the declaration resolves end to
  end (mol builds ECP-honest: nelectron 14, ncore 5; the sector and
  the label<->patterns round trip; the primary grid = 33 cells, the
  planned tiers = 132). NO SCF anywhere.
- Part 2 EXECUTING (2026-11-05, the door plan; commits 408387b,
  bca1057):
  1. THE KEY RENAME: `CAMPAIGN['U14e']` (symbol stays `U`; tag stays
     `U_stuttgart78_dyall`); `elements/U14e/` born, `elements/U/`
     frozen. CaseSpec gained `key` (the element-dir contract, one
     lookup: `campaign_spec`); CellKey admits campaign keys.
  2. CACHE-WARM: the record's chk copied to
     `elements/U14e/rungs/rohf.chk`; the trusted load reads
     E = -50.38185233 (the record's value).
  3. THE CENSUS (five arms, Newton class — the door plan's flag
     SPENT on evidence: all five plain-DIIS arms diverge,
     E ~ +2.4e3..+3.0e3, so the census learned
     `ref.conv_accelerator`): FOUR arms identical to 1e-6 —
     E = -50.38881343, conv, <S^2> = 6.2289 (exact quartet 6;
     mildly broken), sigma_min = 0.8835; huckel non-converged (not a
     second basin). SINGLE BASIN, and 6.96 mEh BELOW the ROHF — the
     open physics question ANSWERED for the UHF class: the
     mean-field minimum is the mildly-broken determinant, not the
     ⁵S° eigenstate-determinant. `occ_e_ref_uhf = -50.38881343`
     declared. (The ROHF-class broken-landscape legs stay gated.)
  4. THE FRAC REFERENCE HAS NO CONVERGED BASIN on the ECP atom
     (measured): DIIS diverges (+137); Newton collapses (-12.74 vs
     the aufbau RHF -49.04); Newton-from-ROHF does not leave the
     seed. The `frac-reference-unconverged` scenario voids the five
     frac-fold cells; MOM-anchored selection is the recorded fix
     (map Issues). The regen head builds the frac reference only
     when a non-void frac fold runs; both head guards now unlink a
     diverged reference's chk (the trusted-cache poison class).
  5. THE ANCHOR-FIRST REGEN: the gate voids the whole optimizer
     class (145 cells, all windows incl. the NO forms) + the frac
     class → 15 CASCI anchors run. IN FLIGHT (detached; priced on
     the real primary cell: 272 s, kernel-only payload class).
- Part 2b: THE RESTRICTED OPTIMIZER BATCH (2026-11-05, user-directed:
  "145 too much, restrict"; commits 2b3624c, f3fa9ce): three
  promoted-pair windows declared over the band probe's identities
  (elements/U14e/logs/probe_bands.log -- the valence f/d character
  sits on 10f/13d; the CASCI pair ranking picked 10f+13d, the only
  pair hosting the f3d1 manifold, ground 5I at <L^2> = 42 exact):
  `10f13d` (4,0) the record's minimal f+d; `10f13d7s` (5,1) the
  record's six13 -- the 7s pair promoted, the fold follows the new
  per-window sector (`CaseSpec.planned_sectors` -- window_grid reads
  it; ncore_of(14,6) = 4 = the full core minus the 7s pair, the
  record's own fold); `10f13d14p` (4,0) the p-entry probe -- the
  physical 6p trio's dominant band in-window (the Ti 3d4s4p4d
  pattern; band labels are exponent-grid coordinates -- the label
  "6p" is a compact core-like band the CI never touches, caught by
  its exactly-zero spectrum shift; the trio lives on 13p/14p/15p).
  The batch's cells: free CASSCF + free SA per window + the PRIMARY
  window's free SA (the 17x5L term cell -- the minimal windows read
  5I ground, the 5L root +22 mEh mixed at <L^2> = 64.7; the record's
  AO 1f+1d row read the same 42) -- materialized through filtered
  --full passes after the anchors are read.  The m-complete 6p
  promotion (CAS(10,15)) is NOT expressible -- the fold-tag grammar
  assumes one band per shell (map Issues: character-selected folds).
- Part 2c: THE REVIEW ROUND (2026-11-06; commits 29a6aad, cb81663,
  86c7f5e, 034b4b1): the ROHF RECHECK (fresh SCF = the same basin,
  |dE| = 4e-09; the open-shell gauge rotates per draw — the
  determinant ⟨L²⟩ = 40.85, neither ⁵S nor ⁵L; the identity
  construction certified against the detbasis matrix at 2.2e-16
  first; `logs/rohf_recheck.log`); the WINDOW-COMPOSITION AUDIT (the
  optimizers stay ON-CLASS in l-space — the SA windows exactly
  12.000/12, 13.000/13, 15.000/15; the "escape" folklore was same-l
  RADIAL mixing, gauge not leakage; `logs/window_audit_ao.log`);
  start=`moc` (select_valence_mo: the character-selected canonical
  window — the SCF's own valence + the virtual complement ranked by
  the named bands' S-projection; the energy-ordered slice stays as
  the control); the PINNED-ROUTE PROBE (the irreps AUTO-DERIVED from
  term=(8,"u"); the MO windows STRUCTURALLY OFF the pinned route —
  irrep-mixed columns, the dominant-irrep repair collapses the span
  (overlap 0.647); the AO windows need a 1e-3 repair (overlap
  0.9999); the D2h substrate = C1-SEEDED NEWTON — from-scratch
  Newton lands −41.2, the seeded restriction reproduces −50.38185233
  exactly; `logs/pinned_route_probe.log`).
- Part 2d: THE CHEMISTRY + THE AO-PUR ROUND (2026-11-06; commits
  52cc9db, d1b89a7, 230f08d): the render chemistry on CaseSpec
  (window_names — the table speaks 5f6d / 5f6d7s / 5f6d6p;
  fold_chemistry — the folds read 6s6p7s / 6s6p, not "Ne";
  anchor_note — "forced filling: 6s²6p⁶7s² closed, 5f³6d¹ open";
  render_primary=False — the primary family suppressed; the tier
  families retired; the heading trimmed, the mapping sentence lives
  in README's reading key).  THE UHF ANCHOR WARMED (the census ran,
  the anchor script had not): E = −50.38881343 (occ_e_ref_uhf
  exactly), ⟨S²⟩ = 6.2289, σ_min = 0.8835.  THE AO-PUR VINDICATION
  (the user's diagnosis — the unclosed MO gauge breaks the pur
  optimizer start): the 5f6d AO-pur SA reads −50.423810 at
  ⟨L²⟩ = 72.0000 exact in 50 CYCLES (the fastest start of all: AO-pol
  56, MO-pol 160, MO-pur capped); the 5f6d7s AO-pur SA RE-REACHES
  THE L=9 STATE (89.9987, 114 cycles, 1.1e-08 — the MO starts miss
  the basin within the cap; the pol start needed 190).  The MO-pur
  rows stay as the recorded gauge finding.  The 6p window's MO rows
  landed (the echo CASCI + the broken optimizer class — the full
  AO/MO comparison on all three promoted windows).  METHOD COVERAGE
  CLOSED (user-reviewed): frac void (no converged reference), NO
  unrunnable at k=2/3 (recorded), the masks = the diagnostics class,
  SO/TM/g-factors = the future waves; the pinned route was the one
  missing piece.
- Part 2e: THE LANE FORK (2026-11-06): the pinned-lane declaration
  exposed the fold-grammar limitation — the chemical fold (6s6p7s)
  is 5 PAIRS spread over THIRTEEN radial columns on this
  uncontracted basis; spelled in band labels the electron-count
  contract refuses it (13 "pairs"); spelled chemically the shell-span
  fold breaks the 1:1 closed_shell_core_energy bookkeeping (it
  doubly-occupies every column — 26 frozen electrons).  The first
  route cells ran on the compact-band fold (the tag "6s6p7s" matched
  BANDS — compact core-like primitives) and read +1234 Eh — killed,
  the payloads deleted, the lane_windows declaration reverted
  (term=(8,"u") stays).  KEPT: the C1-seeded lane_substrate fix (the
  seeded D2h reproduces −50.38185233 exactly; plain DIIS diverges,
  from-scratch Newton lands −41.2 — measured).
- Part 3a: THE SHELL-LAYER WAVE EXECUTED (2026-11-06, the session
  after the plan; commits 34d9542, a64d8dd): the machinery —
  (1) `windows.shell_realization`/`realized_patterns` (the layer; the
  identity default = the contracted behavior, bit-exact; empty and
  other-l realizations refuse) + `CaseSpec.shell_bands`, threaded
  through `_fold_span`/`run_cell`/`analyze_cell`/`kernel_cell_payload`/
  `run_lane_cell`/regen; (2) `core.confined_core_orbitals`: the
  closed core ON a fold span = the SPAN'S OWN closed-shell Roothaan
  filling the lowest n_pairs (the first design — the span-truncated
  projection of the SCF core — DIED on the measurement: the S-weight
  capture is a NORM statement blind to the nuclear cusp, 99.99%
  capture with +13 Eh core-energy error on Sc-unc; the Roothaan is
  draw-free, monotone in the span, the span itself on 1:1, and on the
  whole AO space it reproduces the (2·n_pairs)-ion RHF EXACTLY — the
  pinned external anchor); the convergence certificate = the
  density-Fock commutator (the energy settles quadratically earlier);
  (3) the fold span = the UNION of its bands' rows (shared bands
  across shells duplicated rows = an exactly singular Gram — pinned
  live); (4) the pinned Davidson can deterministically SKIP
  near-degenerate lowest roots (the Sc-unc Ag ground pair: the
  solver's second root four eigenvalues up, reproducible) — the lane
  readout now DENSELY ADJUDICATES every CI block under the
  adjudication budget (12k dets at this commit, raised to 40k in
  Part 3b) with our
  exact H (the conventions certificate promoted to the correction
  source; `anchor_corrected` records the replaced blocks); (5) the
  per-root meters ride the apply family — matrix-free at any sector
  size (S²'s ladder terms are TWO-body — pinned).  THE DECLARATIONS —
  keys `Scunc`/`Tiunc` (the CellKey grammar is alnum; hyphens collide
  with the capitalize canonicalization), the split bases
  data/{sc,ti}_ccpvdz_unc.txt (scripts/make_uncontracted.py,
  self-verifying), the shell layer as CONTIGUOUS BAND RANGES
  partitioning each l-family (disjoint; the fold families FULL down
  to the cusp — the Roothaan fills the lowest, so a fold span must
  exclude every skipped/active shell's bands; the reference is the
  PLAIN aufbau ROHF — the composition log certifies its filling; the
  l-span forcing DIVERGES on this basis, measured).  THE SC VERDICT
  (elements/Scunc/): the free pur-fold row ²D exact at closure 8e-14
  (CAS(3,43), E = −759.7467, the certified contracted record
  −759.7274 below it variationaly); THE LANE 3d4s@ar = the
  multi-column fold exercised: 9 pairs over the 65-column family
  span, the Roothaan core 18 cycles at E_core = the free-ion value
  (−758.2147), ALL FIVE ²D components degenerate at −1.502958 (ptp
  3.7e-12, ⟨L²⟩ = 6.0, ⟨S²⟩ = 0.75 exact), total −759.7177,
  anchor_corrected ['Ag'] honestly recorded.  THE BAND MEASUREMENT —
  scripts/band_composition.py productized: per-band Mulliken
  populations (the display meter — WILD (±2) on near-dependent
  primitives, partition-unity-checked) + the greedy projector-capture
  suggestion (the certificate meter, bounded); the logs are the
  declaration's evidence (elements/{Scunc,Tiunc}/logs/).
- Part 3b: THE TI VERDICT (2026-11-07): the (3,1) sector over the
  full 8-band d family is 530k determinants (hours per cell — the
  first attempt killed at 75 min on one cell); the Ti window
  realization redeclared as the DOMINANT-5 d bands (the greedy's
  0.999-capture set; 28 columns, 92k dets) — the FOLD families stay
  full (the cusp is the fold's business).  _ANCHOR_CAP raised to
  40k dets (the box holds 3 TB — the 12k guess was timid).  THE LANE
  3d4s@ar: the Roothaan core 17 cycles; ALL SEVEN ³F components
  degenerate at −3.088057 (ptp 1.3e-13, ⟨L²⟩ = 12.0, ⟨S²⟩ = 2.0
  exact), total −848.2786 — and the adjudication corrected ALL FOUR
  blocks (Ag/B1g/B2g/B3g, dev 0.127): the Davidson skip is not an
  Sc accident but the pinned solver's CLASS behavior at these block
  sizes; the correction machinery carried the lane.  THE COMPARISON
  vs the certified contracted records — the CONSTRUCTIONS exact by
  certificate (²D/³F exact on both bases; the fold electron counts
  chemical); the TOTALS model-class (the unc lanes sit above their
  contracted twins: the confined ION-core vs the contraction core;
  the free pur rows sit BELOW the contracted free records —
  variational, bigger windows).  BOTH ATOMS PASS.
- Part 3c: THE U14E PORT (2026-11-07; the review round first — both
  axes, the findings fixed/declined/disposed in map.md; the moc
  layer bypass was the one real bug).  The declaration:
  shell_bands from the probe_bands record (5f on 9f-11f, 6d on
  12d-14d, 6s on 13s-15s, 6p on 13p-15p, 7s on 18s-20s) + the lane
  leg 5f6d@6s6p7s at (4,0); the band-labeled windows keep their
  identities (an undeclared pattern passes through the layer
  unchanged — the port is purely additive).  TWO findings on the
  first run: (1) the Part-2e FORK RESOLVED — the chemical fold
  6s6p7s = 5 pairs over the 15-column family span, the frozen core
  the span's own Roothaan (17 cycles, comm 2.6e-11), the C1-seeded
  substrate reproducing −50.38185233 exactly; (2) the declared ⁵Lᵤ
  (L=8) is NOT this window's block-ground — an L=6 quintet sits
  6 mEh below (the block-lowest readout interleaved the terms; every
  root a clean symmetry eigenstate at floor residuals) — the
  adjudication now TARGETS: the dense spectrum holds every root,
  the term meter (⟨L²⟩ per eigenvector) SELECTS the block's lowest
  `mult` components OF THE DECLARED L (the fallback records itself
  in the payload's l2).  THE VERDICT: all SEVENTEEN ⁵Lᵤ components
  degenerate at −3.106617 (ptp 3.0e-14, ⟨L²⟩ = 72.0, ⟨S²⟩ = 6.0
  exact), total −49.0374, anchor_corrected all four blocks.  The
  store-anchored pin carries all three validation lanes.
  THE USER'S DESIGN PRINCIPLE (2026-11-06, do not re-open): AOs are
  never tracked by labels directly — the chemical shell with its
  DECLARED AO realization is the primary object — FOR ALL BASES.
  MEASURED: cc-pVDZ Sc/Ti carry three d bands and six s bands (the
  chemical 3d = at least the "3d"+"4d" contractions) — contracted
  bases do NOT automatically hold the 1:1 shell-to-block match; the
  certified Ti/Sc lanes were 1:1 because their folds (ar/semi) named
  single-block inner shells.  The shell layer + the multi-column
  fold + the projected core energy are the GENERAL construction; the
  1:1 case is the certification special case.  QUEUED REVIEW (the
  user, 2026-11-06): the CASCI calculations (the window-poverty
  readings, the echo roots, the raw-fold anchor mixing) get their
  own review pass later — record the verdicts before any
  construction rides them.

## References

- The campaign record: `elements/U/` + the committed `run_u*` output.
- AGENTS facts 22 (the campaign verdicts) and 15 (the open-shell
  rules); the architecture snapshot: `wayfinder/map.md`.
