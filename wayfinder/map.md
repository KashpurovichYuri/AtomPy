# Wayfinder map — AtomTerm

## Destination

A small, tested Python project that takes an atom (start: carbon, cc-pVDZ)
through ROHF → CASCI → CASSCF → SA-CASSCF in PySCF with checkpoint caching,
classifies every obtained state by exact Lz / L² / S² matrices built in the
determinant basis (first quantization), and quantifies symmetry preservation
via ||[H, Lz]|| and ||[H, L²]||. Done = pipeline runs end-to-end, 21 tests
green, reports written.

## Notes

- Skills: `codebase-design` (deep modules, small interfaces) and `wayfinder`
  live in `.opencode/skill/`.
- Conda env `pyscf` (PySCF 2.13, Python 3.12). AGENTS.md is the canonical
  conventions/facts doc — update it when a new hard-won fact appears.
- Default plan: TDD at the seams (detbasis is pyscf-free and fully testable).
- This map is deliberately lightweight: the first slice fit one session.

## Decisions so far

- [Grilling: scope & configuration](tickets/01-grilling-decisions.md):
  cc-pVDZ, CAS(2,3) (2p shell only), SA over 3 triplet roots (Ms=1), operators
  validated through libcint integrals available inside PySCF.
- [Operator construction](tickets/02-operator-construction.md):
  first quantization — full matrices in the determinant basis via explicit
  Slater–Condon rules; no second quantization, no transition RDMs.
- [Angular momentum operators](tickets/03-angular-operators.md):
  libcint `int1e_irp` decode as primary construction; exact tesseral algebra
  as independent cross-check; both pinned by tests.
- [Active space](tickets/04-active-space.md):
  active orbitals = raw 2p shell AO columns (exactly rotation-closed), core
  folded with the spherically averaged density; CASSCF stages analyzed on
  their optimized orbitals (symmetry breaking is a result, not a bug).
- [Second quantization](tickets/05-second-quantization.md):
  second independent construction via explicit a+/a application
  (`secondquant.py`), wired as `build_operators(method=...)`; validated on
  Carbon in three tiers (elementwise, core-explicit CAS(6,5) vs PySCF
  ncore=0, fixed-core block vs folded+ecore).
- [Generalize pipeline](tickets/06-generalize-pipeline.md):
  `run_atom(cfg, shell_patterns, nelec_act, ...)` drives any atom; Ms
  sectors derived from spin + active-electron parity; sector tags label
  M_S = (na−nb)/2 ("0p5" for doublets); `run_carbon` is a thin wrapper;
  table-driven `tests/test_atoms.py` (new atoms = new data rows).
  Equivalence vs the carbon code proven with pinned ROHF input.
- [Si](tickets/07-si.md): 3p² works identically to 2p² through CASCI — first
  p-type core ([Ne] + 3s) folds exactly (textbook ³P/¹D/¹S, commutators
  ~1e-14). (Its "CASSCF/SA disabled" caveat was later retired — see
  ticket 09.)
- [Carbon ladder](tickets/08-carbon-ladder.md): all orbital-treatment rungs
  verified on carbon with an all-vs-all matrix. Raw core breaks symmetry
  measurably ({0}/{±1} split 4.4e-5 Eh, [H,Lz] 5.6e-6); averaged and
  fractional-occupation cores are exactly invariant (but differ by 6.6 mEh —
  different references); core relaxation gains −8.8 mEh at exact symmetry;
  free CASSCF breaking refined to m-incompleteness (window l-pure, fact 8);
  SA restores ⟨L²⟩ = 2.000000. Optimized stages fold with their own core.
- [Si analysis](tickets/09-si-analysis.md): the "³S multi-solution well"
  was a build_full_mo construction bug (non-orthonormal MO matrix from
  core↔window overlap, fact 13) — fixed block-safe (fact 14); Si now
  behaves like carbon end to end (CASSCF ³P −0.127 Eh, SA ⟨L²⟩ = 2.000000).
  Machinery generalized: run_ladder for any atom, virtmix rung, ΔE column,
  configuration labels, ortho_err divergence guard.
- [Sc](tickets/11-sc.md): first doublet, first TM. Two real bugs found and
  fixed (half-integer S labeling; frac-occ must select its block by AO
  projection, not energy order — the aufbau block runs away to 4p). The
  machinery itself held: CAS(3,6)/{3d,4s} is exactly closed (1e-13) and
  reproduces the full 90-state free-ion term table. TM phenomenology
  mapped (fact 15): percent-level raw-core polarization; free CASSCF
  collapses to an ²S-like basin (variational, not a bug); shell-start SA
  restores ⟨L²⟩ but not [H,Lz] (window 4s ↔ core 3p swap, m-incomplete);
  relax is the honest TM rung.
- [Review protocol](tickets/15-review-protocol.md): two-axis /code-review
  (fixed point f5cfb6f) caught the silently loosened ticket-10 merge bands
  (restored per-case) and added the missing Sc pins. Architecture survey
  (sub-agent + HTML report) confirmed three deepening candidates; user
  approved order 2→1→3, executed as ticket 14's pre-item commit:
  `stages.stage_facts` (one fold-and-analyze owner), `drivers.
  constrained_casscf(mo_full=, mask=, weights=)` over public
  `uniq_var_indices` masks, `ladder.Rung` registry (report iterates it —
  the drop-a-rung bug class is gone).
- [Constrained-optimization & start-variant axis](tickets/14-start-variants.md):
  the radial mask (same-(l, member) rotations over an (l,m)-labeled
  basis) answered all five questions — see facts 16/17. Headlines: the
  p² m-breaking is an energetically flat member GAUGE (radial rung ==
  free rung); the Sc escape is ~25 mEh below the radial ceiling (needs
  window-identity change); the complete {3d,4s,4p} window still escapes
  (optimizer attractor, not 4p hunger); sa_radial anchors the gauge but
  lands on the 1e-3 input-window baseline; the relax-rung ~1e-3 commutator
  is that same baseline, present before optimization (attribution
  resolved; motivates ticket 13). AVAS cross-validates the window; the
  guard negative control fires.

## Not yet specified

- Post-atoms follow-ups from the user ("next tasks" beyond Si/Sc/Ti) —
  expected to build on the operator algebra (candidates: spin–orbit
  coupling, transition moments, g-factors).
- Non-spherical/many-center generalization of the L operators (molecules):
  the current construction assumes a single atom at the origin.
- Cosmetic: the per-stage report line prints `Ms={na−nb}` (i.e. 2·M_S —
  "Ms=+2" for the carbon triplet sector). Found during ticket 06; left
  untouched for output stability, fix belongs in its own decision since it
  changes every report.

## Live tickets (updated 2026-10-24, ticket-27 wave)

- [25: ROHF fold-convention diagnostic](tickets/25-rohf-fold-convention.md)
  — RESOLVED + EXTENDED: the convention is machine-exact on all five
  atoms (1.8e-15…1.7e-14); the cross-check fully decomposed (the
  averaging effect per atom: 0.03/0.12/1.3/2.8/6.2 mEh); the retired
  raw-span ruler's window-definition offsets quantified
  (0/94.7/88.0/85.4/463.4 mEh).
- [26: elements tables T3a/T3b](tickets/26-elements-tables.md) — mostly
  DONE: T3a/T3b landed for all five atoms with verified values;
  README · Results = summary + pointer. Remainder: the study markers in the
  element files (needed by ticket 28's inserts), C/Si T3a rows, the
  per-atom study provenance notes.
- [27: launch & naming polish](tickets/27-launch-naming.md) —
  conv_accelerator DONE; sa-casscf naming DONE (want keys, treatment
  values, tables); ladder.py → treatments.py + sweep.py → study.py
  DONE (shims in place, tests renamed); the v4 diagram in README · The pipeline
  DONE. Remainder: none blocking — close at the wave review.
- [30: The project walk](tickets/30-project-walk.md) — OPEN (2026-10-28):
  the file-by-file review of every remaining src and test module per
  the AGENTS module-walk protocol (Refactor + Verification rules;
  facts homed and indexed as modules are walked). Phase 1 next:
  operators -> analysis. T2/M/U measurement runs stay GATED inside it.
- [29: Ti-first](tickets/29-ti-first.md) — OPEN, the agreed plan:
  T1 code-only (the build_full_mo column generalization IS the mo-pur
  second consumer; CaseSpec-driven key-set exported as spine/grid/
  diagnostics/planned groups, renderers consume; byte-diff guard),
  T2 the gated Ti batch (the 4 mo-pur cells, the (4,10)/(4,14)
  anchors, pass curves), then C/Si/Sc migration, then the AGENTS
  protocol round + the treatments/study retirement ticket (the
  test-pin migration is its deliverable). Supersedes the cellkey-wave
  queue ordering (archived: `archive/cellkey-wave.md`); deep
  refactor/generalize is mandated WHILE Ti work
  runs, with the second-consumer guardrail.
- [28: completion runs](tickets/28-completion-runs.md) — OPEN, the
  runner is ready: U SA-24 / relax-24 / SA-(6,13); Sc/Ti NO-start
  ladders; the approved stop rules (certificate / saturation / cost).
- [20: new-atom ergonomics](tickets/20-new-atom-ergonomics.md) —
  partially delivered (windows/core/reference/casscf split, cases.py,
  study.py, treatments.py, needle guard); the axes acceptance list in
  §6. 25–28 are its decomposition; close 20 when they land.

(The pre-2026-10-24 ticket list follows for history.)

Wave order decision (post-Sc): run the review and the constraint axis
BEFORE Ti -- Ti would only repeat Sc's phenomenology on unsettled
machinery; tickets 14/15 answer the open machinery questions first.

- [15: review protocol](tickets/15-review-protocol.md) — RESOLVED: step 1
  (two-axis review, findings fixed at 26cbd25), step 2 (survey, candidates
  2→1→3 approved and executed as ticket 14's pre-item), step 3 (protocol
  recorded in AGENTS.md Planning).
- [14: constrained-optimization & start-variant axis](tickets/14-start-variants.md) —
  RESOLVED: pre-item refactor (a024d26) + the experiment matrix; findings
  in facts 16/17 and the ticket's Resolution section.
- [16: deep understanding of C/Si/Sc](tickets/16-deep-understanding.md) —
  RESOLVED: plain docs (RUNGS/GLOSSARY), frozen radial scans (fact 18:
  full gain at exact symmetry on C/Si, exactly-flat member tilts, the
  Lz blind spots, the Sc JT-like inversion), SA-deviation decomposition
  + one-particle closure meter + AVAS verdict (fact 19: deviations are
  100% window; free SA heals below baseline on Si; span overlap ≠
  closure; AVAS inherits the mean-field relaxation — bad on Sc).
- [13: purified core / closed core span](tickets/13-purify-general.md) —
  RESOLVED (fact 20): `purified_core_orbitals` (shared radial per closed
  shell, E ⊗ members L-closed for ANY radial E — the ticket's own-shell
  premise was wrong, corrected in implementation); input windows close
  at machine precision (Si 3.4e-3 → 1e-15, Sc 2.1e-2 → 9e-15); THE
  DELIVERABLE: Sc CASSCF/SA-CASSCF with unbroken Lz and L² — relax_pur
  1.3e-13, sa_relax_pur 6.7e-14/3.6e-13 (90-state term table intact);
  free SA machine-exact on C/Si; Sc's residual sa_radial_pur gauge
  floor (~1e-6, no restoring force) and the surviving 4s↔3p swap
  attractor close the attribution questions. Phase 4 not needed.
- [12: Ti](tickets/12-ti.md) — RESOLVED (fact 21): options A+B in one
  session. A: d² tables textbook-exact on first pin (³F×7+³P×3; ¹G —
  the first L=4 term — in Ms=0). B: TI_EXACT_TERMS over 120 dets; the
  reference needed TWO fixes beyond the pre-flight — character-forced
  occupations AND an SCF-basin energy gate (occ_e_ref; a +11.2 mEh
  metastable member-gauge basin with a 1.6e-2 window defect exists even
  with forced occupations); the ³F tilt verdict: JT class (descends,
  ~½ Sc) so sa_radial_pur sits at its floor; a THIRD radial anatomy
  (raw rotations rise; the rungs' ~11 mEh gain lives on the
  orthogonalized substrate's s slot); relax_pur/sa_relax_pur
  MACHINE-EXACT (6.2e-14/1.4e-13) — the fact-20 deliverable transfers;
  the Sc basins transfer (escape +46 mEh, exact 4s↔3p swap on sa,
  partial swap on sa_pur). Suite 91 green.
- [17: drivers line-review](tickets/17-line-review.md) — RESOLVED
  (a3b09a1): _lowdin dedup x5, docstring diet (944 → 855; the −120
  estimate was optimistic, no content lost), _run unified, __all__
  fixed. Zero interface change; suite 92 green.
- [18: ECP pairing + p² validation](tickets/18-ecp-pair.md) — RESOLVED
  (code+tests; docs drafted for review): ccECP + ccecpccpvdz (the ONLY
  bundled PP family covering C/Si/Sc/Ti — cc-pVDZ-PP skips Si/Sc/Ti,
  bfd lacks Sc/Ti). C/Si = 4 valence e (ncore=1), Sc/Ti = 11/12 e
  (ncore=4, 3s3p explicit; ticket 19). MEASURED: term tables exact
  (³P/¹D/¹S, ⟨L²⟩ = 2/6/0, commutators ~1e-15 through SA); the
  ncore=0 full-valence fold is an exact no-op; raw detector LOUDER
  (1.1e-4 vs 6.1e-6); E_CASCI = E_ROHF = E_CASSCF and θ* = 0.00
  exactly — the all-electron "confinement" class (+250.8/+221.5 mEh)
  is a core-polarization artifact, not an intrinsic raw-window cost;
  spacings: one-parameter F₂ structure (¹S:¹D = 5:2 exactly in all
  four tables), ccECP shrinks F₂ by ~34%. Suite 96 green.

Ticket 19 (Sc/Ti ECP) — RESOLVED (tickets/19-ecp-tm.md): gated
references (Ti: good basin −57.60649802, defect 1.4e-8, default lands
the +11.1 mEh basin — fact-21a with the default REVERSED; Sc: forced
basin −46.12052897 from every guess, +0.30 mEh above the guess-dependent
unforced minimum — reproducibility chosen; unforced atom/huckel run away
to 4s²4p¹). Free-ion tables textbook-exact (Sc ²D×5, 90 dets, comm
1.5e-12; Ti ³F×7 over 120 dets, 5.7e-13). Ladder: escape on both
(−51.6/−57.3 mEh; Ti shell grabbed 4f, −0.88 Eh), Sc sa = the 4s↔3p
swap (fold illusion), Ti sa healthy this draw (bimodal — branch), pur
family MACHINE-EXACT at ncore=4 (relax_pur/sa_relax_pur 1.5e-12/5.3e-13;
closure certificate 4.9e-15/7.2e-15 vs 1.4e-2/9.5e-8; the pur fold kills
even the 1e-3 bookkeeping class). Test lesson: per-root labels via
analysis.analyze only — a wrong-axis einsum on exactly-degenerate
multiplets returns reproducible garbage (cost one diagnosis loop). Suite
99 green. ECP wave (18+19) COMPLETE. Next: U (u_14e) per the plan-mode
agreement — frozen-7s² CAS(4,12) tier first, sparse O·v mode deferred.

- [21: U — CASCI campaign](tickets/21-u-casci.md) — RESOLVED (both
  parts). Part 1: gated Newton reference (−50.38185233, 4/4 guesses);
  the label trap quantified (uncontracted labels = tight-end counters);
  three-strike canonical-order bug fixed at all defaults; fact 4
  triggered and handled (orthonormalized windows); pur mandatory
  ([H,L²] 2.2 → 5.1e-13); VECTORIZED + real-dtype slater_condon_matrix
  (identity-pinned, 18k dets hours→22 s; BLAS livelock diagnosed,
  128-thread rule); ⁵L° ×17 odd ground via the window series
  (1f: f¹d³ → 2f: f²d² → 3f: odd ground). Part 2 (2026-10-21):
  select_primitives (exponent bands from MEASURED open-MO weights);
  standard per-configuration term tables (u_run_part2.log, 10 sectors
  2h17m); CORRECTION: ⁵L° is configuration-PURE f³d¹ (⟨N_f⟩=3.000;
  part 1's 2.33 mix reading was an artifact — L=8 unreachable from
  f¹d³); the two-axis race: f-side stable (4f point, +339 mEh, f⁴ at
  +25 mEh), d-side FLIPS back to even ⁵L f²d² by 19.1 mEh (3f2d) —
  resolution = 4f+2d RESOLVED (⁵L° f³d¹ by > 37 mEh, −50.3920 BELOW the
  ROHF — the flip was an asymmetric-corner transient; ground term
  stable at every saturated corner); engine: PySCF Davidson
  straggles in degenerate clusters (38 mEh; pspace_size=ndet exact
  fix, fci_opts hook), dense LAPACK eig wall → ARPACK Lanczos
  (40,920 dets ~35 min/sector, validated vs exact).

- [22: window methodologies — AOs→NOs](tickets/22-window-methodologies.md)
  — OPEN (user frame, 2026-10-21): methodology-plural window
  construction compared on one axis (primitive bands / replicated SCF
  radials incl. virtuals / NO-occupation selection / hybrid NO
  compression; core/window split as a parallel axis). Hard constraints
  only: orthonormality + m-complete replication (fact 20 licenses any
  radial subspace; fact 19 forbids per-member thresholds).
- [24: matched-cost window comparison](tickets/24-matched-cost-windows.md)
  — OPEN (user frame, 2026-10-21): AO-exact (A) vs raw-ROHF-broken (B)
  at EQUAL det counts, U + Sc/Ti, contracted AND decontracted, short
  windows only (≤30 min/sector; heavy runs postponed by user rule);
  the Sc plateau-vs-escape question (fact 16b's +25 mEh must melt in a
  wide-enough m-complete window — full-space intuition); delivers the
  README · Results (window lattice) section (table + per-block ΔE matrix with
  explicit references + the representability comment, brief).
- [23: U SA wave](tickets/23-u-sa-wave.md) — RESOLVED 2026-10-24
  (README §6.4 ticket-23 block): SA over the 17×⁵L block from NO-12 =
  −50.4238 at exactly quantized ⁵L — 42 mEh below the ⁵S° ROHF
  reference, the first U object below the mean field at certified
  symmetry. Mechanism: RADIAL relaxation (l-content preserved, weight
  moves across exponent bands); only the ensemble holds the gauge
  (free single-root destroyed by the ⁵S° attractor; SA-17 m-incomplete).
  7s²-as-valence inert (−0.58 mEh); inner core pairs = the next tier
  (accidental 6p promotion: −551 mEh). Exact-eigh labels mandatory —
  PySCF ci straggle mislabels ⟨L²⟩ in the 17-fold cluster (fact
  candidate for the next review).
  Order: 24 (light) → 23 → 22 → review wave (two-axis +
  codebase-design over 4f3cf27...HEAD) → architecture discussion →
  ticket 20 (the axes-product ergonomics) last.

The atom campaign (C, Si, Sc, Ti) is complete in BOTH reference frames
(all-electron + ccECP); U part 1 done. Beyond: SA wave, ticket 20
(ergonomics), then user candidates: PHYSICS.md narrative (deferred),
spin–orbit coupling, transition moments, g-factors.

## Post-wave survey (2026-10-19, deep-module pass + code review)

Reviews run: two-axis code review of `96e3c61..HEAD` (0 hard standards
violations; spec verified — fixes applied at d30f366: TI_OCC
single-sourced, Lz-spectrum + ms0-ordering pins landed, bands aligned)
and the deep-module survey. Survey verdicts: all modules deep or
correctly arranged; tests cross public seams only. Findings carried
forward:

- **ECP design (pre-verified)**: the seam already exists — PySCF's
  `get_hcore()` includes `int1e_ecp` when `mol.ecp` is set, and h1 is
  built in exactly two places (`operators.active_space_integrals`,
  `operators.full_space_integrals`). The ECP ticket = `AtomConfig.ecp`
  field + `build_mol` wiring + verification. TWO WRINKLES: (a)
  `pipeline._ncore` recomputes nelectron from nuclear charge — WRONG
  under ECP (core electrons removed); use `mol.nelectron`; (b) a named
  `h1_ao(mol, mf)` helper in operators.py would make the seam explicit.
- **Deepening candidates** (approved order): DONE (post-switch session):
  `ReferenceSpec` dataclass owns the occ/frac clump + its validation,
  threaded as `ref=` through run_rohf/run_ladder/run_atom (conftest TI_REF;
  pipeline._ncore deleted — `system.ncore_of(mol.nelectron, …)` and
  `system.sector_occupations(…)` are the shared sector geometry, and the
  mol.nelectron switch pre-closes the ECP (a) wrinkle). Deferred: none
  blocking. DONE earlier:
  `StageFacts.analysis` typed + the three read fields exposed as
  properties (the one interface leak).
- secondquant.py survives the deletion test as a *witness* (independent
  construction), not leverage — keep, documented as such.

Docs survey (pending, first item after context switch): AGENTS physics
facts (18–21) compress to agent-facing content + README pointers;
GLOSSARY.md folds into a README appendix; wayfinder/tickets stay.
[2026-10-19 update: DONE — commits 35cd327/cf369f3/2e0f5ad; deepening
pre-item DONE at 9a22ab0 (ReferenceSpec + shared sector geometry, the
_ncore-from-Z ECP trap closed); test tiers fast/slow at 345351a.]

NEXT WAVE: ticket 17 (drivers line-review: docstring diet ~−120 lines,
Löwdin dedup ×4, __all__ fix, _run branch unify) then the ECP ticket
(design pre-verified above; wrinkle (a) already closed by ncore_of).

[2026-10-19 late update: docs survey + deepening + test tiers DONE
(35cd327…345351a); 7.3 reference axis measured and card corrections
landed (f38fb6d) — including the falsified SA-flatness draft, replaced
by measured tilt-mean curves (C/Si descend, Sc/Ti rise — §3.4); the
handoff is in message.txt.]

[2026-10-25: README §1 (Goal and scope) DELETED — its content is not
redundant only in three pieces, now placed: the "certified" definition
→ §1.1; the three-part deliverable (terms / comparisons / loss
accounting) + the non-relativistic scope → Introduction's closing
paragraph; the per-atom-data-not-code claim → the §3 launch path.
README sections renumbered −1 (§2→§1 … §9→§8); live pointers updated
in AGENTS.md, RUNGS.md, elements/*, and this map. RESOLVED tickets keep
their original numbering — their §refs are historical.]

## Out of scope

- Performance tuning of the determinant-basis builder beyond CAS sizes used
  here (explicit Slater–Condon is O(ndet²); fine for ndet ≤ a few thousand).
- Cartesian basis support (angular analysis requires pure spherical functions).

## Post-wave survey (2026-10-23, tickets 18/19/21/24 + ticket-22 round 1)

Two-axis review of `4f3cf27...HEAD` (22 commits). Fixed now (commit
pair with this entry): (1) `_L_LETTERS` reused P at L = 12 — the
spectroscopic sequence skips J then the second P/S, so O(11) is
followed by Q(12); table extended to L <= 20 with an explicit
beyond-table ValueError, pinned in test_analysis. (2)
`occupied_core_idx` — the fold's core-selection semantics — had no
test (the one hard standards finding); the U-style scattered-block,
contracted-coincidence, and too-few-closed branches are now pinned.
(3) `_slater_condon_reference` removed from detbasis `__all__`
(private name; tests import it directly). Recorded, not fixed:

- **Duplicated `closed_shell_core_energy` / `_ti_ref`** across
  run_u_no.py / run_tm_no.py / run_tm_methb.py (scan-stage scripts;
  library-promotion deferred until a second consumer needs the core
  energy — candidate for the ticket-20 axes pass).
- **Lazy `from atomterm.drivers import ...` cycle dodges** repeated in
  4+ functions across operators/drivers (layering wrinkle; a shared
  core-resolution seam was considered and declined this wave — the
  call sites are settled and the churn risk outweighs the gain).
- **`fci_opts: dict` sets arbitrary fcisolver attributes by string
  key** (primitive obsession; a named `pspace_size` kwarg would be
  type-checked). Declined: the dict is the documented extension seam
  in run_casci's docstring and callers pass exactly one key.
- **`slater_condon_matrix` dtype probe checks one operand**; a mixed
  real-h1/complex-g2 call would raise a casting TypeError rather than
  mis-execute (fail-loud, no in-repo caller). Accepted as-is.
- **Ticket-22 open items** (tracked in the ticket): methodology B
  proper (replicated SCF radials) — decide after this wave; scf_mo
  selector has no certificate run (only the slicing pin); the U
  f3d1-vs-f2d2 ordering comparison inside one NO window (the pilot
  pinned the 5L ground + the next block instead).
- **Ticket-20 scope creep, self-declared**: ergonomics fields
  (AtomConfig.tag, raw basis/ECP passthrough, occupied_core_idx)
  landed during the write-only ticket-20 — the remaining ticket-20
  items (axes product, --fresh, intent manifest) still open.
- **Suite flake census (2026-10-24, ticket-23 wave)**: two full-suite
  runs, two different single-test failures, both in the documented
  chaos classes and neither reproducing file-scoped: (a) Sc sa_pur
  swap-basin commutator fresh-drew 0.046 vs the `> 0.1` pin — band
  widened to `> 0.03` per the 2026-10-15 precedent (p-weight stays
  the discriminator); (b) one session threw ValueError ERRORs in
  purified/radial rungs (C/Si) — not reproduced, traceback lost to a
  `tail -3`; consistent with the needle/Löwdin sensitivity whose fix
  is ticket-20's needle linear-dependence helper. Suite hygiene rule:
  never pipe pytest through `tail` — full log to a file.

Design survey (2026-10-23, wave delta over src/): natural.py is DEEP
at a real seam (four selector adapters; the 1-RDM/sign/transition-table
machinery behind a 4-arg interface; the m-scalar invariant lives in the
module docstring, the certificate is the enforcement tool).
occupied_core_idx deepens the fold seam (callers hand mf+ncore; the
closed-block selection + invariant check are inside). No deletions
warranted; prior 2026-10-19 verdicts stand for untouched modules.
Forward note for ticket 20: Window's dual representation (ao_rows |
c_win) is resolved statically by today's callers; the axes-product
runner will need ONE resolution method (Window.columns(mol, mf)
routing ao_rows through build_full_mo) so selectors stay swappable
there.

UPDATE (2026-10-24, docs wave): README refactored per user review
(compression + generalization; method sections carry structure, values
live in elements/); elements/<Atom>/ = TABLES.md + NOTES.md per atom;
README 7 = table-kind enumerator with Ti as the worked example; 3.6
dissolved (deformations math -> 2.4; SA logic -> 3.4); 5 flattened to
table + per-level comments; 8/9/appendix compressed (README 1420 -> 996
lines). Ticket 20 carries the axes-systematization acceptance list.

## EXECUTION LEDGER (2026-10-24, the durable copy — write-off of the
architecture discussion; do not lose)

README: (1) 2.2 gains the fold formula + the orthogonality requirement
block; 3.2 references it; 3.3 gains the leak-guard line + the
orthogonalization-routines enumeration (Lowdin=shape-preserving-core,
S-GS=identity-preserving-window/radial, completion=span-filling);
(2) 7.3 mechanism tables -> T3b (single E as-run + cross-check);
(3) 4: the v4 launch diagram; sa-casscf naming sweep.
CODE: windows.radial_basis strong-only admission
(STRONG_RESIDUAL_WEIGHT = 1e-3, fact 16e gap), demotion helper deleted,
complete_orthonormal docstring (isometry derivation, the {0,1} count,
the flake note); sweep.py -> study.py (StudySpec/plan_study/run_study);
ladder.py -> treatments.py (registry; shim).
TESTS: Sc/Ti ladder re-run bands; planted-needle toy test; study
plan/render.
DOCS: AGENTS conventions += docstring standard (@brief/math/@param/
@note/@example) + naming rule; RUNGS.md staleness check.
THEN: ticket 28 completion runs; two-axis review + AGENTS facts.

## Post-wave survey (2026-10-25, README review wave + explicit ladders)

Two-axis review of the uncommitted working tree vs 0eba2c5 (README
step-by-step review + Sz delivery + explicit ladder operators; no
ticket — spec = the conversation directives). Fixed immediately:
stale test counts (README 91/116, AGENTS 123 → 128), AGENTS Layout
rows for angular/operators (+l_ladder_matrices, +Sz +ladders), the
Conventions exception for the inter-sector a†/a ladder construction,
the surviving L̂ᵔ/Lᵔ glyphs (3×, incl. two my earlier pass missed —
hatless variant), §6.3 heading "ladders" → "rungs", stale
results/Ti_ccpvdz_s2_report.{md,json} regenerated (Sz line + the
Ms=(na−nb)/2 sector-label fix — provenance for the README block).
Recorded, not fixed:

- **`spin_ladder_matrix(..., plus=True)` flag arg** — two operators
  behind one boolean; docstring covers it; split only if a third
  consumer appears.
- **`analyze(h, lz, l2, s2, sector, sz=)` data clump** — stages.py
  unpacks ops[...] field-by-field; passing the bundle is the deeper
  fix (analyze signature churn across 5 call sites; defer).
- **"ladder" vocabulary residue** — glossary keeps "ladder / rung"
  (deliberate: documents the historic name) and artifact names
  `*_ladder.md` + `run_ladder` remain until a naming ticket lands.
- **Full suite before commit** — deferred at the user's request
  (busy threads); fast tier green (56). operators/stages/analysis/
  cases all changed, so the run is mandatory at commit time.

[2026-10-25 Wave 1: README restructured on the Acts narrative (plan:
.opencode/plan/readme-acts-restructure.md; uncommitted, awaiting the
user's paragraph-by-paragraph review — THE spine of this wave). §1
language (1.1 symmetries/certified, 1.2 operators+ladders, 1.3 window
machinery, 1.4 classes+meters); §2 physics NEW (2.1 restriction+fold,
2.2 stabilizer/Λ-branching with carbon as the marked instance, 2.3
causal chain + L-dies-twice with the SU(2) label protection +
cut-feedback-loop + polarization corollary, 2.4 optimizer side +
scan verdicts folded in, 2.5 levels/resolutions); §3 methods no-math;
§4 pipeline diagram; §6 Ti-only (scans section deleted — verdicts in
§2.4; cross-atom tables → element files); §7 module table moved in;
§8 verification program. External pointers remapped (AGENTS facts,
RUNGS, elements, docstrings). Waves 2 (TABLES.md semantic renaming)
and 3 (verifications) wait for the review.]

[2026-10-26: README restructure LANDED (the Acts narrative:
Motivation/Methods/Results, name-based refs, unnumbered) + Wave 2
(TABLES.md semantic renames; pointer sweep; basin map Stage A: C free,
5 starts -> ONE energy, five tilted gauges) + the gated Ti ladder
regenerated in place (elements/Ti/ladder.json) with the value audit:
the SA-free row's E was our fold's partition-inconsistent reading
(-849.27); the honest as-run total is -848.4166 (dE -10.1, xchk 851 =
the partition break). CASCI-above-ROHF = representability, normal.
The 2.77 mEh cross-check = the raw fold's Lambda branching (solved).
Real AOs are ±|m| combinations — M_L bookkeeping needs the
real->complex transform. HANDOFF: message.txt (task list: multiplet
M_L labels, README 4-axis Methods, frozen-window naming sweep, series
layout d2 2,5/8/13 + 3d2s2 4,6/9/14, basin Stages B/C, the 0.9 mEh
trace, full suite). Next-wave candidates: the linear-molecule
crystal-field test (the Lambda-branching cross-check).]


## Post-wave survey (2026-10-27, README review round + TABLES repair + survey fixes)

Fixed now (commits da4e08c..HEAD of this wave): the round-8 TABLES
defects (folds header + anchor-out, stale duplicate tail, planned-grid
start labels, C/S legends); the 27-point README round (starts/windows
unified — window = span, start = the columns defining it, one
definition; injection-point failure inventory; treatment-grid block
with spine/diagnostics/redundant lists + gates + direct-runs; lattice
vs treatment-grid split; <Lz^2> block -> operator row + diagnostics
row; epsilon single home; NOT-CONVERGED verdict token, plumbed
converged -> StageFacts -> payload -> rungs report, study-consistent
predicate); purified_s_core DELETED (retargeted test_secondquant to
purified_core_orbitals — the l=0 "reduction" is span-exact but
column-gauge-free (fact 16e); the span pin now scales its bound with
the measured contamination); StageFacts duplicate ground_term/ground_l2
properties deleted (later defs silently overrode — Duplicated Code);
cases.py stale workdir="results" defaults -> "elements" (results/ died
in round 1); no_group_windows MOVED study.py -> natural.py (the NO
window factory belongs with the NO mathematics; lazy Window/
build_operators imports — windows<->natural is a module-level cycle).

Recorded (candidates + declinations):

1. CELL-KEY SEAM (the wave's main design item, user-proposed): a
   (window, method, start, mask, substrate, fold) key as the single
   name of every run — payload field, chk store name, table row
   address. StudySpec.cell_key(wtag, treatment) is the germ; promote
   the grid coordinates into the declaration layer and make the runner
   callable per cell. Subsumes the spine/run-order runner (the
   user's run-queue generalization) and the ticket-20 intent manifest.
   Design it DESIGN-IT-TWICE style before the U(14) wave.
   ==> RESOLVED 2026-10-27 (b) as the CellKey wave, plan
   `archive/cellkey-wave.md`: cells.py landed P0-P2 (commits
   5f17df2 RUNGS.md dissolved, 09bdc28 key+bridge+load, 2a6132e
   run_cell/analyze_cell; 145 tests green). The bridge found a REAL
   bug on the way (fact 24: PySCF FCI ordering != detbasis ordering
   for norb>=4 — configuration labels and study diagnostics were
   scrambled on Ti/Sc/U windows). REMAINING: P3 the gated queue
   through run_cell + cells.{json,md} rename; P4 the strangler
   retirement (Rung registry tags -> axis values, StudySpec over
   key-sets, orphaned U scripts). start="no" waits for the
   parent-prefix policy (which NO group a nominal label selects).
   ==> P3 ALSO RESOLVED 2026-10-27 (b) (6f22dfe + 63adf28): the
   no-start policy landed (window grammar <source>no<k>); regen_cells
   regenerates C/Si/Sc/Ti through run_cell (Sc gained a
   character-forced ReferenceSpec after a broken-basin draw);
   build_methods_tables renders the Ti manifest from cells.json (the
   table/provenance cannot diverge; unified method-column tables per
   the user's spec). Two run bugs found+fixed en route: the
   budgets->weights leak (CASSCF cells ran as SA) and the ci
   object-array serialization (ci_stack now). Suite 148 green.
   REMAINING: P4 only; plus the campaign continuation through the
   seam (lattice cells incl. NO starts, optimizer manifest cells,
   pass curves) -- gated per run as before.
   ==> QUEUE STEPS 1+2 EXECUTED 2026-10-27 (c), commits aa21071..1d51de7:
   the d2 lattice (9 cells, CAS(2,5/8/13) x AOs/MO/NO; section rendered
   from cells.json by scripts/build_lattice_tables.py; the seeded anchor
   row was wrong on both counts: 10 dets not 50, -848.3973 measured not
   the (4,6) campaign number; {3d,4p,4d} crosses below the anchor) and
   the FULL Ti manifest grid (33/33 cells + verdicts: the escape
   attractor is start-invariant; the frozen-MOs rows are a new
   'defective window' class with the energy-vs-certificate inversion).
   THREE more seam bugs found by the queue, fixed+pinned: the
   occupation-grouping m-set split (fact 23 amendment: scale-aware
   tolerance; Ti d-shell leakage 2.7e-7), the SA store completion
   (store_name=; the roots store is the loader's preferred one), and
   the mask wiring (frozen/virtmix are start-independent; radial is
   ao-only validated at the key -- mo/no-masked cells had silently run
   FREE). NO-window selection is by group LABEL, not compact-list
   index (k counts groups; sector-small prefixes are skipped). Suite
   150 green. REMAINING QUEUE: C/Si/Sc grid completion (their cell_set
   now declares the grid), the (4,10)/(4,14) planned windows, pass
   curves, P4.
2. stage_io loader (user-proposed): the atomterm/<name>/{ci,energies}
   chk groups are WRITE-ONLY today; ci is PySCF FCI-ordered, not
   detbasis-ordered. A small loader + ordering bridge unlocks post-hoc
   analysis (RDMs, properties, transition moments, g-factors) from
   stored runs without re-running optimizers. Prerequisite for U(14)
   property work.
   ==> RESOLVED 2026-10-27 (b): absorbed by entry 1 (load_cell/
   analyze_cell + pyscf_ci_permutation; the write-only gap is closed
   cell-addressed; old rung chks are not cell-addressed and stay
   provenance-only).
3. LADDER_CASES fixture rename in test_treatments.py (14 sites; the
   rungs rename is mid-flight there).
4. TABLES.md repeats the reading key in four sections — data
   duplication is fine, prose duplication is not; one key + pointers.
5. scripts/ audit: keep README-cited anchors (diag_fold, regen_rungs,
   multiplet_split); retire/fold the orphaned U-campaign one-offs
   (diag_s5, diag_mdecomp, run_u_* adapters) once the cell-key runner
   exists.
   [SUPERSEDED 2026-10-30, ticket 30 phase 6: diag_fold and
   multiplet_split RETIRED with the debug-residue + one-shot waves --
   their findings institutionalized (facts 25/24, committed tables),
   their cites now read 'retired -- git history'; run_u* kept pending
   the ECP-wave U fork (declare U in CAMPAIGN). The audit's 'keep'
   list presumed the scripts were the only anchor owners; the pins and
   tables are.]
6. User-authored README sentences flagged by review, deferred to the
   user's next pass: "The certified variants are similar to CASSCF
   theory." (garbled; the pur variant no longer named), the three
   "Carbon atom" hedges, the deleted free-ion term-table paragraph.
7. Declined: table-structure regression tests (the safe-rebuild
   assertions are process, not product; the markdown is regenerated,
   not coded); making _report public (the fast test pinning it is
   acceptable; revisit if a second renderer appears).

Suite state: 148 green (was 130 at the survey; +18 = test_cells.py.
Two earlier full-suite failures —
test_atom_term_tables[Sc], test_ti_free_casscf_escapes_the_shell —
did not reproduce in isolation nor in rerun; recorded as fresh-draw
basin flakes, facts 10/15g; no band change: no observed values to
cite). FLAKE MECHANISM IDENTIFIED 2026-10-27 (b): the rungs fixtures
are SESSION-scoped — one chaotic SCF draw per suite run fails
whichever tight band grazes it (rotating failures across runs: Sc
tables, the ladders fixture, Sc complete-window, Ti sa_pur 2.2e-2
captured+widened, Ti casci closure, C virtmix uncaptured transient).
Protocol recorded in AGENTS · Tests.

## Post-wave survey (2026-10-28, CellKey wave two-axis review + design round)

Review span: d8cd792..df6ec4f (the CellKey wave P0-P3 + queue steps
a/b + the reading passes), two-axis per protocol (Standards/Spec
sub-agents over AGENTS.md + the plan file + the handoff specs).

Fixed now:
1. (Spec, the substantive one) The substrate axis was a silent no-op
   for start="mo": run_cell's mo branch sets mo_full = mf.mo_coeff
   regardless of substrate (documented "moot" in the branch), yet
   cell_set declared four (mo,pur) optimizer cells -- the same
   experiment under two addresses (energies/window digit-identical;
   the meter differences between the twins were pure rerun chaos --
   two independent runs of one cell). Fix: the key REJECTS mo x pur
   (void, like radial x mo); cell_set + the manifest renderer pruned;
   the four stale chks deleted; Ti regenerated cache-hit (29/29
   manifest cells; the diff is PURE DELETIONS -- the analyze path
   reproduced everything else byte-stable). Door left open: a
   pur-core MO start is a real experiment needing the
   window-projection machinery (one new builder + relaxed validator)
   -- not built, not wanted yet.
2. (Spec) Two of the three seam-bug fixes were unpinned. Pinned on
   carbon: test_sa_run_completes_the_roots_store (the store_name=
   WRITER side -- completing <address> instead of <address>_roots
   pairs ensemble energies with roots CI and the preferred-store
   load fails loud) and test_frozen_mask_applies_to_the_mo_start
   (mask start-independence: the converged window columns == the
   canonical MO slice up to phase; a free run rotates them away).
3. (Standards, hard) AGENTS analysis.py row missing shell_occupations;
   the store-contract docs (cells.py docstring + fact 24) named `ci`
   where the key is `ci_stack`. Both fixed.
4. (Standards, judgement) regen's render parsed addresses positionally
   (ruler.split("_")[2], twice) -- now through CellKey.from_address.

Declined:
5. AUDIT hardcodes address strings (Primitive Obsession flag): the
   AUDIT dict is human verdict DATA; addresses are its natural keys.
   The renderer's key() helper covers constructed addresses.
6. REFS re-states the conftest TI_REF numbers: known-declined class
   (2026-10-23 precedent) -- the runner stays self-contained
   gate-wise; accepted duplication.

Candidates (recorded, not now):
7. THE manifest grid is declared twice: regen_cells.cell_set vs
   build_methods_tables' loops + CASCI_FOLDS. Run-side drift renders
   as pending (caught); policy drift between the two stays silent.
   The fix -- export the grid declaration beside spine_keys in
   cells.py -- belongs to P4 (StudySpec-over-key-sets gives it a
   natural home). Flagged again by the design round below.
8. run_cell's nroots=min(9, ndet) magic 9 (the term-table root cap):
   name it when cells.py is next touched.
9. verbose is both a run_cell kwarg and a budgets key -- two knobs
   for one thing; collapse at P4.

Suite state: 153 collected (+2: the SA store-writer pin, the
frozen-mask wiring pin); one rotating-draw flake in the full
run (test_ti_casci_exact_closure_and_term_table, isolation-
green on rerun -- the session-draw class, AGENTS Tests; no
band change, no observed value captured).
Design round (same date, the codebase-design vocabulary over src/):

1. FIXED NOW -- THE INWARD FACADE: operators/stages/analysis/treatments
   imported through drivers -- their own compat facade -- keeping it
   mid-cycle (drivers -> windows -> natural -> operators -> drivers,
   held together by lazy imports). Owners import owners now
   (windows/core/reference/casscf directly, lazy sites stay lazy);
   drivers is OUTWARD-ONLY (scripts/tests), as its docstring always
   claimed. Deletion test on the inward usage: the complexity never
   lived in the facade -- pure misdirection.
2. THE CAMPAIGN DECLARATION HAS TWO OWNERS (sharpens P4 + review
   candidate 7): cases.CaseSpec declares the per-atom campaign;
   regen_cells re-declares it (per-atom window/sector + REFS), and
   the manifest grid is declared twice more (cell_set vs the
   renderer's loops). The product declaration belongs behind ONE
   interface -- a declared key-set exported beside spine_keys,
   consuming CaseSpec's references; runners and renderers CONSUME
   it. P4 = "one declaration, many renderers", not just "fold the
   scripts".
3. treatments.py FAILS THE DELETION TEST POST-P4, by design: once
   cell keys carry the axis values, deleting it makes the rung loop
   reappear as cells' key-set iteration (which already exists).
   study.py likewise becomes a render contract over key-sets. Both
   already the plan; the survey confirms the direction.
4. windows.py carries the widest interface in src/ (22 defs): window
   naming/selection vs substrate construction (radial_basis/
   lm_sectors/lowdin) is a visible internal seam -- but no second
   adapter justifies splitting today (one adapter = hypothetical
   seam). DECLINED; revisit if a second consumer of the substrate
   machinery appears (the pur-core MO builder would be one).
5. cells.py is the wave's deep module done right: a 5-name external
   seam, the store format hidden, axis resolution through owners,
   and the test surface == the interface (run == direct, analyze ==
   in-session, store round-trip). The budgets dict (five knob names)
   is the one interface leak -- documented, acceptable. stages.py is
   exemplary depth (2 public names).

## Post-wave survey (2026-10-30, ticket-30 phase 6 — code-review + design survey)

Review (fixed point 30ab209, the phase-6 scripts+tests commits): the
R6 lockstep misses fixed in-commit (diag_fold/basin_map dangling
cites, stale README counts, the map item-5 supersession, the sys.path
bootstrap idiom unified — radial_scan had NO bootstrap: its own usage
line failed). The suite caught a REAL straggler the fast tier cannot
see: test_atoms' `results["sa_casscf"]` (underscore) left by the
0909b64 sa-casscf rename — KeyError, deterministic, green again in
isolation. LESSON: a rename round that touches stage tags owes a FULL
suite run, not the fast tier (slow-only lookups are invisible to -m
fast).

Design survey findings (vocabulary: module/interface/depth/seam/
adapter, the deletion test):

1. **The window-kind registry is a HYPOTHETICAL SEAM** (verified):
   `windows.SELECTORS` + `window_shell`/`window_primitive`/
   `window_scf_mo`/`window_no` have ZERO production consumers — only
   tests/unit/test_natural.py; production resolves windows through
   cells.run_cell's start-axis branches and the select_* family. The
   name->window dispatch the registry was meant to own REAPPEARED in
   run_cell's branches instead. DECISION NEEDED (phase 7 or the ECP
   wave): retire the constructors (tests re-home to natural's own
   builders) or consummate (route run_cell's start axis through the
   registry — needs substrate/pur handling inside, likely forced).
   Candidate declined for now: recorded, not executed.
2. `analysis.save_report` RETIRED (zero consumers, the deletion test);
   `run_case`'s missing `return res` FIXED (declared `-> dict` lied —
   the phase-5 study-arm deletion took the return with it but left
   the assignment).
3. `casscf.run_casscf` (one production caller, cases) — kept: the
   stage-name contract; the phase-5 verdict stands (borderline,
   revisit only if cases' stage list changes).
4. `stage_facts`'s `c_act` parameter: every production caller passes
   the same derived slice — deepening candidate (derive inside,
   override for scans), medium confidence, deferred.
5. The rest of src/ reads deep-and-right (cells.CellKey, the mask
   seam with three live adapters, operators' twin PySCF boundaries,
   stages' 2-name surface called exemplary by the 10-28 round).
