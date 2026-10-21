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

## Live tickets (the current wave, reordered after ticket 11)

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
- [23: U SA wave](tickets/23-u-sa-wave.md) — OPEN: 17-root SA over
  the ⁵L° M_L set; relax/pur-relax; the open question — orbital
  relaxation vs window growth (ticket 22's comparison in disguise).
  Order: 21p2 → 22/23 → review wave (two-axis + codebase-design over
  4f3cf27...HEAD) → architecture discussion → ticket 20 (the axes-
  product ergonomics) last.

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

## Out of scope

- Performance tuning of the determinant-basis builder beyond CAS sizes used
  here (explicit Slater–Condon is O(ndet²); fine for ndet ≤ a few thousand).
- Cartesian basis support (angular analysis requires pure spherical functions).
