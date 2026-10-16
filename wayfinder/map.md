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
- [12: Ti](tickets/12-ti.md) — NEXT: d² tables (A: CAS(2,5)/{3d}, first
  L=3 ground term and first L=4 term, first f-block load) + the Ti-atom
  ladder (B: CAS(4,6)/{3d,4s}, pur rungs, the ³F tilt-stability
  question). Pre-flight verified 2026-10-16: default ROHF is an aufbau
  runaway (+123 mEh) — the ticket ships a verified character-forced
  occupation recipe.

## Out of scope

- Performance tuning of the determinant-basis builder beyond CAS sizes used
  here (explicit Slater–Condon is O(ndet²); fine for ndet ≤ a few thousand).
- Cartesian basis support (angular analysis requires pure spherical functions).
