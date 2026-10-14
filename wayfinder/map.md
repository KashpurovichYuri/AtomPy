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

## Live tickets (the current wave)

- [10: Sc](tickets/10-sc.md) — 3d¹4s² ²D; first doublet, first TM, ROHF
  occupation risks documented. Full ladder via the ticket-09 machinery.
- [11: Ti](tickets/11-ti.md) — 3d² ³F; L = 3/4 terms, core = [Ar]+4s²
  (apply 10's learnings).
- [12: purified_core_orbitals](tickets/12-purify-general.md) — optional,
  per-shell core purification for future core-explicit TM tests.
- [13: review protocol](tickets/13-review-protocol.md) — /code-review after
  the wave (fixed point 507b796), then /improve-codebase-architecture.

## Out of scope

- Performance tuning of the determinant-basis builder beyond CAS sizes used
  here (explicit Slater–Condon is O(ndet²); fine for ndet ≤ a few thousand).
- Cartesian basis support (angular analysis requires pure spherical functions).
