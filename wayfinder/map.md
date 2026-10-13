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

- [07: Si](tickets/07-si.md) — 3p² ³P; first p-type core (data row only).
- [08: Sc](tickets/08-sc.md) — 3d¹4s² ²D; first doublet, first TM, ROHF
  occupation risks documented.
- [09: Ti](tickets/09-ti.md) — 3d² ³F; L = 3/4 terms, core = [Ar]+4s²
  (apply 08's learnings).
- [10: purified_core_orbitals](tickets/10-purify-general.md) — optional,
  per-shell core purification for future core-explicit TM tests.
- [11: review protocol](tickets/11-review-protocol.md) — /code-review after
  the wave (fixed point 507b796), then /improve-codebase-architecture.

## Out of scope

- Performance tuning of the determinant-basis builder beyond CAS sizes used
  here (explicit Slater–Condon is O(ndet²); fine for ndet ≤ a few thousand).
- Cartesian basis support (angular analysis requires pure spherical functions).
