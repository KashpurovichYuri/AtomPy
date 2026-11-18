# SPEC — The perturbation analyzation (`perturb.py`)

The wave's contract. Tickets 42–46 execute it rung by rung; each
ticket references this file and carries its own slice.

## Problem Statement

The campaign's tables certify *static* correlation only: a CAS span,
its exact symmetry, its commutators. But the deeper question is the
**non-saturation of CASSCF itself**: enlarging the active window keeps
buying large, non-converging amounts of energy (the poverty class:
+250.8 mEh on C, +221.5 on Si between the compact and @r1 windows, at
identical determinant counts and exact certificates). Two readings
compete: the window is missing genuinely STATIC correlation
(near-degenerate radial directions that must be variational), or the
CAS is absorbing DYNAMIC correlation — badly, at growing cost, never
saturating. Only perturbation theory can tell: PT2 measures what the
external space supplies on top of each window. If the PT2-corrected
ladder saturates where the CAS ladder does not, the "missing" energy
was dynamic correlation in fact; if the total keeps moving, the
missing part is static and belongs in the window.

## Solution

One new deep module, `perturb.py`, in the analysis row: second-order
perturbation theory over a certified CAS state, built as a ladder of
H₀ choices from simplest up (Epstein–Nesbet → Dyall → NEVPT2 strongly
contracted; CASPT2 planned, not built). First atoms: **C and Si**.
The analyzation is the **window-ladder diagnostic**: apply PT to the
different CAS windows the campaign already stores, per declared
solution (root, term-labeled), and read whether dynamic correlation
mixes with the static one — quantified by how much of the CAS spread
the corrected ladder flattens. A dedicated validation script
certifies correctness BEFORE any table renders; results land in a
per-atom `PERTURBATION_TABLE.md`, separate from TABLE.md. No
campaign/CellKey integration — post-processing only, wave 1.

## User Stories

1. As the researcher, I want to declare which solutions (stored CI
   roots) to perturb, so that the multi-solution class stays under my
   control and no draw is silently chosen.
2. As the researcher, I want each rendered row term-labeled (from the
   analysis readout), so that a perturbed root's identity is visible
   without cross-referencing.
3. As the researcher, I want a separate validation script that
   certifies the PT2 constructions accurately (twins, analytic
   identities, cardinalities) at realistic sizes, so that no table
   ever renders from unvalidated math.
4. As the researcher, I want Epstein–Nesbet first, so that the
   universal PT2 baseline exists and every fancier H₀ must face it.
5. As the researcher, I want the Dyall H₀ and NEVPT2 strongly
   contracted classes next, so that a parameter-free, intruder-free
   PT2 is available (no level-shift dials).
6. As the researcher, I want per-atom `PERTURBATION_TABLE.md` files,
   so that the dynamic-correlation axis reads beside TABLE.md without
   overloading it.
7. As the maintainer, I want the module pure (no optimizer, no
   chkfiles, no draws), so that reproducibility is the whole
   computation and the trusted-cache conventions stay untouched.
8. As the maintainer, I want the module's import edges downward only,
   so that the architecture stays cycle-free.
9. As the tester, I want the frozen carbon fixture to carry the unit
   tier, so that no SCF runs in tests and pins are seconds-class.
10. As the tester, I want degenerate anchors (CAS = whole space → 0;
    one active electron → the closed answer; empty active → the MP2
    limit under Dyall H₀), so that any deviation is a bug, not physics.
11. As the tester, I want two independent constructions of the
    Epstein–Nesbet correction (sum formula vs first-order solve), so
    that term-class errors cannot cancel.
12. As the tester, I want FOIS cardinality checked against the
    brute-force excitation count at pinned sizes, so that
    coupling-only enumeration is proven complete.
13. As the tester, I want cumulant reconstruction identities (4-RDM
    cumulants reproduce the lower RDMs), so that the NEVPT2 class
    inputs are certified before the solves.
14. As the researcher, I want honest meters per row (denominator
    floor per class, class dimensions, FOIS completeness), so that the
    intruder surface is measured, never assumed.
15. As the maintainer, I want CASPT2 (rung 4) specified but unbuilt,
    so that the level-shift dial class stays out until certified rungs
    demand it.
16. As the researcher, I want the results stored as a per-atom
    payload (not a CellKey axis), so that the campaign grammar is not
    widened by a post-processing analyzation.
17. As the researcher, I want the saturation view (the CAS spread vs
    the PT2-corrected spread across the window ladder), so that the
    static/dynamic mixing question has a quantified answer surface.
18. As the researcher, I want PT applied only to certified-symmetric
    references, so that the diagnostic reads the static part cleanly
    and never the broken-symmetry controls.

## Implementation Decisions

- **Module**: `src/atomterm/perturb.py` — analysis row
  (`analysis · stages · scans · perturb`). Docstring owns the
  contract; import edges downward only (system, detbasis, secondquant,
  operators); `natural.py`'s lazy-driver pattern is the model.
- **Interface** (grows rung by rung, one function + one result type
  per rung): `epstein_nesbet(sub, ci, roots=…)`, `dyall_e2(…)`,
  `nevpt2_strong(…)` — with `caspt2(…)` specified, unbuilt. Shared
  `PT2Result`: correction, per-class breakdown, meters.
- **Substrate**: the signature `(mol, mf, c_act, ncore, na, nb)` +
  the stored `ci_stack` (PySCF order; the bridge lives in operators).
- **Solution selection**: explicit `roots` (indices into the stored
  stack), default the ground root; term labels via the analysis
  readout. Declared, never inferred. Multi-state (quasi-degenerate)
  PT2 is out of scope — per-state corrections only.
- **Rung ladder**: 0 FOIS (coupling-only enumeration, Slater–Condon
  native) → 1 Epstein–Nesbet (H₀ = diag H, external) → 2 Dyall H₀
  (generalized Fock, inactive/active/external range separation) → 3
  NEVPT2-SC (7 classes, RDM stack to 4-RDM, cumulants + partial
  integrals, small per-class solves) → 4 CASPT2 (internally
  contracted, level shift — SPECIFIED, NOT BUILT: its level shift
  masks exactly the intruder information the diagnostic reads;
  NEVPT2's per-class floors are signal, CASPT2's shift would hide
  them. Joins later only on demand, with shift-sensitivity disclosed
  per row).
- **Validation first**: `scripts/validate_perturb.py` — the
  three-layer chain printed as certificates on C and Si at realistic
  sizes; runs BEFORE the table exists. This script is the run
  instrument (like `run.py`), not a pytest tier.
- **Artifacts**: `elements/<Atom>/perturbation.json` (payload store,
  keyed window × sector × core × root × variant) +
  `elements/<Atom>/PERTURBATION_TABLE.md` (the render). Separate from
  TABLE.md/cells.json — no CellKey axis, no regen grammar change.
- **PERTURBATION_TABLE.md schema — ONE table** (the correction
  rows; no separate summary table — the per-window read is derivable
  from the ground-root rows and lives in the README key instead):
  one row per declared (window, sector, core, root, variant):
  `window | sector | core | root (term) | ndet | E_ref | ⟨L²⟩ |
  variant | E_PT2 | E_total | ΔE_total | floor | w_ref | notes` —
  `root (term)` is the declared index + the analysis term label;
  `E_ref` on the analysis ruler (disclosed); `variant` ∈ EN / Dyall /
  NEVPT2 (every variant reads beside the others — the twin grammar);
  `floor` = min |E₀ − H₀μμ| (per-class min for NEVPT2) — the
  intruder surface; `w_ref` = 1/(1 + Σ|C⁽¹⁾μ|²) — the reference
  weight, the expansion-quality gate (below ~0.8 conventionally
  suspect); `notes` carries FOIS completeness, class dims, refusals.
  RAW columns only. Three meters total: R over rows, floor + w_ref
  per row — one decomposition read, two validity gates. (‖C⁽¹⁾‖∞
  stays an internal debug meter; no table contract.)
- **The README reading key** (README · Results, one home, brief):
  R = (ΔE_CAS − ΔE_total)/ΔE_CAS over the window ladder — R ≈ 1
  dynamic-dominated (the poverty was dynamic; the enlargement was
  absorbing it), R ≈ 0 static (the missing radial is a true
  near-degeneracy), intermediate = mixed (each window's E_PT2 share
  is the dynamic it already holds — read off the E_PT2 column's
  decline directly). Guards: R > 1 or R < 0 is the
  intruder/overshoot signature, flagged, never celebrated — read
  with the floor. Every verdict CONDITIONAL on the row's gates
  (floor, w_ref) — the ladder-internal honesty: no external ruler
  exists in-env, so "R = 0.9 at w = 0.95, floors clean" is a
  finding, "R = 0.9" alone is not. (NEVPT2 class-resolved shares
  stay payload-level: surfaced only if the ladder read is
  ambiguous.) The mechanisms' home; tables stay raw.
- **Render policy (the declared scope)**: the module accepts ANY
  stored solution; the table renders the DECLARED scope, whose
  default = the SA-CASSCF roots that deliver identifiable terms
  (clean ⟨L²⟩/⟨S²⟩ → a term symbol) across the window ladder. The
  m-incomplete controls (MO/MOC, broken ⟨L²⟩) are excluded from the
  default — a PT2 correction on a broken-symmetry reference answers
  nothing about the static/dynamic separation. Other classes (CASCI
  frozen-orbital rows, excited terms) join by declaration in the
  payload's scope list. Mirrors the NO-window refusal discipline.
- **Renderer**: `scripts/build_perturbation_tables.py` reading the
  payload store; byte-determinism + declared-roots pins.
- **No chkfiles, no caching**: everything derives deterministically
  from trusted substrates + stored CI.
- **External anchor**: none in env (PySCF 2.13 ships no NEVPT2; no
  forge package). Wave 1 stands on twins + analytic identities +
  literature class structure; installing pyscf-forge stays an open,
  user-gated option (recorded, not decided).

## Testing Decisions

- Unit tier: `tests/unit/test_perturb.py` on the frozen carbon
  fixture (`carbon` fixture — NO SCF) + synthetic substrates. Good
  test = one invariant, named for it; degenerate anchors first.
- The three-layer chain: (a) elementwise twins (EN sum vs solve;
  FOIS column vs brute-force `slater_condon_matrix` slice; cumulant
  identities) on inputs that exercise every term class; (b) external
  anchor deferred (env fact above; open option); (c) exact identities
  (whole-space → 0; one electron; empty active → MP2 limit under
  Dyall H₀ — the textbook known answer).
- Renderer pins: byte-determinism, declared-roots policy (the
  lane/lattice emitter pattern).
- The validation script carries the campaign-size correctness
  evidence (printed certificates), priced at doors — not a test tier.

## Out of Scope

- CASPT2 rung 4 construction (specified, unbuilt).
- Multi-state / quasi-degenerate PT2 (per-state corrections only).
- Campaign integration: CellKey axis, regen grammar, TABLE.md rows.
- Atoms beyond C and Si (the twins, Sc/Ti, U14e ride later waves).
- Spin–orbit, transition moments, g-factors.
- Installing pyscf-forge (open user-gated option, not this wave).

## Further Notes

- **The diagnostic logic** (the wave's physics goal, stated as
  falsifiable readings on the saturation view):
  1. *Dynamic dominates the spread* — E_PT2(compact window) is large
     and the corrected ladder FLATTENS (ΔE_total ≪ ΔE_CAS): the
     poverty was dynamic correlation in fact; the CAS enlargement was
     absorbing it variationally, never saturating.
  2. *Genuinely static* — E_PT2 is small or window-independent while
     the corrected ladder keeps the spread: the missing radial is a
     true near-degeneracy; only in-window variational treatment gets
     it (the @r1/@s1 verdicts stand as static facts).
  3. *The mixed regime* — E_PT2 shrinks monotonically as the window
     grows (the PT2 budget eaten by the CAS) with partial flattening:
     each window's PT2 share quantifies how much dynamic it already
     holds — the mixing, measured per window.
  On C/Si the poverty class (+250.8/+221.5 mEh) is the concrete
  instance; the U14e +923 mEh poverty and the 1 µH floor coincidence
  are the later-atom confrontations (out of wave-1 scope, the same
  instrument).
- **The physics payoff on the live record**: the window-ladder
  diagnostic on the stored C/Si CAS; a dynamic axis for the p² terms
  (³P vs ¹D differential correlation — the selection interface has
  real content).
- Rung-by-rung shipping keeps every commit certified: each rung lands
  with its twins/anchors in the same commit (lockstep rule).
- Intruder-freeness of NEVPT2 is a theorem our meters VERIFY
  (per-class denominator floors), never an assumption.
