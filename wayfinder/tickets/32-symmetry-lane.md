# Ticket 32 — The symmetry lane: D2h-pinned CI as the second route

## Goal

Add the independent route to L²-symmetric states: PySCF's
symmetry-pinned CI (D2h `wfnsym` solvers, `state_average_mix_` over a
term's components) evaluated and certified by our exact meters.
Two routes, one comparison — they never merge:

- the free route (ours): diagonalize the symmetry-blind Ms sector, MEASURE
  symmetry with the exact operators, restore by SA over roots;
- the pinned route (theirs): CONSTRUCT symmetry by group theory — one solver
  per D2h component, weights from the SO(3)→D2h branching.

The shared surface is the existing seam: the substrate signature
`(mol, mf, c_act, ncore, na, nb)`, the full-space CI-vector bridge,
the analysis meters. No new kernel path.

## Why D2h (probed 2026-10-31)

- An atom with `symmetry=True` gets SO(3), but PySCF's CI symmetry
  machinery is Abelianized there (m-sum blocks) — `wfnsym` cannot pin
  L. D2h pins the |M_L| classes properly; it is the working dialect.
- PySCF's D2h irrep ids are the Z2³ XOR encoding (Ag=0 … B3u=7;
  irrep product = XOR) — a determinant's irrep is the XOR over its
  occupied spatial orbitals; doubly-occupied pairs cancel.
- `wfnsym` solvers return FULL-space CI vectors (zeros on forbidden
  determinants) — the detbasis bridge consumes them unchanged.
  Demonstrated live: C 2p² triplet sector, the three g components
  degenerate to 1e-9, our meters reading L(L+1)=2 and S(S+1)=2
  exactly on every one.
- A configuration's inversion parity gates which irreps are reachable
  at all (p² → g only; s¹p¹ → u) — reachability is a contract, not a
  guess (PySCF raises `WfnSymmetryError` on the rest).

## Design

- New module `irreps` — the SO(3)→D2h dictionary, three functions:
  - `branch_d2h(L, parity)`: which irreps a term's 2L+1 components
    occupy, with multiplicities (character algebra; anchored at real
    shells). For L≥2 an irrep can host TWO components (e.g. Ag holds
    both z² and x²−y²) — the multiplicities are the SA weights and
    per-solver `nroots`.
  - `det_irreps(orbsyms, na, nb)`: per-determinant irrep id in
    detbasis order (the block identification for our own matrices).
  - `orb_irreps(mol, mo_coeff)`: the PySCF symm boundary (labels →
    XOR ids).
- The lane evaluator (pyscf-native CASCI/CASSCF +
  `state_average_mix_`), the cells method axis, and the cases
  term-targets come in the NEXT parts — after the dictionary is
  pinned. `AtomConfig` gains its symmetry knob when the production
  path needs it, not before (the fixture builds its own mol).

## Verification ladder (unit first, integration next, campaign gated)

1. Branch multiplicities sum to 2L+1 (L=0..8, both parities);
   hand-known p/d/f patterns; agreement with `label_orb_symm` on real
   shells. [unit]
2. The flip-derived character table equals PySCF's; the irrep
   product structure equals id XOR (all pairs). [unit]
3. `det_irreps` support equals the nonzero support of the
   `wfnsym=Γ` CI vector, per reachable irrep. [unit]
4. Doubly-occupied spatial orbital = the identity irrep (the (1,1)
   cancellation). [unit]
5. THE ANCHOR, on a frozen D2h-C fixture (repr-round-trip source, no
   SCF in the unit tier): per-irrep PySCF CASCI energies ==
   eigenvalues of our H restricted to the irrep block (to solver
   tol); our meters on their states read exact L(L+1), S(S+1); the
   term's components degenerate to ~1e-9; the full sector spectrum
   == the union of blocks. [unit]
6. SA head-to-head: our two-pass SA vs `state_average_mix_` CASSCF —
   ensemble energies agree, both densities pass closure certificates
   at the same floor. NEVER bit-compare orbitals across the routes
   (degenerate gauges — facts 10/15). [integration, part 2]
7. The d-shell stress case (Sc/Ti): the pinned route's guaranteed component
   set vs our nroots-based SA. [integration, part 2]
8. Campaign tables via the lane + the TABLES.md contract redesign,
   both records in hand. [gated on the user's enable]

## Parts

- Part 1 DONE (2026-10-31, commits e22260e + 16c3d4d): the `irreps`
  dictionary (branch_d2h/det_irreps/orb_irreps, character algebra
  anchored at PySCF's table and at real shells), the frozen D2h-C
  twin (generator `d2h` mode, double-run sha256; the C1 fixture
  untouched), and the unit pins 1–5 — THE ANCHOR at machine
  precision: pinned roots == our restricted-H eigenvalues (4e-16),
  term purity exact (2e-15), components degenerate (0.0), the
  spectrum == the block union. Review round done (two-axis); full
  suite 198 GREEN.
- Part 2 DONE (2026-10-31, the lean scope): the character-algebra
  and XOR math land in the docstrings (Notes carry the four chi
  formulas + the projection; det_irreps Notes carry the Z2^3
  argument); the lane evaluators in casscf — `run_irrep_casci` (one
  pinned block solve per component, the shared folded h1) and
  `run_irrep_sa` (state_average_mix_ + the two-pass readout); pins:
  the runner == the anchor on frozen C, and THE HEAD-TO-HEAD (C in
  the unit tier, Si in integration).
  - The head-to-head FINDING: ensembles identical (C 2e-12, Si
    ~3e-11 over three fresh draws); the pinned route's optimized windows
    close orders tighter (C: lz 1e-15 vs 2.3e-6; Si: 5.1e-8 vs
    1.5e-5) — the irrep-blocked rotations cannot leak the m-set.
  - PROBED PySCF fact: state_average_ on a symmetry-adapted molecule
    forces one wfnsym per solver (its wrapper passes wfnsym through,
    and None is GUESSED to one irrep) — the symmetry-blind SA exists
    only on C1, so the head-to-head CROSSES SUBSTRATES: the free route
    native (C1, the production driver) vs the pinned route native (D2h twin).
    Recorded in the test's comment; the lane evaluator docstrings
    carry the never-bit-compare-orbitals rule.
- Part 3: the evaluator lane is DONE (lean); what remains of the
  original part-3 scope: declaration (cases term-targets, cells
  axis) + step 8 (gated campaign tables + the TABLES contract
  redesign), plus ladder step 7 (Sc/Ti stress) whenever the TM
  references are enabled.

## Step 7 spec (the next session -- agreed 2026-10-31)

DONE 2026-10-31 (commits faf1adb + ff62a3d, review round included).
`tests/integration/test_symm_tm.py`, Sc + Ti parametrized, the Si
pattern (native substrate pairs: C1 the free route + a D2h-twin the pinned route,
fresh ROHFs under pytest tmp; any converged reference works -- the
whole active set sits in the window, so the fold sees only the
closed noble-gas core).  Windows/sectors/weights from the CAMPAIGN
rows (sector_na_nb; spec.sa_weights); the ground-term table is test
data (the test_atoms kind of row).  The spec as written:

One file `tests/integration/test_symm_tm.py` (unique basename), Sc +
Ti sections, the Si pattern (native substrate pairs: C1 the free route + a
D2h-twin the pinned route, fresh ROHFs under pytest tmp; the D2h constraint
makes the degenerate-shell reference deterministic, and any converged
reference works -- the m-scalar fold + the m-complete window absorb
the broken sphere).  Windows/sectors derive from the CAMPAIGN rows
(the one-declaration rule).  The stress content, beyond C/Si:

- **The multiplicity-2 branching goes live**: Sc's ground term ²D
  (L=2, g) branches `{Ag: 2, B1g: 1, B2g: 1, B3g: 1}` -- the first
  per-solver `nroots=2` through the runner; Ti's ³F (L=3, g) goes
  deeper (7 components).  The anchor pin covers it via block sizes.
- **nroots-hoped vs guaranteed components**: the free route's Ti SA hopes
  the 7 lowest (3,1)-roots are the ³F set (¹D can interleave); route
  B guarantees each.  The head-to-head is CONDITIONAL and honest:
  measure the free route's roots with our meters FIRST; assert the ensemble
  identity iff the root sets coincide; a differing ensemble is
  itself the finding (document, never force).
- **Certificates on the hostile atom**: the pinned route's mix window at/near
  the floor vs the free route at tolerance -- the C/Si finding re-tested
  where member tilts break free runs (fact 18).
- Pins per atom: runner anchor (blocks == our H, meters exact),
  component degeneracy, the conditional head-to-head, closures.

Declined up front: pinned-route single-irrep runs (fact-18 territory
already covered for the free route by test_treatments).  Suite work, NOT
gated; full suite at the boundary; the two-axis review after.

## Step 7 record (what was found)

- THE ANCHOR holds at machine precision on both atoms, and the
  multiplicity-2 branching went live: Sc's ²D Ag block hosts TWO
  components, Ti's ³F each B_g two (per-solver nroots=2, the first
  live exercise).  Pinned lowest-per-block roots == our restricted-H
  eigenvalues; components degenerate to ~1e-14 on the m-scalar fold.
- THE HEAD-TO-HEAD SPLITS BY ATOM, and the finding refines the
  spec's gate: MEMBERSHIP ALONE IS NOT ENOUGH.
  - Sc (d¹): the routes agree as on C/Si -- ensembles to 1.1e-8,
    multiplet formed (ptp 3.1e-7).  Bonus reading: the symmetric
    free-route window's certificate runs 5.2e-3 while the ensembles
    agree to 1e-8 -- the one-particle certificate sees f-direction
    drifts the energy cannot feel (hence _CERT_A = 2e-2).
  - Ti (d²): all 7 free-route roots still READ ⟨L²⟩=12 while the
    free optimizer had performed a CORE<->ACTIVE SWAP (see the
    mechanism addendum below): certificates O(1), multiplet split
    3.1e-3, span identity vs the raw shell 1e-4, ensemble 5.2 mEh
    BELOW the pinned route's.  The production Ti table records the
    same free-route basin (the swap row: comms 1.1/6.0, Elec
    s0 p2.0 d2.0).  The pinned route's D2h-blocked rotations FORBID
    the swap (it is cross-irrep: 4s is Ag, 3p is u): declared
    identity kept (s1.0 d5.0), components degenerate to 3.8e-7,
    certificates at 7e-6..1e-5 -- the guaranteed component set at
    ~5 mEh.  The symmetry-vs-variation trade, measured.
- THE HONEST GATE is the multiplet DEGENERACY (member AND formed);
  the divergence branch asserts ens_a < ens_b -- rigorous in BOTH
  shapes (free variation in the escape; displacement-below in the
  interleave: the intruder enters below the term root it pushes
  out) -- and never forces the finding's shape.  The free route's
  certificate in the divergence branch is measured, not pinned
  (basin chaos).
- Suite gates: TM file green on repeated draws; full suite 206/206
  GREEN (three documented rotating classes across the first three
  runs, each isolation-green).  Counts 206/132 in AGENTS; the
  seam-file exception now names the test_symm_* family.

## Mechanism addendum (2026-10-31, the user's 'why Sc but not Ti')

THE ESCAPE IS A CORE<->ACTIVE SEMI-CORE SWAP, measured on the
inactive block: the Ar start reads l-content s:3 p:6; the free
SA's converged inactive reads s:4 p:5.  The 4s pair rotated INTO
the fold; one semi-core 3p orbital -- with its electron pair,
0.998 pure 3p by S-metric AO content -- rotated INTO the active.
The gain is the 3p-3d semi-core correlation (the TM M-shell
effect); the price is a torn m-set (only one 3p m-function fits in
ncas=6).  The Sc/Ti asymmetry is the semi-core appetite scaling
with the open-d count (d1: the swap does not pay; d2: 5.2 mEh).
The pin (test_symm_tm head-to-head) measures the swap detector
and asserts the mechanism in every branch arm.

Lessons from the hunt (both caught by the meters/scales):
- The first hypothesis ('the optimizer wants the p shell; the
  m-complete version should pay') DIED: raw {3d,4p}/{3d,3p}
  probes came out ~300 mEh above with scrambled terms -- AND were
  malformed (fold-9 dropped the 4s pair: 20-electron problems).
  e_cas is fold-relative; only totals compare across folds.
- 'A full 3p shell in' (my first record fix) was also wrong --
  one column, not a shell; and the pair in it is the SEMI-CORE
  pair moving in with its electrons, not the 4s pair relocating.
- THE STEP-8 INPUT (gated): the m-complete symmetric version of
  the swap is the semi-core window CAS(6,8) {3d,3p} with fold
  {1s..3s,4s} on Ti (all six 3p electrons active) -- a candidate
  planned window alongside (4,10)/(4,14).  The radial-transplant
  (β) and NO-laundry (α) experiments stay queued; the cross-route
  window-intake decision is recorded in map.md Issues (D2h-born
  sources only).

## The CAS ladder (2026-10-31, the user's 'test different CAS' -- Ti, pinned, CI level, frozen shells)

All totals = component e_cas mean + closed_shell_core_energy (the
bookkeeping constant; fold-independent by construction -- fact 9's
fixed-core-block identity -- so cross-fold comparable).  One
substrate (the D2h ROHF), one fold rule; (4,·) windows PLACED
through build_full_mo (raw columns carry the documented ~85 mEh
Ti window-definition error -- the raw (4,14) run even read
<L^2>=13.66, a non-eigenvalue: the overlapping fold breaks [H,L^2]);
the (8,8) window is overlap-free by construction (its fold excludes
3p).

| window | E, Eh | vs declared | term |
|---|---|---|---|
| (4,6) {3d,4s} declared | -848.3912 | 0 | 3F exact |
| (4,10) {3d,4d} planned | -848.2045 | +187 mEh | D/H scramble |
| (4,14) {3d,4s,4p,4d} planned | -848.4349 | -43.7 mEh | 3F exact |
| **(8,8) {3d,3p} semi-core** | **-848.4520** | **-60.8 mEh** | **3F exact** |

THE VERDICT: the m-complete symmetric semi-core correlation (the
CAS(8,8) {3d,3p} window, fold {1s,2s,2p,3s,4s}, sector (5,3),
branch(3,g) -- an EXACT degenerate 3F, ptp 5e-14) is worth ~61 mEh
-- the free SA's torn escape (5.2 mEh) was reaching for a small
corner of a large pie, and the semi-core window beats even the
valence (4,14) window by 17 mEh.  Readings for the planned-cells
design: (4,10) as-declared (fold Ar(9), the 4s pair in the dead
zone) is NOT a 3F window -- its lowest g-states scramble D and H
components; the campaign's planned (4,10) cells have their own
constrained semantics and should not be read against this ladder.

NEXT for this thread (suite work + gated): CI-level pins of the
ladder (deterministic; the 17/61 mEh margins are robust), the SA
level on the winners ((8,8) first), and the campaign-declaration
question for the semi-core window (a new planned-window row).

REMAINING for the lane: step 8 + the declaration wave (cases
term-targets, cells axis, campaign tables via the lane, the Ti
TABLES compression with both records) -- GATED on the user's
enable; the ticket-31 interleave alongside.

## References

- The live probes and the mini-anchor of 2026-10-31 (session record;
  the demo script is the anchor test's skeleton).
- Block's with-pyscf docs (the `state_average_mix_` pattern); PySCF
  `fci.direct_spin1_symm`, `symm.param` (the XOR id convention).
- Ticket 31's gated U batches re-scope in this lane's light once the
  evaluator lands.
