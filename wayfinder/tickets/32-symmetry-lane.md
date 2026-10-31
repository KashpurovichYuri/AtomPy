# Ticket 32 — The symmetry lane: D2h-pinned CI as the second route

## Goal

Add the independent route to L²-symmetric states: PySCF's
symmetry-pinned CI (D2h `wfnsym` solvers, `state_average_mix_` over a
term's components) evaluated and certified by our exact meters.
Two routes, one comparison — they never merge:

- route A (ours): diagonalize the symmetry-blind Ms sector, MEASURE
  symmetry with the exact operators, restore by SA over roots;
- route B (theirs): CONSTRUCT symmetry by group theory — one solver
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
7. The d-shell stress case (Sc/Ti): route B's guaranteed component
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
    ~3e-11 over three fresh draws); route B's optimized windows
    close orders tighter (C: lz 1e-15 vs 2.3e-6; Si: 5.1e-8 vs
    1.5e-5) — the irrep-blocked rotations cannot leak the m-set.
  - PROBED PySCF fact: state_average_ on a symmetry-adapted molecule
    forces one wfnsym per solver (its wrapper passes wfnsym through,
    and None is GUESSED to one irrep) — the symmetry-blind SA exists
    only on C1, so the head-to-head CROSSES SUBSTRATES: route A
    native (C1, the production driver) vs route B native (D2h twin).
    Recorded in the test's comment; the lane evaluator docstrings
    carry the never-bit-compare-orbitals rule.
- Part 3: the evaluator lane is DONE (lean); what remains of the
  original part-3 scope: declaration (cases term-targets, cells
  axis) + step 8 (gated campaign tables + the TABLES contract
  redesign), plus ladder step 7 (Sc/Ti stress) whenever the TM
  references are enabled.

## Step 7 spec (the next session -- agreed 2026-10-31)

One file `tests/integration/test_symm_tm.py` (unique basename), Sc +
Ti sections, the Si pattern (native substrate pairs: C1 route A + a
D2h-twin route B, fresh ROHFs under pytest tmp; the D2h constraint
makes the degenerate-shell reference deterministic, and any converged
reference works -- the m-scalar fold + the m-complete window absorb
the broken sphere).  Windows/sectors derive from the CAMPAIGN rows
(the one-declaration rule).  The stress content, beyond C/Si:

- **The multiplicity-2 branching goes live**: Sc's ground term ²D
  (L=2, g) branches `{Ag: 2, B1g: 1, B2g: 1, B3g: 1}` -- the first
  per-solver `nroots=2` through the runner; Ti's ³F (L=3, g) goes
  deeper (7 components).  The anchor pin covers it via block sizes.
- **nroots-hoped vs guaranteed components**: route A's Ti SA hopes
  the 7 lowest (3,1)-roots are the ³F set (¹D can interleave); route
  B guarantees each.  The head-to-head is CONDITIONAL and honest:
  measure route A's roots with our meters FIRST; assert the ensemble
  identity iff the root sets coincide; a differing ensemble is
  itself the finding (document, never force).
- **Certificates on the hostile atom**: route B's mix window at/near
  the floor vs route A at tolerance -- the C/Si finding re-tested
  where member tilts break free runs (fact 18).
- Pins per atom: runner anchor (blocks == our H, meters exact),
  component degeneracy, the conditional head-to-head, closures.

Declined up front: route-B single-irrep runs (fact-18 territory
already covered for route A by test_treatments).  Suite work, NOT
gated; full suite at the boundary; the two-axis review after.

## References

- The live probes and the mini-anchor of 2026-10-31 (session record;
  the demo script is the anchor test's skeleton).
- Block's with-pyscf docs (the `state_average_mix_` pattern); PySCF
  `fci.direct_spin1_symm`, `symm.param` (the XOR id convention).
- Ticket 31's gated U batches re-scope in this lane's light once the
  evaluator lands.
