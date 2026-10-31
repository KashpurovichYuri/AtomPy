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
- Part 2: the evaluator lane + pins 6–7.
- Part 3: declaration (cases term-targets, cells axis) + part 8.

## References

- The live probes and the mini-anchor of 2026-10-31 (session record;
  the demo script is the anchor test's skeleton).
- Block's with-pyscf docs (the `state_average_mix_` pattern); PySCF
  `fci.direct_spin1_symm`, `symm.param` (the XOR id convention).
- Ticket 31's gated U batches re-scope in this lane's light once the
  evaluator lands.
