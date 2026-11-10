# AGENTS.md — AtomTerm

Working agreement for coding agents on this repository.

## What this project is

**AtomTerm** classifies electronic states of atoms by their atomic term
symbols (`^(2S+1)L`) using *exact* angular momentum operators. Pipeline:
ROHF → CASCI → CASSCF → SA-CASSCF (PySCF), then a first-quantized analysis:
`Lz`, `L^2`, `S^2`, `N` built as plain matrices in the Slater-determinant
basis (second quantization appears only as the elementwise
cross-check twin), symmetry quantified by
`||[H, Lz]||_F` and `||[H, L^2]||_F`.

## Environment

```bash
source ~/miniconda3/etc/profile.d/conda.sh && conda activate pyscf
# or call ~/miniconda3/envs/pyscf/bin/python directly
pytest                      # from repo root; pythonpath=src is preconfigured
python scripts/run.py [atoms]  # term tables per declared atom
```

Python 3.12, PySCF 2.13, numpy 2.5, scipy 1.18 (note: `scipy.special.sph_harm`
is gone; use `sph_harm_y(n, m, polar, azim)`).

**Filesystem rule (user-set)**: using ANY directory outside this project
directory is restricted without the user's explicit permission — no scratch
files, tool caches, or artifacts elsewhere. Temporary files live inside the
repo (e.g. a `*.tmp` beside its target, removed once used). Plans live
in-repo as `wayfinder/tickets/` (plan-mode DRAFTS are confined to
`~/.opencode/plan` by the harness and move into the repo as soon as mode
restrictions lift). Mechanical enforcement: the global
`external_directory: ask` permission rule (user-side).

## Layout (deep modules, small interfaces)

| Module | Interface | Role |
|---|---|---|
| `system.py` | `AtomConfig`, `build_mol`, `sector_na_nb`, `ncore_of` | validates: single atom, spherical basis (origin by construction); `ecp` pairs basis+pseudopotential; the sector/core helpers raise on parity violations and take `mol.nelectron` (ECP-honest — never recompute from Z) |
| `angular.py` | `l_matrices`, `l_ladder_matrices`, `l2_matrix`, `l_matrices_analytic`, `real_to_complex` | one-electron `Lx/Ly/Lz/L²` AO matrices — TWO independent constructions certifying each other (libcint `int1e_irp` primary, exact tesseral algebra the cross-check, per-slot through `windows.shell_slots` so fused shells work); `real_to_complex` = the per-shell real→\|l,m⟩ transform |
| `irreps.py` | `IRREP_NAMES`, `IRREP_IDS`, `branch_d2h`, `det_irreps`, `orb_irreps` | the SO(3)→D2h dictionary — which irreps a term's 2L+1 components occupy (character algebra; a multiplicity can exceed one for L≥2 and is the state-average weight / per-solver root count) and which irrep every sector determinant carries (XOR of PySCF's D2h ids over occupied orbitals — a doubly-occupied pair cancels; aligned with detbasis order); `orb_irreps` = the ONE PySCF symm boundary; D2h is the working dialect (the SO(3) CI mode is Abelianized m-sums) |
| `detbasis.py` | `determinant_basis`, `sector_dimension`, `one_body_matrix`, `two_body_matrix`, `slater_condon_matrix`, `operator_apply`, `commutator`, `frobenius_norm` | **pure numpy, no pyscf** — Slater–Condon rules over explicit determinants; `slater_condon_matrix` = the fused one-pass production path, the name-twinned parts cross-certify with secondquant; `operator_apply` = the matrix-free twin (the same signed couplings accumulated into the output vector; one-body-only closures ride a singles-only enumeration — no N² container, no pairwise classification); ordering ≠ PySCF FCI strings (the bridge lives in operators; detbasis stays numpy-only) |
| `secondquant.py` | `annihilate`, `create`, `apply_one_body`, `apply_two_body`, `one_body_matrix`, `two_body_matrix`, `secondquant_matrix` | **pure numpy, no pyscf** — the same matrices via explicit a+/a application on occupation vectors; the elementwise twin of detbasis (random non-particle-symmetric-G pins certify the string contraction terms physical integrals are blind to) |
| `operators.py` | `build_operators`, `active_space_integrals`, `full_space_integrals`, `window_closure`, `core_anisotropy`, `ground_block`, `pyscf_ci_permutation`/`ci_in_detbasis_order`, `pyscf_apply` | the first-quantized assembly: H, Lz, L², S², Sz, N, the EXACT parity (exterior power), population labels; owns BOTH PySCF↔ours convention boundaries (the integral transpose + the CI-ordering bridge) and the measured-block family (closure certificate, core anisotropy, `ground_block` = the SA-anchoring/NO-averaging rule, SCALE-AWARE: the term identified by ⟨L²⟩ within 0.1 with the energy window spanning one percent of the different-term gap — the block is the WHOLE split multiplet under a residual anisotropy); `pyscf_apply` = the Hamiltonian through PySCF's C contractions at any sector size (the absorb-`.5` pairing; the one-body C leg does not exist — the kernel packs the triangle and would corrupt Lz); core folded into h1 (m-scalar density; `c_core` for optimized stages); `method=` switches the evaluator |
| `reference.py` | `ReferenceSpec`, `run_rohf`, `run_uhf`, `uhf_basin_census`, `uhf_contamination`, `run_frac_occ`, `run_folded` | the SCF reference engine: character-forced occupations + basin gate thread every runner (the cache/gate skeleton and the acceleration dispatch are shared; the UHF gate anchors on `occ_e_ref_uhf`, the census-measured UHF basin — `uhf_basin_census` measures it, any atom one call); the UHF class is a REFERENCE-LAYER product, not a cells axis (the folds consume one closed-shell core density — the broken-core fold is the declared future step); the folded reference (`run_folded`: the confined core FROZEN, the valence ROHF re-relaxed in its field, complement-local, MOM-gated, seeded from the substrate's valence — the lane's reference class; the pick rule DECLARED per atom-gauge (`picks=`), the draw judged by the converged-filling certificate — the uncontracted twins' folded landscape is basin-boundary draw-chaotic, measured 2026-11-09); `uhf_contamination` = the spin-honesty meter (the overlap formula, certified against the exact determinant-space S²); module Conventions own the caching-trust (cached chks TRUSTED, never re-gated) and the S-projector occupation-metric contract (`span_projector` is the one home) |
| `casscf.py` | `run_cassci`, `constrained_casscf(mo_full=, mask=, weights=, mc_opts=)`, masks, `run_irrep_casci`/`run_irrep_sa`, `CASResult` | the ONE deep optimizer driver — a mask is any `uniq_var_indices`-shaped callable; every recipe is a `constrained_casscf` call (the cases stages included — the one-line stage adapters retired at the ECP wave's close), never a new kernel path; SA = the two-pass pattern (ensemble kernel → exact roots readout in the frozen window); the symmetry-lane evaluators (`run_irrep_*`) are pinned-route twins on the substrate signature -- ONE pinned solver per `irreps.branch_d2h` component, no mask/checkpoint vocabulary; module Conventions own the mcscf mutation semantics and the full-MO contract |
| `windows.py` | `Window(.columns)`, `select_*` (incl. `select_span`, `select_primitives`), `shell_realization`/`realized_patterns`, `member_chains`, `distilled_ncas`/`window_distilled`, `primitive_chain`, `fold_shells`/`closed_electrons`/`NOBLE_CORES` (the fold vocabulary), `shell_weights`, `window_avas`, `span_overlap`, `span_projector`, `s_metric_roots`, `shell_slots`, `lowdin(tol=)`, `complete_orthonormal`, `build_full_mo`, `lm_sectors`, `l_per_row`, `radial_basis`, `mask_radial_rotations` | window construction + selection + the kind vocabulary + the SHELL LAYER (a chemical shell's AO realization is DECLARED band-by-band — basis-specific, certified by measurement; the identity default = the contracted behavior) and its WINDOW DISTILLATION (`window_distilled`: each shell as its declared k radial directions distilled from the reference's occupied MOs by the member-summed radial Gram, replicated on the m-set — m-complete at every k, rotation-closed whatever the capture; the SELECTION is a significance floor, not an ownership test — a hybrid-open valence column carries its shell's radial; a multi-band realization's chemically-sized window, the CAS label reads the delivered records' shapes; the FOLD stays full-family) (module docstring owns: purpose, the m-completeness contract, the five kinds, the internal layering — `shell_slots` the one shell walker, the S-metric algebra beneath, the two substrate builders; `member_chains` = the m-replication substrate's one home — the chains + the shared radial overlap block every distillation walks); `build_full_mo` takes AO rows OR column vectors through ONE project-renormalize-assemble treatment; `lowdin` raises on near-linear dependence |
| `core.py` | `purified_core_orbitals`, `confined_core_orbitals`, `closed_shell_core_energy` | the exactly rotation-closed core span (the pur substrate fixing the input-window baseline) + the closed core ON a fold span (the span's own closed-shell Roothaan filling the lowest `n_pairs` — the span itself on a 1:1 realization; `ncore` follows the PAIRS, never the span width) + the closed-shell core energy (the pur-fold bookkeeping constant; over confined columns it IS the projected-density form) |
| `basis.py` | `LadderSpec`, `run_ladder`, `valence_block`, `metric_radials`, `series_window` | the BASIS-PREPARATION layer (stage 0b, the contraction recipe): the configuration LADDER as basin-gated forced-occupation references (the head = the case's own trusted reference, never re-run — the gate MANDATORY, the wandering measured) + the prior axis's radial signal (the occupation-weighted member-summed valence density block, averaged per-configuration over the ladder — the density forms are rejected, the Gram/block average is the one that works) + `metric_radials` = the CAPTURE-MAXIMIZING metric eigenproblem (the `(block, s_r⁻¹)` transport, `eigh(LᵀBL)` with vectors through `L⁻ᵀ`; the `(block, s_r)` conjugation returns s_r⁻¹-distorted radials — the planted anchor pins it) + `series_window` = the tagged windows' construction (k from the tag, the kept-signal tail certificate, the discarded-fraction disclosure, per-config captures); the prior axis (@r/@s source tags, cells' `source_tag`) is the density-block family beside the incumbent `window_distilled` reference-source span-signal family — two declared families, never mixed in one comparison row |
| `stages.py` | `StageFacts`, `stage_facts` | ONE analysis pass over a finished stage: fold + analyze + energies vs PySCF + window purity + ortho divergence guard; module Conventions own the two ruler readings (the reconstruction is the ANALYSIS ruler; max\|dE\| decomposes into anisotropy + window definition + bookkeeping) |
| `scans.py` | `ScanSpec`, `scan_window`, `scan_state`, `radial_scan`, `member_tilt_map` | frozen-orbital scans, optimizer-free and DETERMINISTIC given the reference: m-symmetric radial curves (exactly closed windows) + member-tilt maps; `ScanSpec` = the declared per-atom anatomy carried by `CaseSpec.scan` |
| `analysis.py` | `analyze`, `configuration_summary`, `shell_occupations`, `format_report` | term symbols from ⟨L²⟩/⟨S²⟩ (half-integer S grid), commutator norms, configuration labels per root (PySCF-order CI in — the bridge runs inside; l-grouping via `lm_sectors`, no label parsing) |
| `cases.py` | `CaseSpec`, `Scenario`, `CAMPAIGN`, `scenarios_for`, `run_case`, `run_atom`, `run_carbon` | the DECLARED campaign (the intent manifest): `CaseSpec` fixes atom/reference/window/sectors/term-target/lane-legs/scenarios — and the render chemistry (window/fold display names, the filling note, the primary-family suppression) — and `shell_bands` (the SHELL LAYER: each chemical shell's declared AO realization, measured-certified) + `shell_radials` (the WINDOW-SIDE distillation depth, one k per shell, coverage validated at the declaration) before anything runs; `CAMPAIGN` is the ONE declaration — the cell grid, regen, and the renderers all derive from it; new atom = a registry row |
| `cells.py` | `CellKey`, `run_cell`, `load_cell`, `analyze_cell`, `cell_payload`, `spine_keys`, `window_grid`/`grid_keys`/`diagnostics`/`planned`/`lattice_keys`, `lane_keys`, `run_lane_cell`/`load_lane_cell`, `void_reason` | THE addressing layer: one frozen `CellKey` per grid point, the address grammar as chk name/payload key (a lane window carries its `@fold` tag — part of the identity; the reference field rides as an OPTIONAL ninth slot, appended only when `folded` — an eight-slot address is a neutral cell, the pre-axis store names untouched); VOID cells validated at the key; the declared key-set derives from `CaseSpec` (`window_grid` = the one render order); the scenario RUN GATE (`void_reason` — first enabled match, `--full` materializes); the pinned lane = the CI-level ladder legs as cells under the uniform construction (folds resolve through the shell layer — the chemical partition; the frozen core is the reference's OWN core — the occupied columns S-projected onto the fold span, joining the neutral zero by construction: the retired confined-Roothaan core measured +124.3/+38.7 mEh off on the twins, the repair left +6.1/+5.1 = the projected-density bookkeeping form; the lane payload records `core_pairs` + the in-span certificate); the self-owned store (`ci_stack` stays PySCF order — the bridge lives in operators); SA weights default to the measured ground block |
| `natural.py` | `ci_rdm1`, `natural_orbitals`, `occupation_groups`, `select_no_window`, `no_group_windows` | the NO mathematics: spin-summed 1-RDM from CI roots, NOs = density eigenvectors, occupation-GROUP selection never splits an m-set (scale-aware tolerance); symmetry-safe only from a term-AVERAGED root set — the source rides the PURIFIED span (substrate-independent, like the distilled columns) and every group prefix's closure is VERIFIED, the m-unequal density class REFUSED at the meter; the Window object + selectors live in `windows.py` |

## Hard-won facts (do not re-derive)

General claims + owners. Measured values, atom specifics, and campaign
history live in the test pins, `elements/`, and git history — a fact
carries a number only when the number IS the fact.

1. **AO labels are truthful but carry member suffixes** (`2px` is really
   px): prefix-match, never regex; fused shells (`nctr > 1`) unpack
   through the shell walker — `windows` (`select_shell_ao`, `shell_slots`).
2. **`int1e_irp` decoding**: `T[3a+b] = <mu| r_a d_b |nu>`, `L = -i` ×
   the antisymmetric part — `angular` (module docstring). Pinned by tests.
3. **Two integral conventions collide**: PySCF `int2e`/`ao2mo` are
   chemist, the project tensor physicist `<pq|rs>` = chemist `(pr|qs)`,
   transposed exactly ONCE at the boundary — `detbasis` (the tensor),
   `operators` (the boundary), the FCI-comparison test (the arbiter).
4. **Matrix products are exact only in an orthonormal basis** — the
   `S^{±1/2}` transport (`windows.s_metric_roots`); AOs of one contracted
   shell are self-orthonormal (`angular.l2_matrix` Notes; the
   m-completeness contract in `windows`).
5. **The core density entering any fold must be m-scalar** (the physical
   l=0 component — not a hack): the SCF core carries numerical
   anisotropy that otherwise splits multiplets at the closure floor —
   `operators` (the fold).
6. **mcscf objects mutate** -> `casscf` (module Conventions: the
   mutation semantics + the pass-mo-explicitly/fresh-object discipline;
   the SA roots pass is the worked instance).
7. **CASCI/CASSCF take the FULL MO matrix** and slice
   `[ncore:ncore+ncas]` -> `casscf` (module Conventions; a truncated
   matrix fails `nvir >= 0` in `check_sanity`).
8. **Known physics, not bugs**: single-root CASSCF breaks rotational
   symmetry by m-incompleteness of the optimized window — it stays
   l-pure (radial member mixtures lacking the full m-set, NOT l-leakage;
   l-purity is necessary, not sufficient); state-averaging over the
   degenerate M_L members restores it. Only the commutators measure
   closure itself.
9. **Second-quantized two-body convention**: `O2 = ½ Σ G[p,q,r,s] a+_p
   a+_q a_s a_r` (G physicist; annihilators reverse order) — elementwise
   identical to the Slater–Condon matrices (the twin tests). MO
   orthonormality lives in the AO metric (`c†Sc = 1`), never the
   Euclidean norm; folded H + explicit core energy == the fixed-core
   block of the full-space H.
10. **From-scratch runs are not bit-reproducible** -> `reference`
    (module Conventions: cached chkfiles TRUSTED, never re-gated; never
    judge a refactor by diffing two from-scratch runs — optimized stages
    amplify the draw chaotically).
11. **Generic pipeline** (`cases.run_atom`): Ms sectors derive from
    `cfg.spin = 2S` + active-electron parity — primary
    `((n+2S)/2, (n−2S)/2)`, plus `(n/2, n/2)` when n is even and S > 0;
    sector chkfile tags label M_S = (na−nb)/2 (na−nb itself would
    collide). New atoms = a new row in the atom-table test, no new code.
12. **Fold and reference classes**: a raw (un-averaged) core field splits
    multiplets in the quadrupole pattern — the averaged fold is exact;
    the fractional-occupation (average-of-configuration) reference is a
    DIFFERENT reference, not the ε→0 limit of averaging; SA kernels
    return only the ensemble energy (per-root energies from the roots
    pass — `constrained_casscf`); optimized stages must fold with their
    OWN inactive columns (`build_operators(c_core=)`), else
    E-bookkeeping breaks at O(Eh); in a complete-single-irrep sector
    [H, L²] = 0 trivially — use [H, Lz] and the splitting as
    fingerprints.
13. **An apparent multi-solution may be a construction bug**: a
    NON-orthonormal MO matrix fed to an optimizer diverges below the
    variational bound silently — check C†SC − I of anything fed to an
    optimizer (every stage carries `ortho_err`; the divergence guard
    fires in `stages`).
14. **`build_full_mo` is block-safe**: core Löwdin'd within its span,
    actives projected + renormalized; the correction is same-l and
    m-diagonal, so l-purity and m-sets survive. Naive complement filters
    admit spurious directions on non-S-orthonormal input — silently
    wrong shapes. max|dE| vs PySCF decomposes into core-field
    anisotropy + window definition + bookkeeping (exactly 0) — `stages`
    Conventions.
15. **Open-shell campaigns bring new physics, never new code paths**: S
    can be half-integer (`analysis` snaps to the half-grid; integer
    rounding mislabels doublets as singlets at noise-level ⟨S²⟩); the
    fractional-occupation block is selected by CHARACTER, never energy
    order (energy-ordered filling runs away on diffuse bases; the metric
    is the S-projector — `reference`); free CASSCF has real
    multi-solution basins incl. window-identity swaps — tests branch on
    window weights, never energies; SA anchors the member gauge but
    cannot repair an m-incomplete window; raw [noble-gas] cores polarize
    — commutators detect it, ⟨L²⟩ alone does not on every atom;
    quasi-degenerate ground-term orderings flip across SCF chaos —
    never pin them.
16. **The radial constraint axis** (`radial_basis` +
    `mask_radial_rotations`): rotations allowed only between
    same-(l, member) slots — the mask is a PAIR predicate ("same-l"
    alone would mix members). The member gauge of a closure-flat
    direction is chaotic on single roots (even term labels wobble); SA
    anchors it reproducibly. An escape not radially reachable is a
    window-identity change, not correlation hunger; a complete START
    window is not a complete CONVERGED window. Fused shells are
    contraction-outer/member-inner (unpack via `bas_nctr`, not label
    lists); the Löwdin gauge within a degenerate member set is
    span-preserving, but column identity needs member-preserving
    Gram–Schmidt; keep the strong radial directions labeled, pad the
    near-linear needles unlabeled.
17. **The CASSCF-rung input-window commutator baseline is the start
    object's core-orthogonalization defect** — mechanism + honest
    re-orthogonalization: `windows.build_full_mo` Notes; the fix (the
    rotation-closed purified span): `core`.
18. **Frozen-scan verdicts** (`scans`; optimizer-free, DETERMINISTIC —
    never pin optimizer output where a scan exists): on p-shell valence
    the full CASSCF gain is a one-parameter m-symmetric radial rotation
    at exact symmetry (a single-member tilt is exactly flat; tilted
    pairs cost energy); open d shells INVERT this — the m-symmetric
    curve only rises, member-asymmetric tilts descend (the Jahn–Teller
    drive behind broken single-root solutions; SA averages it away).
    Sector commutators have BLIND SPOTS ([H,Lz] cannot see an m=0 or
    equal-±pair tilt; maximal-Ms single-irrep sectors see nothing) —
    the one-particle projector pair (fact 19) has none.
19. **The one-particle closure meter is the universal certificate**
    (`operators.window_closure`): ‖[P,Lz]‖ and ‖[P,L⊥]‖ of the window
    span projector (S-metric) — [P,L⊥] sees EVERY member asymmetry;
    both ≈ 0 certifies closure completely. Span overlap is the IDENTITY
    axis, independent of closure — it says nothing about it. SA
    deviations decompose into window (all) + core-fold (none) — the
    attribution rides every cells row.
20. **The closed core span** — E_l ⊗ members (per l: top-k span of the
    member-resolved radial Gram, every radial replicated on ALL members,
    one global Löwdin) makes input windows close at machine precision
    and the core density exactly m-scalar -> `core`.
21. **Degenerate-shell references are layered chaos**: energy-ordered
    (aufbau) filling collapses from every start while a lower shell is
    open (the runaway), and at fixed occupations the open shell is
    MULTI-BASIN (guess-dependent) — character-forced occupations + the
    basin gate land the intended reference (`reference`); the gate also
    fixes downstream substrates. Leave real headroom on
    degenerate-but-numerical test bands.
22. **ECP-atom campaign verdicts** (the f-shell record; the `run_u*`
    scripts retired at the wave's close — git history): the gated-Newton ROHF reference
    determinant can be a numerically pure L=0 mean field — "below ROHF"
    then means beating a pure state, and windows below the
    representability threshold all sit above it; raw-ROHF MO spans
    NEVER rotation-close (m-delocalized open orbitals — the closure
    certificate excludes symmetric roots from raw-MO CASCI without
    running any CI); m-complete AO lattices win the large pairs on
    representability + correlation.
23. **The γ ladder and NO windows** -> `natural` (the m-scalar criterion
    + the scale-aware grouping tolerance).
24. **PySCF FCI ordering ≠ detbasis ordering** -> `operators` (the
    bridge) + `cells` (the store contract).
25. **Atomic-ROHF CASCI selection rules**: CASCI over a window
    containing the ROHF determinant is variational with EXACT equality
    unless a live coupling exists — singles are Brillouin-dead
    (stationarity), doubles dead unless same-l (the same-center
    different-l exchange integral vanishes angularly; parity kills
    more), so the only live channel is the same-l double shell. The
    lattice E column is AS-RUN for this reason (the reconstruction
    stays the ANALYSIS ruler — `stages`). NOs-from-ROHF ≡ the MO window
    (the ROHF density is idempotent); fractional occupations need the
    term-averaged CASCI density.
26. **The shell layer and the multi-column fold** (the Sc-unc/Ti-unc
    wave): a chemical shell's AO realization is BASIS-SPECIFIC and
    DECLARED (`CaseSpec.shell_bands` through
    `windows.shell_realization`), never read off labels — contracted
    bases do not automatically hold the 1:1 shell-to-block match
    (cc-pVDZ Sc/Ti: three d bands, six s bands; the certified lanes
    were 1:1 because their folds named single-block inner shells). A
    fold's ELECTRON count is chemical (`closed_electrons`), never its
    span width: the PRODUCTION frozen core (the cells' lane + free
    folds) is the reference's OWN core — the occupied columns
    S-projected onto the fold span (`cells._fold_core`), joining the
    neutral zero by construction (2026-11-09: the confined-Roothaan
    core measured +124.3/+38.7 mEh off the anchors on the twins —
    the bare core, solved without the valence, over-relaxes on a
    wide span; basis-stable per the contracted control).  The
    confined Roothaan itself stays `core.confined_core_orbitals` —
    `run_folded`'s construction and the certified class (on the
    whole AO space it reproduces the (2·n_pairs)-ion RHF exactly,
    pinned). The span carries the shell
    grammar: family-covering down to the CUSP, excluding every
    skipped or active shell's bands (the Roothaan fills the lowest —
    a norm-capture minimal list misses the tight primitives:
    measured 99.99% S-capture with +13 Eh core-energy error; the
    S-projector capture is a NORM statement, and Mulliken band
    populations run wild (±2) on near-dependent primitives —
    `scripts/band_composition.py` prints both but only the projector
    is bounded). The convergence certificate is the density-Fock
    commutator (the total energy settles quadratically earlier).
    Windows/fold partitions are CHEMICAL; the shells' band lists may
    overlap where a partition is declared non-disjoint, and the
    window projection against the confined core owns the shared rows.
27. **The second-order (Newton) SCF path freezes occupations from its
    early cycles** — a converged solution can therefore VIOLATE a
    character-forced occupation hook (a 3d-forced configuration
    converging onto a different configuration's solution): energy
    agreement alone cannot certify the configuration — the
    occupation post-check (the hook re-applied to the CONVERGED MOs
    must reproduce the run's own occupations) is mandatory wherever
    the accelerator runs under a forcing hook. The census instrument
    (`reference.rohf_basin_census`) carries it as the `occ_ok`
    column; `basis.run_ladder` post-checks every gated ladder run.
    Corollaries: the plain-DIIS arms can diverge on exactly the
    configurations that need forcing (the accelerator is declared
    per LadderSpec), and a ladder's forced configurations are
    multi-basin at draw-chaos scale — the gate anchors a CENSUSED
    basin (the occupation-safe rows), never a fresh draw's energy.
28. **The two distillation families and the radial-prior metric**:
    the incumbent `window_distilled` (the reference's span-signal
    Gram) and the prior axis's `basis.series_window` (the
    occupation-weighted density blocks, averaged over the
    configuration ladder) are DIFFERENT signal definitions sharing
    the member chains (`windows.member_chains`) — never mixed in one
    comparison row; the density-block radial eigenproblem must
    maximize the CAPTURE (the generalized problem `(block, s_r⁻¹)`,
    `eigh(LᵀBL)` with vectors through `L⁻ᵀ`): the `(block, s_r)`
    conjugation mixes the Euclidean coefficient pairing into the
    numerator and returns s_r⁻¹-distorted radials (the
    planted-radial anchor is the pin; both families agree on a
    single-configuration source — the consistency certificate). At
    k equal to the chain width the two conventions coincide (the
    span is the whole space — the CASCI total is convention-blind).

## Conventions

- Per-atom tables live in `elements/<Atom>/TABLE.md` (committed,
  log-provenanced): the atom heading + the ONE ladder table between
  the `ladder` markers + the `TABLE.csv` mirror (the generator:
  `scripts/build_methods_tables.py`, which refuses atoms without
  TABLE.md); README carries one worked example; cross-atom
  comparisons live in the element files. New atom = new element file + a
  data row, not README growth.
- Atomic units, ℏ = 1. Spin-orbitals: alpha = `0..norb-1`, beta = `norb..2norb-1`.
- Determinants: sorted tuples of occupied spin-orbitals.
- The substrate signature `(mol, mf, c_act, ncore, na, nb)` — the
  atom, its converged reference, the active columns, the core size,
  and the sector's (n_alpha, n_beta) — recurs through the analysis
  seam; functions document only their own knobs beyond it.
- Every operator enters the determinant space as (one-body matrix,
  two-body integral tensor) — never creation/annihilation algebra.
  Two documented exceptions: the explicit ladder operators
  (`spin_ladder_matrix`) apply a†/a directly — an inter-sector map has
  no sector-symbol representation; its signs are pinned by the
  twin-construction tests — and the many-body parity (`parity_operator`)
  is multiplicative — the exterior power, pinned by the hand-derived
  mixed-window matrix and the involution on closed spans.
- New physics code ships with a test that would catch a sign or phase error
  (spectra, commutators, or independent-construction agreement).

## Tests

`pytest` runs 303 tests. The fast tier: `pytest -m fast` = 209
tests in seconds (machine-light selection). The heavy tail is
HOURS-class now (the `test_symm_tm` fixture carries the (8,13)
SA-leg optimizer pairs, 1–2 h each — the old "~15 min" record
predated their admission and was never re-priced; caught 2026-11-05).
**The per-commit gate** (user-directed 2026-11-04): the fast tier
plus the diff's targeted files, together 5–10 minutes MAXIMUM — the
suite's purpose is correctness confidence, not endurance.
**The full-suite run is RETIRED** (user-directed 2026-11-05): no
wave-close full run — the gates close a wave. A heavy pin runs
ISOLATION, on evidence (a targeted red, a touched file, an open
physics question), detached (`setsid`) with `OMP_NUM_THREADS`
capped (8): uncapped, 255 BLAS threads thrash the box. A rename
round touching stage tags owes the slow-tier FILES carrying the
renamed lookups — targeted, not the whole suite. A red alone is
never the verdict — isolation reruns are.
**The boundary the suite guards** (user-directed 2026-11-05):
tests verify CONSTRUCTIONS at pinned small sizes — seconds to
minutes; heavy compute belongs to the ATOMIC CALCULATIONS (the
campaign runs — the project's goal), priced per batch and gated at
doors. A test that needs hours is a run in a test costume: the one
known offender (the `test_symm_tm` SA legs) re-scopes to
store-anchored verification — fresh CI legs stay, the SA verdicts
load the committed route-batch stores and re-verify their meters.

Layout and tiers:

- `tests/unit/` = the seconds tier: pure numpy/libcint, synthetic
  substrates, and the FROZEN carbon ROHF (`tests/carbon_reference.py`,
  repr-round-trip-exact source, loaded by the root `carbon` fixture —
  NO SCF runs in the unit tier; downstream bit-stable, double-run
  sha256-certified; regenerate only on basis/default changes via
  `scripts/generate_fixture_carbon.py`) plus its D2h twin
  (`tests/carbon_reference_d2h.py`, the `carbon_d2h` fixture — the
  symmetry lane's substrate; the tier also runs small
  symmetry-pinned FCI solves, no SCF).
- `tests/integration/` = the campaign chemistry (TM references,
  CASSCF/SA, term tables); conftest split the same way.
- The fast/slow MARKS are the machine axis; the directories the human
  one. Every file marks at module level except `test_cells`, which mixes
  tiers per function as stated policy (1:1 with its module).
- Shared DATA and helpers live in plain modules with UNIQUE basenames
  (`integration/refs.py`, `unit/helpers.py`) — never import the
  `conftest` basename directly: pytest claims it per-directory, and
  direct imports resolve order-dependently across tiers.
- Test files map 1:1 to src modules; tests are named for the invariant
  they pin (one invariant per test). Cross-module SEAM files are the
  stated exception (`test_ecp_wiring`, the `test_symm_*` lane files —
  the anchor spans two modules by design). Each file's own docstring
  carries its content — there is no catalog here.

Suite-flake protocol (the environment findings):

- Session-scoped rungs fixtures mean a chaotic SCF draw can fail
  whichever tight band grazes it — rotating failures across runs, each
  isolation-green on redraw.
- FREE-MACHINE AMPLIFICATION: with the machine unloaded, full-suite
  sessions intermittently crash a purified rung inside PySCF's Davidson
  (NaN CI overlaps; one fixture death cascades) while isolated reruns
  are green — load-dependent chaos, NOT a regression.
- **A red suite alone is never the verdict — isolation reruns are.**
  Widen only a grazed band, and only with the observed fresh-draw value
  in a comment.

### The module walk

The file-by-file walk (ticket 30, completed 2026-10-30) reviewed every
src module + test file: purposes, contracts, reuse, naming, per the
Refactor/Verification rule lists below. For any future walk: read the
module and its test file together; restructure per Refactor; certify
per Verification; collapse facts whose home a module now owns; update
the layout row in the same commit.

## Refactor

The rules every restructuring of src/ follows. Certifying the result is
the Verification section's job.

1. **Contract over tolerance.** The goal-sentence comes before any
   math: state what we want and why, then the equations as the
   realization. Evaluate every line: conceptual or incidental --
   would a caller or tester miss it? Docstrings carry what a caller
   must know: purpose and goal first -- the module docstring owns the
   purpose, the shared contract, and the internal layering, each
   contract stated once and referenced elsewhere; per function, a
   one-line brief + Parameters / Returns / Raises / Notes / Examples.
   Defaults are contract: state what the default physically is.
   Tolerance envelopes and arrival history are comments at the point of
   use, or commit-message material. Definitional slots (module
   docstrings, function contracts, section comments) carry the RULE
   stated generally — no atom names, no measured campaign values, no
   fact/ticket tags; instances live in Examples/Notes sections, and the
   measured record lives in tests and committed artifacts (`elements/`,
   `results/`) — never in definitional prose (amended 2026-10-30,
   user-directed).
   The standard is concrete: measure against the exemplar passes
   (`system.py`, `windows.py`, `secondquant.py`, `detbasis.py`) --
   one sentence carries one fact, every fact sits in its slot,
   formulas that carry a convention get their own line, the
   load-bearing word is highlighted sparingly (ONE, NEVER, EXACTLY),
   and zero `Parameters` sections in a walked file is a failed pass.
2. **Names join their family.** Selectors take `select_*`, window
   constructors take `window_*`, pure linear algebra takes plain
   nouns. A name that misstates its function renames together with
   its consumers in one commit.
3. **One home per idea.** A formula or idiom written twice gets one
   owner. The deletion test admits the extraction: complexity that
   vanishes was a pass-through; complexity reappearing across several
   callers earns the keep. A seam with one consumer stays
   hypothetical.
4. **Real separation over wrappers.** A wanted split divides the
   code -- separate evaluators over one shared classification, not
   forwarding shells. Twin modules expose identically-named parts so
   the cross-checks are name-symmetric. A nested helper that is not
   itself the product moves to module level; a closure that IS the
   product stays.
5. **Illegal states become unrepresentable.** Delete the fake degree of
   freedom rather than guard it; guard what the interface can still
   express, failing fast at the seam that owns the invariant, the
   message stating why.
6. **Lockstep.** A surface or behavior change updates the docstrings,
   `__all__`, the layout row, consumers, and test imports in the same
   commit -- documentation that lies is fixed in the commit that
   broke it. Map the blast radius first (grep the consumers); prefer
   whole-file moves to function surgery.
7. **Sectioned grouping.** Functions sit under a marker naming the
   role; new functions join their section.

A restructuring pass ends by walking this list rule by rule; an
unapplied rule is either applied or explicitly declined in the commit
message.

## Verification

The rules every change to physics-bearing code is verified by. The
stance: test as much as possible -- a doubt is a test waiting to be
written, and the degenerate case is the strongest anchor.

1. **Test shape.** Tests are named for the invariant they pin -- one
   invariant per test; name + content state the purpose. Test
   docstrings are one line or absent; narrative is commit-message
   material. Explanations that are content live as code comments at
   the point of use; fixture docstrings may carry the construction.
2. **Interfaces first.** Tests exercise the owning module's public
   surface; private functions get micro-pins only when they are the
   module's product or the public seam cannot reach them (a forced
   parameter). Test helpers that repeat get one home in the tier
   conftest.
3. **Every code path reached.** A construction is verified only by
   inputs that reach each branch: when production defaults hide a
   path at test sizes, force the parameter and pin it.
4. **The three-layer chain** (numerical constructions): (a) elementwise
   identity of independent constructions, on inputs that exercise
   every term class -- random non-symmetric integrals, not only
   physical ones; (b) an external anchor for the conventions the
   constructions share, in more than one symmetry class (open-shell
   and closed-shell references); (c) exact analytic identities
   (spectra, commutators, sum rules). A shared-convention error
   cancels in (a); only (b) and (c) catch it.
5. **Degenerate anchors and known answers.** Seek degenerate cases
   actively -- as many as constructible. Pin answers known from
   reference knowledge (textbook spectra, Clebsch factors) and
   coincidences expected across paths (two representations, two
   constructions). Exact-by-construction cases (one electron; the
   empty core; the whole space as window) carry answers without any
   reference: any deviation is a bug, not physics.
6. **Certified, not assumed.** A behavior-preserving refactor proves
   it: bitwise A/B against pre-refactor arrays on pinned references
   (deterministic constructions only; same-process for the CASCI
   kernel — cross-process bitwise A/B is invalid there, gauge chaos);
   elementwise pins for algorithm restructures; AST comparison with
   docstrings neutralized for docstring-only passes
   (`scripts/check_docstring_only.py`); the full suite before
   physics-bearing commits. Mathematically idempotent steps still
   count -- they perturb degenerate gauges and break bit-exactness.
7. **Measured pins.** Bands cite the observed fresh value in a
   comment; docstring examples run once before being written; hand
   derivations defer to the machine certification.

A verification pass ends the same way: rule by rule, applied or
explicitly declined.

## Planning

`wayfinder/map.md` is the navigator: live **Tickets**, open **Issues**
(dispositions, watchlist), and **Examples** (pointers into `elements/`).
Update it when a decision is made; keep tickets one-session sized.
Resolved tickets are DELETED (git history is the archive); each live
ticket is self-contained (no conversation context needed).

State: the atom campaign (C, Si, Sc, Ti), the project walk (ticket
30), and the ECP wave (ticket 31: the U14e campaign, the shell layer
+ the multi-column fold validated on Sc-unc/Ti-unc and ported, the G5
retirement — `run_u*` and the stage adapters deleted; git history is
the archive) are COMPLETE.  Gated measurement batches (T2/M) run
only on the user's enable.  The queued CASCI review (map.md Issues)
precedes the next physics wave; beyond it the candidate list:
spin–orbit, transition moments, g-factors.

Review protocol after each ticket wave, then a design survey:

- **Two-axis review** — run the `code-review` skill over
  `git diff <fixed-point>...HEAD`, fixed point = the last reviewed
  commit. Spec source = the wave's wayfinder ticket; docs commits have
  no file spec — review them only against this working agreement.
  Standards source = this file; the skill's smell baseline applies as
  judgement calls only, and this file overrides it.
- **Design survey** — apply the `codebase-design` skill vocabulary
  (module/interface/depth/seam/adapter, the deletion test) over `src/`.
- **Findings** — fix what is worthwhile now and commit; record the rest
  (candidates, declinations, and why) in map.md under Issues with the
  date. Widen a grazed test band only with the observed fresh-draw
  value in a comment.
