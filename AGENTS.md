AGENTS.md — AtomTerm

Working agreement for coding agents on this repository.

## What this project is

**AtomTerm** classifies electronic states of atoms by their atomic term
symbols (`^(2S+1)L`) using *exact* angular momentum operators. Pipeline:
ROHF → CASCI → CASSCF → SA-CASSCF (PySCF), then a first-quantized analysis:
`Lz`, `L^2`, `S^2`, `N` built as plain matrices in the Slater-determinant
basis (second quantization appears only as the elementwise
cross-check twin), symmetry quantified by
`||[H, Lz]||_F` and `||[H, L^2]||_F`. Atoms delivered: C (2p²), Si (3p²),
Sc (3d¹4s²), Ti (3d²4s²), the uncontracted twins Sc-unc/Ti-unc, and
U14e (f³d¹; ECP) — certified end-to-end (the ⁵L ground term at
machine-precision closure).

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
files, tool caches, or artifacts elsewhere. Temporary files live inside
the repo (e.g. a `*.tmp` beside its target, removed once used). Plans live
in-repo as `wayfinder/tickets/` (plan-mode DRAFTS are confined to
`~/.opencode/plan` by the harness and move into the repo as soon as mode
restrictions lift). Mechanical enforcement: the global
`external_directory: ask` permission rule (user-side).

## Running things

```bash
pytest -m fast                                  # the gate's base: 258 tests, seconds
pytest tests/unit/test_X.py tests/integration/Y.py   # + the diff's targeted files
python scripts/run.py <atoms>                   # term tables -> elements/<Atom>/
python scripts/regen_cells.py <atom> \
    [--group grid|spine|planned|lane|lattice|route] \
    [--method casci|casscf|sa-casscf] [--window W] [--mask M] \
    [--start ao|mo|moc|no] [--core raw|avg|frac] [--sector NA-NB] [--full]
python scripts/build_methods_tables.py <atom>    # re-render TABLE.md + TABLE.csv
python scripts/u14e/u14e_vv_ladder.py [--cases TAG ...] [--routes casci,casscf,sa,pinned] [--ground-only]
                                                  # the U14e V&V table's one writer (scripts/u14e/VV.md)
```

- **The per-commit gate** (user-directed): the fast tier plus the diff's
  targeted files, 5–10 minutes maximum. The full-suite run is RETIRED.
  A red alone is never the verdict — isolation reruns are.
- **One regen at a time; no src edits while a chain runs.** Re-render +
  commit per batch/section. Threads: OMP 16–32 default; the large ECP
  legs at 128.
- **Store semantics**: filtered runs MERGE into the store (untouched
  payloads survive); only the unfiltered run replaces it (rebuild =
  grid pass + `--group planned` + `--group lane`). `--full`
  materializes scenario-voided cells — never give it a frac class on an
  atom whose fractional reference does not converge (the head build
  refuses and exits).
- **Per-cell refusals** (`[refused]` in the log) are declared classes —
  the radial mask's construction refusal, `natural.NoCapExceeded`: the
  chain goes on, nothing stored, the pending row holds the place.
- **No calculations without the user's approval.** Strange results get
  reported and revisited, never masked.

## Layout (deep modules, small interfaces)

The interface column is the index; each module's docstring owns the full
contract (purpose, shared conventions, layering — read it before working
in the module).

| Module | Interface | Role |
|---|---|---|
| `system.py` | `AtomConfig`, `build_mol`, `sector_na_nb`, `ncore_of` | validation: single atom, spherical basis; `ecp` pairs basis+pseudopotential; the sector/core helpers raise on parity violations and take `mol.nelectron` (ECP-honest — never recompute from Z) |
| `angular.py` | `l_matrices`, `l_ladder_matrices`, `l2_matrix`, `l_matrices_analytic`, `real_to_complex` | one-electron Lx/Ly/Lz/L² AO matrices — TWO independent constructions certifying each other (libcint `int1e_irp` primary, exact tesseral algebra the cross-check, per-slot through `shell_slots`); `real_to_complex` = the per-shell real→\|l,m⟩ transform |
| `irreps.py` | `IRREP_NAMES`, `IRREP_IDS`, `branch_d2h`, `det_irreps`, `orb_irreps` | the SO(3)→D2h dictionary — which irreps a term's 2L+1 components occupy (a multiplicity can exceed one for L≥2 and is the SA weight / per-solver root count) and which irrep every sector determinant carries; `orb_irreps` = the ONE PySCF symm boundary |
| `detbasis.py` | `determinant_basis`, `sector_dimension`, `one_body_matrix`, `two_body_matrix`, `slater_condon_matrix`, `operator_apply`, `commutator`, `frobenius_norm` | **pure numpy, no pyscf** — Slater–Condon rules over explicit determinants; `slater_condon_matrix` = the fused production path, name-twinned with secondquant; `operator_apply` = the matrix-free twin; ordering ≠ PySCF FCI strings (the bridge lives in operators) |
| `secondquant.py` | `annihilate`, `create`, `apply_one_body`, `apply_two_body`, `one_body_matrix`, `two_body_matrix`, `secondquant_matrix` | **pure numpy, no pyscf** — the same matrices via explicit a+/a application; the elementwise twin of detbasis (random non-particle-symmetric-G pins certify the terms physical integrals are blind to) |
| `operators.py` | `build_operators`, `active_space_integrals`, `full_space_integrals`, `window_closure`, `core_anisotropy`, `ground_block`, `pyscf_ci_permutation`/`ci_in_detbasis_order`, `pyscf_apply` | the first-quantized assembly: H, Lz, L², S², Sz, N, the EXACT parity (exterior power), population labels; owns BOTH PySCF↔ours convention boundaries (the integral transpose + the CI-ordering bridge); the measured-block family (`ground_block` = the SA-anchoring/NO-averaging rule, SCALE-AWARE); core folded into h1 as an m-scalar density; `pyscf_apply` = the Hamiltonian through PySCF's C contractions at any sector size |
| `reference.py` | `ReferenceSpec`, `run_rohf(dm0=, mom_seed=)`, `mom_self_seed`, `run_uhf(dm0=)`, `uhf_seed_density`, `occupation_postcheck` (character-of-occupied), `projected_source`, `rohf_basin_census`, `uhf_basin_census`, `uhf_contamination`, `run_frac_occ(dm0=)`, `run_folded`, `subspace_rohf` | the SCF reference engine: character-forced occupations + basin gate thread every runner; the SEEDED-MOM mode (`seed_basis` — one knob: the projected source carries the configuration, `_mom_occ` follows it, the character post-check verifies; `source_closed/source_open` declare the source basis's own vocabulary; `dm0=` the caller-supplied form — the masters refuse each other; `mom_seed=` the EXPLICIT-seed master: `mom_self_seed` re-assigns the reference's own columns to a declared configuration — majority floor, ε-first); a refused draw's chkfile is deleted — the trusted cache is never poisoned; the UHF seed is `uhf_seed_density`'s spin-split density, hooks unattached — aufbau maintains, the amplified tail sorts last); UHF is a reference-layer product, not a cells axis; `run_folded` = the lane's reference class (MOM-gated, per-atom-gauge pick rule); `subspace_rohf` = the confined-valence certificate (the `run_folded` mirror inverted: valence confined to a distilled window, core exact in its family; MOM assignment on the clean spin-summed Fock — the kernel's Roothaan gauge is mirror-degenerate; Newton arms the certificate legs); module Conventions own caching-trust (cached chks TRUSTED, never re-gated) |
| `casscf.py` | `run_cassci`, `constrained_casscf(mo_full=, mask=, weights=, mc_opts=)`, masks, `run_irrep_casci`/`run_irrep_sa`, `irrep_solvers`, `CASResult` | the ONE deep optimizer driver — a mask is any `uniq_var_indices`-shaped callable; every recipe is a `constrained_casscf` call, never a new kernel path; SA = the two-pass pattern (ensemble kernel → exact roots readout in the frozen window); the symmetry-lane evaluators are pinned-route twins (ONE pinned solver per `branch_d2h` component); `irrep_solvers` = the pin's one construction (the runner's blocks AND the lane readout twin's pricing share it) |
| `windows.py` | `Window(.columns)`, `select_*` (incl. `select_span`, `select_primitives`), `pattern_projector`, `shell_realization`/`realized_patterns`, `member_chains`, `distilled_ncas`/`window_distilled`, `core_shells`/`closed_electrons`/`NOBLE_CORES`, `shell_weights`, `window_avas`, `span_overlap`, `span_projector`, `s_metric_roots`, `shell_slots`, `lowdin(tol=)`, `complete_orthonormal`, `build_full_mo`, `project_orbitals` (the cross-basis projector — the seed carrier), `lm_sectors`, `l_per_row`, `radial_basis`, `mask_radial_rotations` | window construction + selection + the kind vocabulary + the SHELL LAYER (a chemical shell's AO realization DECLARED band-by-band; the identity default = the contracted behavior) + the WINDOW DISTILLATION (k radial directions per shell by the member-summed Gram, replicated on the m-set — m-complete at every k, rotation-closed whatever the capture; the selection is a significance floor, not an ownership test); `build_full_mo` takes AO rows OR column vectors through ONE project-renormalize-assemble treatment; `lowdin` raises on near-linear dependence |
| `core.py` | `purified_core_orbitals`, `confined_core_orbitals`, `closed_shell_core_energy` | the exactly rotation-closed core span (the pur substrate fixing the input-window baseline) + the closed core ON a span (a 1:1 realization; `ncore` follows the PAIRS, never the span width) + the closed-shell core energy bookkeeping constant |
| `basis.py` | `SeriesSpec(open_counts=, seed=)`, `run_series(seed_basis=)`, `valence_block`, `metric_radials`, `series_window` | the BASIS-PREPARATION layer (stage 0b): the configuration SERIES as basin-gated forced-occupation references (the head = the case's own trusted reference, never re-run; a member's `seed="self"` runs MOM-anchored on the declaring head's own columns re-assigned by the per-span `open_counts` — the flat ranking's vocabulary gap closed; `source_closed/source_open` are the seed_basis thread's own and refuse on a self-seed member) + `metric_radials` = the CAPTURE-MAXIMIZING metric eigenproblem (the `(block, s_r⁻¹)` transport) + `series_window` = the tagged windows' construction with certificates; the density-block prior axis beside the incumbent span-signal family — two declared families, never mixed in one comparison row |
| `stages.py` | `StageFacts`, `stage_facts` | ONE analysis pass over a finished stage: the core-field build (the swap-basin detector reads the RAW field below the core-span tolerance) + analyze + energies vs PySCF + window purity + the ortho divergence guard; module Conventions own the two ruler readings (the reconstruction is the ANALYSIS ruler; max\|dE\| decomposes into anisotropy + window definition + bookkeeping — exactly 0) |
| `scans.py` | `ScanSpec`, `scan_window`, `scan_state`, `radial_scan`, `member_tilt_map` | frozen-orbital scans, optimizer-free and DETERMINISTIC given the reference: m-symmetric radial curves + member-tilt maps; `ScanSpec` = the declared per-atom anatomy |
| `analysis.py` | `analyze`, `configuration_summary`, `shell_occupations`, `format_report` | term symbols from ⟨L²⟩/⟨S²⟩ (the half-integer S grid), commutator norms, configuration labels per root (the CI-order bridge runs inside; l-grouping via `lm_sectors`) |
| `perturb.py` | `ClassRow`, `PT2Row`, `readout_mc`, `assert_stored_roots`, `nevpt2_rows`, `nevpt2_rows_on` | the PT2 analyzation at the PySCF-mrpt boundary (the engine: `pyscf.mrpt.NEVPT`, strongly contracted): `readout_mc` = the multi-root CASCI readout rebuilt on the supplied orbitals (the `run_casci` twin without the checkpoint; `ci0=` the STORE-SEEDED form — a quasi-degenerate manifold's plain re-solve re-partitions the roots); `assert_stored_roots` = the fail-fast refusal, the declared bands (energies strict; CI amplitudes 1e-8; the degenerate block's joint projector at p_tol 1e-4, the measured tail sweep; the re-solve angle o_tol 1e-3; the same-state energy noise e_noise 1e-6 riding the vector certificate); `nevpt2_rows` = per declared root, the kernel's canonicalization prelude + the 8 Koopmans classes through mrpt's public functions, the sum pinned against `NEVPT(...).kernel()` (the plumbing-vs-kernel twin: 1e-10 on small CAS, the engine's own ~3e-7 Eh run-to-run wobble floor at campaign size); `nevpt2_rows_on` = the pricing half on a ready `mc` (the lane blocks' entry); pyscf-only imports (substrate objects arrive built), no chkfiles; the reading discipline: README · Analysis; the contract: `wayfinder/specs/pt2-perturbation.md` |
| `cases.py` | `CaseSpec`, `Scenario`, `CAMPAIGN`, `campaign_spec`, `scenarios_for`, `run_case`, `run_atom`, `run_carbon`, `cell_substrate`, `lane_substrate` | the DECLARED campaign (the intent manifest): `CaseSpec` fixes atom/reference/window/sectors/term-target/scenarios + the shell layer (`shell_bands`, `shell_radials` — declared, measured-certified) + the render chemistry; `CAMPAIGN` is the ONE declaration — the cell grid, regen, and renderers all derive from it; new atom = a registry row; `cell_substrate` = the stored cell's live substrate (the post-processing seam's left side: the trusted ROHF rung rebuilt beside the payload); `lane_substrate` = the pinned lane's D2h substrate (basin-steered; the lane runner's and the readout twin's `mol, mf`) |
| `cells.py` | `CellKey`, `run_cell`, `load_cell`, `analyze_cell`, `cell_payload`/`kernel_cell_payload`, `key_core`/`lane_core`/`window_columns` (the key→partition + window-identity resolution the fresh constructions share), `prior_tag_of`/`no_label`, `spine_keys`, `window_grid`/`grid_keys`/`planned`/`lattice_keys`, `lane_keys`, `run_lane_cell`/`load_lane_cell`/`readout_lane`, `void_reason`, `radial_refused` | THE addressing layer: one frozen `CellKey` per grid point, the address grammar as chk name/payload key; VOID cells validated at the key (the scenario RUN GATE — first enabled match; `--full` materializes); the two per-cell construction refusals at their own homes; the pinned lane (the frozen core = the reference's own core columns S-projected); `readout_lane` = the lane store's readout twin (the runner's own re-solve recipe asserted against the record, kernel-ready per-block `mc` objects out); the self-owned store (`ci_stack` stays PySCF order — the bridge lives in operators) |
| `natural.py` | `ci_rdm1`, `natural_orbitals`, `occupation_groups`, `select_no_window`, `no_group_windows`, `NoCapExceeded` | the NO mathematics: spin-summed 1-RDM from CI roots, NOs = density eigenvectors, occupation-GROUP selection never splits an m-set (scale-aware tolerance); symmetry-safe only from a term-AVERAGED root set — the m-unequal density class REFUSED at the closure meter; `NoCapExceeded` = the declared no_cap budget's refusal type (the runner's per-cell skip) |

## Hard-won facts (do not re-derive)

Standing claims + owners; measured values live in the tests, `elements/`,
and git history. Module docstrings own the details — this list is the
cold-start warning surface.

1. **AO labels are truthful but carry member suffixes** (`2px` is really
   px): prefix-match, never regex; fused shells unpack through
   `windows.shell_slots` — never label lists.
2. **Single-root CASSCF breaks rotational symmetry by m-incompleteness of
   the optimized window** (radial member mixtures lacking the full m-set —
   it stays l-pure; l-purity is necessary, not sufficient); SA over the
   degenerate members restores it. Only the commutators measure closure
   itself. Verified end-to-end on the ECP atom (U14e: single-root ⟨L²⟩
   64–66, SA 72.00 exact).
3. **From-scratch runs are not bit-reproducible** — never judge a refactor
   by diffing two from-scratch runs; optimized stages amplify the draw
   chaotically. Cached chkfiles are the reproducibility unit.
4. **Sector derivation + new atoms**: Ms sectors derive from
   `cfg.spin` + active-electron parity — primary `((n+2S)/2, (n−2S)/2)`,
   plus `(n/2, n/2)` when n is even and S > 0. New atom = a registry row
   + an element file, no new code.
5. **Fold/reference classes**: a raw (un-averaged) core field splits
   multiplets in the quadrupole pattern — the averaged core is exact; the
   fractional-occupation reference is a DIFFERENT reference, not the
   ε→0 limit of averaging; SA kernels return only the ensemble energy
   (per-root energies from the roots pass); optimized stages must freeze
   their OWN inactive columns (`build_operators(c_core=)`), else
   E-bookkeeping breaks at O(Eh); in a complete-single-irrep sector
   [H, L²] = 0 trivially — use [H, Lz] and the splitting as
   fingerprints.
6. **An apparent multi-solution may be a construction bug**: a
   NON-orthonormal MO matrix fed to an optimizer diverges below the
   variational bound silently — check C†SC − I of anything fed to an
   optimizer (every stage carries `ortho_err`).
7. **Open-shell campaigns bring new physics, never new code paths**:
   S can be half-integer (the ⟨S²⟩ half-grid snap); the fractional block
   is selected by CHARACTER, never energy order; tests branch on window
   weights, never energies (multi-solution basins incl. window-identity
   swaps are real); SA anchors the member gauge but cannot repair an
   m-incomplete window; raw [noble-gas] cores polarize — commutators
   detect it where ⟨L²⟩ alone does not; quasi-degenerate ground-term
   orderings flip across SCF chaos — never pin them.
8. **The radial constraint axis** (`radial_basis` +
   `mask_radial_rotations`): rotations allowed only between
   same-(l, member) slots — the mask is a PAIR predicate. The member
   gauge of a closure-flat direction is chaotic on single roots; SA
   anchors it reproducibly. An escape not radially reachable is a
   window-identity change, not correlation hunger; a complete START
   window is not a complete CONVERGED window. The Löwdin gauge within a
   degenerate member set is span-preserving, but column identity needs
   member-preserving Gram–Schmidt.
9. **Frozen-scan verdicts** (`scans`; optimizer-free, DETERMINISTIC —
   never pin optimizer output where a scan exists): on p-shell valence
   the full CASSCF gain is a one-parameter m-symmetric radial rotation
   at exact symmetry; open d shells INVERT this — member-asymmetric
   tilts descend (the Jahn–Teller drive behind broken single-root
   solutions; SA averages it away). Sector commutators have BLIND SPOTS
   (m=0 or equal-±pair tilts; maximal-Ms single-irrep sectors) — the
   one-particle projector pair (fact 10) has none.
10. **The one-particle closure meter is the universal certificate**
    (`operators.window_closure`): ‖[P,Lz]‖ and ‖[P,L⊥]‖ of the window
    span projector — both ≈ 0 certifies closure completely. Span overlap
    is the IDENTITY axis, independent of closure — it says nothing
    about it.
11. **Degenerate-shell references are layered chaos**: energy-ordered
    (aufbau) filling collapses from every start while a lower shell is
    open, and at fixed occupations the open shell is MULTI-BASIN —
    character-forced occupations + the basin gate land the intended
    reference; the seed (not the gate) is the cure when gate cycling
    diverges on the basis. Leave real headroom on degenerate-but-
    numerical test bands.
12. **ECP-atom campaign verdicts**: the gated-Newton ROHF reference
    determinant can be a numerically pure L=0 mean field — "below ROHF"
    then means beating a pure state; raw-ROHF MO spans NEVER
    rotation-close (m-delocalized open orbitals); m-complete AO
    lattices win the large pairs on representability + correlation.
13. **Atomic-ROHF CASCI selection rules**: CASCI over a window
    containing the ROHF determinant is variational with EXACT equality
    unless a live coupling exists — the only live channel is the
    same-l double shell (singles Brillouin-dead; same-center
    different-l exchange vanishes angularly). The lattice E column is
    AS-RUN for this reason. NOs-from-ROHF ≡ the MO window (the ROHF
    density is idempotent); meaningful NO sources are term-averaged CI
    densities.
14. **The shell layer and the frozen core**: a chemical shell's AO
    realization is BASIS-SPECIFIC and DECLARED (`CaseSpec.shell_bands`),
    never read off labels; a core's ELECTRON count is chemical
    (`closed_electrons`), never its span width; the PRODUCTION frozen
    core is the reference's OWN core — the occupied columns
    S-projected onto the core span; the purified rotation-closed span
    (`core.py`) fixes the input-window commutator baseline.
15. **The two distillation families** — the incumbent `window_distilled`
    (the reference's span-signal Gram) and the prior axis's
    `basis.series_window` (the occupation-weighted density blocks,
    averaged over the configuration series) — share `member_chains` and
    are NEVER mixed in one comparison row; the density-block radial
    eigenproblem MAXIMIZES the CAPTURE (the `(block, s_r⁻¹)` form; at k
    equal to the chain width the conventions coincide — the span is the
    whole space).

## Conventions

- Per-atom tables live in `elements/<Atom>/TABLE.md` (committed,
  log-provenanced): the atom heading + the ONE methods table between
  the `methods` markers + the `TABLE.csv` mirror (the generator:
  `scripts/build_methods_tables.py`, which refuses atoms without
  TABLE.md); README carries one worked example; cross-atom comparisons
  live in the element files. New atom = new element file + a data row.
- Atomic units, ℏ = 1. Spin-orbitals: alpha = `0..norb-1`, beta = `norb..2norb-1`.
- Determinants: sorted tuples of occupied spin-orbitals.
- The substrate signature `(mol, mf, c_act, ncore, na, nb)` recurs
  through the analysis seam; functions document only their own knobs
  beyond it.
- Every operator enters the determinant space as (one-body matrix,
  two-body integral tensor) — never creation/annihilation algebra.
  Two documented exceptions: the explicit ladder operators
  (`spin_ladder_matrix`) and the many-body parity (`parity_operator`).
- New physics code ships with a test that would catch a sign or phase
  error (spectra, commutators, or independent-construction agreement).
- Definitional slots (module docstrings, function contracts, section
  comments) carry the RULE stated generally — no atom names, no
  measured campaign values; instances live in Examples/Notes, and the
  measured record lives in tests and committed artifacts.

## Tests

`pytest` runs 460 tests. The fast tier: `pytest -m fast` = 366
tests in seconds (machine-light selection). The heavy tail is
HOURS-class — the per-commit gate is the fast tier plus the diff's
targeted files (see Running things); the full-suite run is RETIRED. A
heavy pin runs ISOLATION, on evidence, detached (`setsid`) with
`OMP_NUM_THREADS` capped (8): uncapped, 255 BLAS threads thrash the
box.

**The boundary the suite guards** (user-directed 2026-11-05):
tests verify CONSTRUCTIONS at pinned small sizes — seconds to minutes;
heavy compute belongs to the ATOMIC CALCULATIONS (the campaign runs —
the project's goal), priced per batch and gated at doors. A test that
needs hours is a run in a test costume.

Layout and tiers:

- `tests/unit/` = the seconds tier: pure numpy/libcint, synthetic
  substrates, and the FROZEN carbon ROHF (`tests/carbon_reference.py`,
  loaded by the root `carbon` fixture — NO SCF runs in the unit tier)
  plus its D2h twin (the `carbon_d2h` fixture).
- `tests/unit/test_u14e_vv_*.py` pin the U14e V&V library
  (`scripts/u14e/vv`: the D2h column algebra, the reference substrate,
  the active spaces, the routes and their birth certificate, the
  pricing adapters, the table, the case table) on that frozen carbon
  reference — the library's contract is its modules' docstrings, and
  the pins are the only tests outside `src/`.
- `tests/integration/` = the campaign chemistry (TM references,
  CASSCF/SA, term tables).
- The fast/slow MARKS are the machine axis; the directories the human
  one. Shared DATA and helpers live in plain modules with UNIQUE
  basenames — never import a `conftest` basename directly.
- Test files map 1:1 to src modules; tests are named for the invariant
  they pin (one invariant per test). Each file's own docstring carries
  its content.

Suite-flake protocol: session-scoped fixtures can fail whichever tight
band grazes a chaotic SCF draw (rotating failures across runs, each
isolation-green); free-machine amplification can crash a purified rung
inside PySCF's Davidson while isolated reruns are green — load-dependent
chaos, NOT a regression. **A red suite alone is never the verdict —
isolation reruns are.** Widen only a grazed band, and only with the
observed fresh-draw value in a comment.

### The module walk

The file-by-file walk (ticket 30, completed 2026-10-30; the same
protocol re-ran as the design survey at each wave close) reviewed every
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
   Tolerance envelopes and arrival history are comments at the point
   of use, or commit-message material.
   The standard is concrete: measure against the exemplar passes
   (`system.py`, `windows.py`, `secondquant.py`, `detbasis.py`) --
   one sentence carries one fact, every fact sits in its slot,
   formulas that carry a convention get their own line, and zero
   `Parameters` sections in a walked file is a failed pass.
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
5. **Illegal states become unrepresentable.** Delete the fake degree
   of freedom rather than guard it; guard what the interface can still
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
   kernel -- cross-process bitwise A/B is invalid there, gauge chaos);
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
(open dispositions / rides-Stage-2 / disposition candidates /
watchlist / settled facts), and **Examples** (pointers into `elements/`).
Update it when a decision is made; keep tickets one-session sized.
Resolved tickets are DELETED (git history is the archive); each live
ticket is self-contained (no conversation context needed).

State: the atom campaign (C, Si, Sc, Ti, the unc twins, U14e) is
COMPLETE — U14e certified end-to-end (⁵L at machine-precision closure,
2026-11-16); the ECP wave and the module walk are complete; the docs
round (2026-11-16) consolidated README/AGENTS/map; ticket 37 (the
basis-repair bench, the @s wave, the subspace-ROHF) COMPLETE
(2026-11-12..16 — the seeded-MOM repair, the twins' @s legs, the
per-contraction ΔE certificate). State, ticket 45 (2026-11-20): the
PT2 wave is DONE on all seven atoms — the lane seam (the store's
readout twin), the ladders collected/rendered/gated (the combined
gate ALL GREEN, 1259 roots; the per-atom green logs committed), the
CAS-size flagship read (the corrected ladder PREFERS (4,12) — the
(4,15) promotion double-counts dynamic correlation).  Declared-out
classes recorded in the ticket: the uncs' lane substrate steer, the
semi legs' Löwdin needles, the same-l distilled legs (the future
seam).  Next: the wave review + the Stage-3 adoption
decision (map Disposition candidates) is the user's call; beyond it
the candidates: spin–orbit, transition moments, g-factors.

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
