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

**Filesystem rule (user-set, 2026-10-28)**: using ANY directory outside
this project directory is restricted without the user's explicit
permission — no scratch files, tool caches, or artifacts in `/tmp`, the
home dir, or elsewhere. Temporary files live inside the repo (e.g. a
`*.tmp` beside its target, removed once used). Plans live in-repo as
`wayfinder/tickets/` (plan-mode DRAFTS are confined to
`~/.opencode/plan` by the harness and move into the repo as soon as
mode restrictions lift). Mechanical enforcement: the global
`external_directory: ask` permission rule (user-side, 2026-10-28).

## Layout (deep modules, small interfaces)

| Module | Interface | Notes |
|---|---|---|
| `src/atomterm/system.py` | `AtomConfig`, `build_mol`, `sector_na_nb`, `ncore_of` | validates: single atom, spherical basis (origin by construction); `AtomConfig.ecp` pairs basis+pseudopotential; the sector/core helpers raise on parity violations and take `mol.nelectron` (ECP-honest — never recompute from Z) |
| `src/atomterm/angular.py` | `l_matrices`, `l_ladder_matrices`, `l2_matrix`, `l_matrices_analytic`, `real_to_complex` | one-electron `Lx/Ly/Lz/L²` AO matrices, TWO independent constructions certifying each other (libcint `int1e_irp` primary — the decoding lives in the module docstring; exact tesseral algebra via overlaps the cross-check — per-slot through `windows.shell_slots`, so fused shells work: the whole-shell-block form never ran on Sc's fused p); `real_to_complex` = the per-shell real→\|l,m⟩ transform every M_L bookkeeping needs |
| `src/atomterm/detbasis.py` | `determinant_basis`, `sector_dimension`, `one_body_matrix`, `two_body_matrix`, `slater_condon_matrix`, `commutator`, `frobenius_norm` | **pure numpy, no pyscf** — Slater–Condon rules over explicit determinants; `slater_condon_matrix` = the fused one-pass production path (10^5 dets), `one_body_matrix`/`two_body_matrix` = the parts, name-twinned with secondquant (the cross-certification is name-symmetric); `_slater_condon_slow` = the pinned reference loop; ordering ≠ PySCF FCI strings (the mapping lives in operators; detbasis stays numpy-only) |
| `src/atomterm/secondquant.py` | `annihilate`, `create`, `apply_one_body`, `apply_two_body`, `one_body_matrix`, `two_body_matrix`, `secondquant_matrix` | **pure numpy, no pyscf** — the same matrices via explicit a+/a application on occupation vectors (ket application + scattering, sign tracking); `one_body_matrix`/`two_body_matrix` twin the detbasis names; the random non-particle-symmetric-G pin certifies the string contraction terms that physical integrals are blind to |
| `src/atomterm/operators.py` | `build_operators`, `active_space_integrals`, `full_space_integrals`, `window_closure`, `core_anisotropy`, `ground_block`, `pyscf_ci_permutation`/`ci_in_detbasis_order` | the first-quantized assembly: H, Lz, Lz^2, L^2, S^2, Sz, N + the EXACT parity (`Pi = C^dag(signed S)C` via `l_per_row`; many-body = the exterior power, Cauchy–Binet minors — sharp ±1 on parity-closed windows, deviations are window data, no tolerance) and population labels as `(ndet, ndet)` matrices; owns BOTH PySCF↔ours convention boundaries (the chemist↔physicist integral transpose — ONE, fact 3 — and the CI-ordering bridge `pyscf_ci_permutation`, fact 24) plus the measured-block family (`window_closure`, `core_anisotropy`, `ground_block` = the SA-anchoring / NO-averaging rule, ops-first so no consumer repeats the eigh); core folded into h1 (mcscf convention, spherically averaged l=0 density; `c_core` contract for optimized stages); `method="slater_condon"|"secondquant"` switches the evaluator; `window_closure` = one-particle closure certificate ‖[P,Lz]‖/‖[P,L⊥]‖ of any window span (fact 19); `core_anisotropy` = the ε of the discarded raw core field; the two convention exceptions: `spin_ladder_matrix` (inter-sector a†/a) and `parity_operator` (multiplicative exterior power) — both pin the L²/S² twins / parity exactness |
| `src/atomterm/reference.py` | `ReferenceSpec`, `run_rohf`, `run_frac_occ` | the SCF reference engine: character-forced occupations + basin gate thread `run_rohf(ref=)`/`run_case(ref=)` (and the tests' `run_ladder(ref=)`); the module Conventions own the caching-trust (cached chks TRUSTED, never re-gated) and the S-projector occupation-metric contracts (`span_projector` is the one home); the improper-convergence toolbox lives here (README · Methods) |
| `src/atomterm/casscf.py` | `run_casci`, `constrained_casscf(mo_full=, mask=, weights=, mc_opts=)`, masks, `run_casscf`, `run_sa_casscf`, `CASResult` | the ONE deep optimizer driver — a mask is any `uniq_var_indices`-shaped callable; `run_casscf`/`run_sa_casscf` are the classic pipeline's one-line stage adapters over it, every other recipe a direct `constrained_casscf` call (as cells does — never a new kernel path); mask names form the `mask_` family stating the constraint (`mask_frozen_window`, `mask_no_core_window_mix`, `windows.mask_radial_rotations`) — the serialized axis tokens (frozen/virtmix/radial) are data identifiers mapped once at the cells layer; the module Conventions own the mcscf mutation semantics (fact 6), the full-MO contract (fact 7), and the mask vocabulary; `mc_opts`/`fci_opts` = per-run knobs; all checkpointed (`.chk`) under `atomterm/<name>/`; labels come from exact diagonalization (PySCF ci vectors straggle in degenerate clusters — never a label source); SA runs return roots-pass energies with the SA kernel's CI (the `<name>_roots` store carries the exact per-root CI); `CASResult` carries results only — the stage identity flows through `stage_facts(tag=)`/store keys |
| `src/atomterm/windows.py` | `Window(.columns)`, `SELECTORS`, `window_shell/primitive/scf_mo/no`, `select_*` (incl. `select_span`, `select_primitives`), `primitive_chain`, `shell_weights`, `window_avas`, `span_overlap`, `span_projector`, `s_metric_roots`, `shell_slots`, `lowdin(tol=)`, `complete_orthonormal`, `build_full_mo`, `lm_sectors`, `l_per_row`, `radial_basis`, `mask_radial_rotations` | window construction + selection + the kind registry (the module docstring owns the statement: purpose, the m-completeness contract, the five kinds incl. custom columns, the internal layering — `shell_slots` is the one shell walker [consumed by angular too], the S-metric algebra [`s_metric_roots` = the S^{±1/2} transport] beneath, the two substrate builders); `build_full_mo` takes AO rows OR column vectors through ONE project-renormalize-assemble treatment; `lowdin` raises on near-linear dependence instead of NaN-ing downstream; `shell_weights` groups by `lm_sectors` (no label parsing; row order fixed for bitwise-stable payloads) |
| `src/atomterm/core.py` | `purified_core_orbitals`, `closed_shell_core_energy` | the exactly rotation-closed core span (fact 20's ONE home: the E ⊗ members construction + traps — the pur substrate fixing the fact-17 input-window baseline) and the closed-shell core energy (the pur-fold bookkeeping constant; J/K formula in its docstring) |
| `src/atomterm/stages.py` | `StageFacts`, `stage_facts` | ONE analysis pass over a finished stage: fold (with the stage's own core) + analyze + energies vs PySCF + window purity + ortho divergence guard; the single owner of the bookkeeping pipeline every stage runner shares — callers only render (report text / flat json payload); the module Conventions own the two ruler readings: the reconstruction (eigvalsh(H)+core_energy) is the ANALYSIS ruler (meters/labels computed on it; tables' E columns carry the run's own totals — never compared), and max\|dE\| vs PySCF decomposes into three contributions (core-field anisotropy, window definition, bookkeeping — the last must be 0; on core-overlapping shells the window-definition term dominates, Ti ~85 mEh); `cas.ci` labels the stage and on SA runs is the SA kernel's CI (the constrained_casscf contract) |
| `src/atomterm/scans.py` | `scan_window`, `scan_state`, `radial_scan`, `member_tilt_map` | frozen-orbital scans, optimizer-free: explicit m-symmetric radial curves (exactly closed windows) + member-tilt maps; deterministic given the ROHF reference (fact 18) |
| `src/atomterm/analysis.py` | `analyze`, `configuration_summary`, `shell_occupations`, `format_report` | term symbols from `<L^2>`, `<S^2>` (half-integer S grid; L-letter table to L=20); commutator norms; configuration labels per root — PySCF-order CI in (the fact-24 bridge runs inside), l-grouping via `lm_sectors` (no label parsing) |
| `src/atomterm/cases.py` | `CaseSpec`, `CAMPAIGN`, `run_case`, `run_atom`, `run_carbon` | the DECLARED per-atom campaign (the intent manifest): `CaseSpec` fixes atom/reference-data/window/sectors before anything runs — `planned_windows` declares not-yet-run windows, `lattice_windows` the free-ion d² lattice spans, and the gated reference specs (Ti/Sc character-forced + basin gate) are `CAMPAIGN` data (the ONE declaration: the cell grid, regen, and the renderers all derive from it — regen_cells' REFS/CASES dicts are retired); `run_case(spec)` executes; `run_atom` = the staged implementation |
| `src/atomterm/cells.py` | `CellKey`, `run_cell`, `load_cell`, `analyze_cell`, `cell_payload`, `spine_keys`, `window_grid`/`grid_keys`/`diagnostics`/`planned`/`lattice_keys` | THE addressing layer: one frozen key per grid point (atom/window/sector/method/start/mask/substrate/fold; `substrate` = the core SPAN, `fold` = the core FIELD; pur canonicalizes fold=avg), the address grammar `<atom>_<window>_<na>-<nb>_<method>_<start>_<mask\|->_<substrate>_<fold>` as chk name/payload key; VOID cells validated at the key (radial needs member-labeled AO columns — ao only; mo×pur is VALID since the projection builder — the canonical columns re-projected on the pur core span); the declared key-set derives from `CaseSpec` (`window_grid` owns the one render order; `grid_keys` = the runnable set, `diagnostics` = virtmix, `planned` = the planned windows' grids, `lattice_keys` = the declared d² lattice — renderers branch on key identity only, grouping lives once in the declaration); axes resolve through their OWNING modules; the cell store is self-owned (this PySCF's CASCI kernel writes no chkfile; `ci_stack` stays in PySCF order — the bridge lives in operators, fact 24); SA weights default to the measured ground block; cache hit ⇒ analyze path |
| `src/atomterm/natural.py` | `ci_rdm1`, `natural_orbitals`, `occupation_groups`, `select_no_window`, `no_group_windows` | the NO mathematics: spin-summed 1-RDM from CI roots (secondquant signs), NOs = density eigenvectors in the orthonormal span, occupation-GROUP selection never splits an m-set (scale-aware tolerance); symmetry-safe only from a term-AVERAGED root set (a single \|L,M_L⟩ density is m-diagonal but m-unequal — pinned); `no_group_windows` = the NO arm (group-prefix windows from the term-averaged source density; the block measurement is `operators.ground_block`, the cost stop raises); the Window object + selectors live in `windows.py` |

## Hard-won facts (do not re-derive)

MIGRATION (2026-10-28): fact text is moving into the owning modules'
docstrings as the file-by-file review proceeds; an entry reduced to
`N. <claim> -> <module>` is index-only -- the full text lives in that
module.

1. **AO labels are truthful but carry member suffixes** (`2px` is
   really px); prefix-match labels, never regex patterns; fused
   shells (`nctr > 1`) unpack through the shell walker →
   `windows` (select_shell_ao, shell_slots).
2. **`int1e_irp` decoding** — `T[3a+b] = <mu| r_a d_b |nu>`, `L = -i` ×
   the antisymmetric part → `angular` (module docstring). Pinned by tests.
3. **Two integral conventions collide** — PySCF `int2e`/`ao2mo` are
   chemist, the project tensor is physicist `<pq|rs>` = chemist
   `(pr|qs)`, transposed exactly once at the boundary; the FCI bridge
   passes `g.transpose(0,2,1,3)`; the textbook hand-check is the
   arbiter -> detbasis (the tensor), operators (the boundary),
   test_hamiltonian_matches_pyscf_fci (the arbiter).
4. **Matrix products are exact only in an orthonormal basis** — the
   `S^{±1/2}` transport (`windows.s_metric_roots`); shell AOs of one
   contracted shell are self-orthonormal → `angular.l2_matrix` Notes +
   the windows m-completeness contract.
5. **Core density must be spherically averaged** (`_spherical_average_density`):
   the SCF core carries ~1e-5 numerical p-contamination that would otherwise
   split the ³P multiplet by ~4e-5 Eh. Physical (l=0) component, not a hack.
6. **mcscf objects mutate** -> `casscf` (module Conventions: the
   mutation semantics + the pass-mo-explicitly/fresh-object
   discipline; the SA roots pass is the worked instance).
7. **CASCI/CASSCF take the FULL MO matrix** and slice
   `[ncore:ncore+ncas]` -> `casscf` (module Conventions; a truncated
   matrix fails `nvir >= 0` in `check_sanity`).
8. **Known physics, not bugs**: single-root CASSCF breaks rotational symmetry
   (`<L^2>` ≈ 1.67, `||[H, Lz]||` ≈ 1e-2..1e-1 across reruns); SA-CASSCF over
   the 3 M_L components restores it (⟨L²⟩ = 2.000000). Textbook 2p²: ³P(×3) <
   ¹D(×5) < ¹S(×1), `<L^2>` = 2, 6, 0. Refined (ticket 08): the CASSCF
   breaking is m-incompleteness of the optimized window — it stays l=1-pure
   (p-weight > 0.95), i.e. radial p mixtures lacking the full m-set, NOT s/d
   leakage. l-purity is necessary but not sufficient for closure; only the
   commutators measure closure itself.
9. **Second-quantized two-body convention**: `O2 = ½ Σ G[p,q,r,s] a+_p a+_q
   a_s a_r` with G in the detbasis (physicist) arrangement; annihilators run
   in reverse order (s before r). Elementwise identical to the Slater–Condon
   matrices (tests/test_secondquant.py, tier 1). For core-explicit tests use
   `core.purified_core_orbitals`: canonical core MOs replaced by an
   exactly rotation-closed span (per-l shared radials on all members,
   one global Loewdin) — close to the SCF solution yet exactly l = 0.
   Raw s AOs work equally well;
   unpurified canonical MOs break orthonormality and the [H, Lz] commutator
   at the 1e-4 level (measured contamination: 1s 1.9e-7, 2s 6.3e-6). The
   MO s-weight lives in the AO metric (`c†Sc = 1`), not the Euclidean norm.
   Folded (active-space) H + explicit core energy == the fixed-core block
   of the full-space H (tier 3, machine precision).
10. **From-scratch runs are not bit-reproducible** -> `reference`
   (the module Conventions own the caching/trust contract: cached chkfiles
   TRUSTED, equivalence checks pin the same chkfile -- CASCI-stage reports
   then byte-identical). Residuals: orbital-optimized stages amplify the
   draw chaotically (root splittings wobble ~1e-2 Eh); signed zeros on
   noise-level Lz eigenvalues; never judge a refactor by diffing two
   from-scratch runs.
11. **Generic pipeline** (`cases.run_atom`): Ms sectors derive from
   `cfg.spin = 2S` + active-electron parity — primary `(na, nb) =
   ((n+2S)/2, (n−2S)/2)`, plus `(n/2, n/2)` when n is even and S > 0.
   Sector chkfile tags label M_S = (na−nb)/2 (`casci_ms1`, `casci_ms0`,
   `casci_ms0p5` for a doublet) — na−nb itself would collide. New atoms =
   a new row in `tests/test_atoms.py` `CASES`, no new code.
12. **Carbon rungs facts** (ticket 08; the pins run through cell addresses in test_treatments.py): raw (un-averaged) core
   field → [H, Lz] ≈ 5.6e-6 and the ³P multiplet split 4.4e-5 Eh in the
   {M_L=0} vs {±1} pattern (quadrupole k=2, m=0; dipole excluded by parity);
   core-density anisotropy norm ≈ 2.4e-3 relative. In a sector that is a
   complete single irrep (3-det Ms=1), L² = L(L+1)·Id identically, so
   [H, L²] = 0 trivially — use [H, Lz] and the splitting as fingerprints
   there. Averaged core: exact zeros. Fractional-occupation (average-of-
   configuration) core: invariant by construction but its energy differs
   from averaged-ROHF by 6.6 mEh — average-of-configuration is a DIFFERENT
   reference (open-shell exchange with smeared vs localized valence), not
   the ε→0 limit. Core relaxation (window frozen via `uniq_var_indices`
   mask, inactive↔virtual only) gains −8.8 mEh with exact symmetry.
   SA kernels return only the ensemble energy; per-root energies come from
   a CASCI pass in the converged window (`run_sa_casscf` does this).
   Optimized stages must fold with their OWN inactive columns
   (`build_operators(..., c_core=...)`), else E-bookkeeping breaks at O(Eh).
13. **The Si "multi-solution well" was a construction bug** (ticket 09;
    retires the ticket-07 interpretation): Si's CASSCF/SA did NOT escape to
    an ³S-like basin — `build_full_mo` silently returned a NON-orthonormal
    MO matrix (raw 3p AOs overlap the canonical 2p core MOs by ~0.3, so the
    complement projector `I − P` was inexact and extra junk columns
    appeared) and the optimizer diverged on it (energies below the
    variational bound, C†SC − I ~ 1). With the block-safe construction Si
    behaves exactly like carbon: CASSCF ³P with −0.127 Eh correlation, SA
    restores ⟨L²⟩ = 2.000000. Lesson: check C†SC − I of anything fed to an
    optimizer — every rung carries `ortho_err`.
14. **`build_full_mo` is block-safe** (ticket 09): core block Löwdin'd
    within its own span (span- and density-preserving); active raw AOs
    projected against the core span and renormalized. The correction is
    same-l and m-diagonal (different-l blocks are exactly S-orthogonal),
    so window l-purity and the m-set structure survive. Two consequences:
    (a) naive complement filters (`sv > 1e-8`) admit spurious directions
    when v is not S-orthonormal — silently wrong shapes; (b) on Si,
    max|dE| for CASCI rungs ≈ 9.5e-2 measures the raw-vs-orthogonalized
    window-definition difference, NOT fold anisotropy — max|dE| has up to
    three contributions: core-field anisotropy, window definition,
    bookkeeping (the last must be 0).
15. **Sc / transition metals** (ticket 11): same machinery, new physics.
    (a) Basis: ccpvdz (Balabanov–Peterson) ships for Sc; AO labels are
    clean (3d = exactly 5 AOs, 4s = 1; the fused d shell labels its two
    contractions 4d/5d, so prefix matching has no traps). (b) S can be
    HALF-integer: `spin_quantum_number` snaps S to the half-grid;
    integer rounding relabeled doublets/quartets as singlets/triplets
    with 1e-8 <S^2> noise deciding which (first doublet = Sc). (c) The
    fractional-occupation block is selected by CHARACTER, never energy
    order (on Sc the aufbau block runs away to a p-containing ensemble;
    the metric must be the S-projector -- the bare quadratic form is
    unbounded on diffuse bases) -> `reference` (`run_frac_occ`'s
    `open_idx` + the module Conventions' metric contract). (d) TM rungs
    fingerprints: the raw
    [Ar] core polarizes at the PERCENT level (anisotropy 1.3e-2 rel,
    [H,Lz] 2.7e-2, 2D quintet split 1.9e-3 in PySCF's raw fold); free
    CASSCF from every start collapses to an 2S-like basin +52 mEh with a
    {s,p,d}-mixed window — orthonormal and variational (a real
    multi-solution, NOT the ticket-09 bug class); the shell-start SA
    restores <L^2> = 6.002 but NOT [H,Lz] (0.83): the optimizer swaps the
    window's 4s for a core 3p (window p = 1/6 exactly, core = 4s+5p,
    anisotropy 0.30) — an m-INCOMPLETE window; the canonical-start SA
    collapses outright (start dependence). The relax rung (window frozen)
    is the honest TM rung: sharp 90-state free-ion term table
    (3d1 4s2 + 3d2 4s1 + 3d3), commutators ~1e-3 from core polarization
    alone. (e) New max|dE| class ~0.8 Eh: our averaged fold of a strongly
    polarized SA core vs PySCF's raw one — the Eh-scale version of the
    1e-5 polarized-core detector. (f) The frac-occ reference (+30 mEh vs
    avg) softens the field so the 2D terms of 3d3/3d1 4s2/3d2 4s1 are
    quasi-degenerate at the bottom — their ordering flips across SCF
    chaos (never pin it in tests). (g) Amendment (ticket 14): the
    shell-start SA is BIMODAL across chaos — usually the 4s↔3p
    m-incomplete swap, sometimes the healthy window (commutators then at
    the fact-17 baseline). The sa/relax "≈1e-3" class and the sa_pyscf
    collapse are basins, not deterministic outcomes; tests branch on the
    window weights.
16. **The radial constraint axis** (ticket 14, `radial_basis` +
    `mask_radial_rotations`): rotations allowed only between same-(l, member) slots
    (PySCF's `uniq_var_indices` indexes orbital PAIRS, so the mask is a
    pair predicate; "same-l" alone would NOT preserve closure — general
    same-l rotations mix members). Measured hierarchy: (a) on p² atoms
    the ENTIRE free-CASSCF gain is radially reachable (radial rung ==
    free rung to ≤3e-7 Eh, C and Si) with exact shell purity — the
    m-breaking buys nothing energetically; it is a member-GAUGE freedom
    of a flat direction (closure needs shared radial vectors per member;
    the energy is invariant under member differences). Which gauge a
    single-root optimizer lands in is CHAOTIC: C 5e-11, Si 1.2e-1 (same
    as free), Sc 2e-3..1.6 with the ground term even relabeling ²D→²P —
    never pin single-root radial commutators/L. SA anchors the gauge:
    sa_radial is reproducible to µEh with exact ⟨L²⟩. (b) On Sc the
    escape is NOT radially reachable: the radial ceiling sits ~25 mEh
    above the escaped minimum (stable: +24.4, +25.1) — the escape
    requires l/window-identity change, not radial correlation. (c)
    sa_radial restores ⟨L²⟩=6 but its commutators land at ~1e-3 = the
    input-window baseline (fact 17): SA suppresses the member asymmetry
    (single-root 4e-2 → ~2e-3) but cannot go below the baseline. (d) The
    complete {3d,4s,4p} window (CAS(3,9), 324 dets): CASCI gains ~25 mEh
    of static 4p correlation at EXACT closure, yet free CASSCF still
    escapes (window grabs even f character) — the escape is an optimizer
    attractor, not correlation hunger; a complete START window is not a
    complete CONVERGED window. (e) Construction lessons: the (l, member)
    sectors are the exactly S-orthogonal decomposition for one center,
    but fused shells (nctr > 1) are contraction-outer/m-inner — unpack
    via bas_nctr (label lists miss it); the symmetric Löwdin GAUGE within
    degenerate member sets is an arbitrary rotation — span-preserving,
    which is why CASCI/CASSCF never noticed, but column IDENTITY needs a
    member-preserving Gram–Schmidt (radial_basis); genuine radial
    residuals are ≥8e-2 of S-weight, core tails ≤3e-6, and the
    near-linearly-dependent needles (Sc 6p: 1.7e-4) sit between — keep
    the strong directions labeled, pad the needles unlabeled.
17. **The ~1e-3 CASSCF-rung commutator class is the input window's
    core-orthogonalization defect** — present on the CASSCF start
    object before any optimization; the mechanism and the honest
    re-orthogonalization: `windows.build_full_mo` Notes; the fix (the
    rotation-closed purified span): `core`; the "baseline + optimizer
    gauge" reading of every Si/Sc CASSCF-family commutator:
    README · Methods.
18. **Frozen-scan verdicts** (ticket 16, `scans.py`; optimizer-free,
    DETERMINISTIC — never pin optimizer output where a scan exists):
    (a) p² atoms: the full CASSCF gain is a ONE-parameter m-symmetric
    radial rotation at exact symmetry (C θ*=0.61 over 2p→3p, −250.8 mEh;
    Si θ*=0.85 over 3p→4p, −126.7 mEh; [H,Lz]=0 along the whole curve);
    single-member tilts cost EXACTLY 0 Eh, ≥2 tilted members cost real
    energy (C +1.6 mEh/pair at δ=0.1). (b) Sc/Ti INVERT this: the
    m-symmetric radial curve only rises (nothing to gain anywhere:
    3d→4d, 4s→5s/6s all harmful), member-asymmetric tilts DESCEND (the
    Jahn–Teller class — Sc −0.3/−0.9 mEh, Ti −163 µEh at δ=0.1 — the
    drive behind their chaotically broken single-root solutions; SA
    averages it away). (c) Sector commutators have BLIND SPOTS: [H,Lz]
    cannot see an m=0-member or equal-±pair tilt; multi-term sectors
    catch those via [H,L²], maximal-Ms single-irrep sectors catch
    NOTHING — the projector pair (fact 19) has no blind spots.
    Numbers/curves: README · Motivation; anatomy: elements.
19. **The one-particle closure meter is the universal certificate**
    (ticket 16 P2; `operators.window_closure`, `window_avas`,
    `span_overlap`): ‖[P,Lz]‖ and ‖[P,L⊥]‖ of the window SPAN projector
    (S-metric) close the fact-18c blind spots — [P,L⊥] sees EVERY
    member asymmetry; both ≈ 0 certifies closure completely. Carried by
    every rungs row (|[P,Lz]|, |[P,Lp]|) and StageFacts payloads.
    Span overlap is the IDENTITY axis, independent of closure (Sc's
    AVAS window: 0.9855 overlap yet O(0.5) closure defects — overlap
    says nothing about closure). AVAS verdict: right size,
    mean-field-inherited span; raw shells strictly dominate. SA
    deviations decompose 100% window / 0% core-fold (witness: refold
    the converged core with the raw shell window; committed artifacts
    results/*_sa_decomposition.md -- the one-shot generator retired
    with the rungs, the attribution now rides every cells row as the
    span-certificate + anisotropy columns). Verdict tables: README ·
    Results (Methods compared).
20. **The closed core span** — E_l ⊗ members (per l: top-k span of
    the S-metric COLUMN Gram of the member-resolved radial blocks,
    every radial replicated on ALL members, one global Löwdin) makes
    both input windows close at machine precision and the core density
    exactly m-scalar -> `core` (`purified_core_orbitals` owns the
    construction + the traps). Rung phenomenology on this substrate
    (pur rungs' machine-exact closure incl. Ti, the exactly FLAT
    member gauge on Sc/Ti and sa_radial_pur's floor, the 4s↔3p swap
    attractor surviving purification — branch on window weights, the
    window-definition offset vs raw-fold CASCI Sc +88.7 mEh,
    relax-family max|dE| ~1e-3 fold bookkeeping — bound, don't pin):
    README · Methods, · Results; the rungs pins: `tests/integration/test_treatments.py` (through cells).

21. **Ti: d², the aufbau runaway, and the basin gate** (ticket 12):
    (a) THE REFERENCE IS TWO-LAYERED CHAOS -- the aufbau runaway
    (+123.07 mEh, a doubly-occupied 3d inside the canonical core;
    energy-ordered filling falls back from EVERY start) and, at fixed
    occupations, a MULTI-BASIN degenerate open shell (+11.2 mEh
    metastable basin, same occupations, guess-dependent) -> `reference`
    (the character-forced hook + the `occ_e_ref` gate land
    -848.40648875 with pure-3d open MOs; the gate also fixes the relax
    rung's substrate: good-basin input defect ~1e-8, casscf_relax at
    2.3e-8 -- not Sc's fact-17 class, that class is basin-specific).
    (b) Tables (pinned in tests): option A CAS(2,5)/{3d}
    — Ms=1 ³F×7+³P×3; Ms=0 ¹S+¹D+³P+¹G×9+³F in NIST order ³F<¹D<³P<¹G<¹S
    (the `ms0_order` pin); first L=3 ground, first ¹G, first f-bearing
    basis (core f-tails ~1e-4). Option B CAS(4,6): TI_EXACT_TERMS in
    conftest (120 dets, S≥1 only). Raw-core caveat: Ti's ground ⟨L²⟩ is
    UNSHIFTED at the percent-class commutators (4.9e-2/2.6e-1) — ⟨L²⟩
    alone is not a raw-core detector on every atom; commutators are.
    (c) The Sc phenomenology transfers (escape, the exact 4s↔3p swap,
    the partial swap under pur — branch on window weights); the radial
    anatomy is a THIRD kind: raw valence rotations all RISE, the ~11 mEh
    radial gain lives on the ORTHOGONIZED substrate's s slot
    (converged sa_radial window: 4s 0.75 + 6s 0.27); never rotate raw
    valence AOs to find it, and a rotated-but-nonorthogonalized column
    is the ticket-09 bug class. d² tilts descend (JT class) so
    sa_radial_pur stops above its Sc floor (1.4e-4/5.3e-4). pur rungs
    machine-exact (6.2e-14). (d) Test-band lessons: leave real headroom
    on degenerate-but-numerical quantities (the frac septet 6.5e-11 vs
    a 1e-10 bound flaked once) and on basin-mixed ⟨L²⟩. Narrative:
    README · Methods, · Results, · Motivation.


22. **U campaign: matched cost, the ⁵S° reference, and the window
    verdicts** (tickets 21/24): (a) the gated Newton ROHF reference
    (−50.38185233) determinant is a numerically pure ⁵S° state —
    ⟨L²⟩ = 0 in three independent constructions, all 3,876 generic
    dets of 1f+1d sit at ⟨L²⟩ ≥ 6; "below ROHF" means beating an ⁵S°
    mean field. Mechanism open (AOC-smearing hypothesis). (b) The f²d²
    ground blocks are ⁵L (L = 8, ⟨L²⟩ = 72) — the 56 labels that
    circulated in conversation were wrong; logs are authoritative.
    (c) Matched-cost verdicts: raw-ROHF MO windows win small pairs by
    REPRESENTABILITY (the reference det's projection onto 3f+2d has
    norm 0.948; every AO window below 4f+2d sits above the mean
    field), then saturate correlation-starved (−2.6 mEh); the
    m-complete AO lattice wins the top pairs (Sc 1,274 dets by
    19.7 mEh, Ti 5,096 by 18.9) and the 4f+2d corner resolves the
    Δf/Δd race (⁵L° by > 37 mEh). (d) Raw-ROHF MO spans NEVER
    rotation-close: ‖[P,L⊥]‖ = 7–9 at every size 4..100 (the SCF's
    open orbitals are m-delocalized; each added virtual adds
    fractional m-directions) — the certificate excludes symmetric
    roots from raw-MO CASCI without running any CI. (e) Both TMs
    cross below their escaped CASSCF basins at exact symmetry
    (fact-16b ceilings are window statements).
23. **The γ ladder and NO windows** -> natural (the module docstring
   owns the m-scalar criterion — term-averaged γ = R ⊗ c·𝟙_{2l+1} ⇒
   whole m-sets; single-root γ is m-diagonal but m-UNEQUAL — and the
   scale-aware grouping tolerance with the defect-leakage numbers;
   measured campaign scale: U groups d/f/d/f = 5/7/5/7, every group
   window ~1e-13 certificate, 12 NO cols within 4.6 mEh of the full
   window) — index-only.

24. **PySCF FCI ordering ≠ detbasis ordering** -> operators (the
   bridge + the Conventions statement; certified by the
   eigen-equation pins in test_cells) + cells (the store contract:
   `ci_stack` stays PySCF order; this PySCF's CASCI kernel writes no
   chkfile, the cell layer owns its complete store) — index-only.

25. **The atomic-ROHF CASCI selection rules** (the 3d² lattice,
    2026-10-27): CASCI over a window containing the ROHF determinant
    is variational (≤ E_ROHF) with EXACT equality unless a live
    coupling exists — singles are Brillouin-dead (stationarity),
    doubles are dead unless same-l (the same-center different-l
    EXCHANGE integral vanishes angularly; parity additionally kills
    d¹p¹↔d²), so the only live channel is the same-l double shell.
    Measured as-run (Ti 3d² Ms=1): MO windows {3d}/{3d,4p} +0.00 µEh
    (the ROHF determinant IS the CI ground state), {3d,4p,4d}
    −39.66 µEh (the 4d channel opens); AO windows are unbounded
    (+6.44 mEh raw-shell offset; the −3.14 mEh crossing = the d′
    channel minus the offset). The lattice E column is AS-RUN for
    this reason — the averaged-fold reconstruction sits +2.8 mEh
    above every as-run Ti CASCI total and bends every anchor
    comparison by exactly that; it stays the ANALYSIS fold (the
    meters are computed on it). NO-source corollary: the ROHF density
    is idempotent, so NOs-from-ROHF ≡ the MO window (occupations
    exactly {0,1,2}); fractional occupations need the term-averaged
    CASCI density (fact 23) or the frac-occ SCF (uncorrelated → no
    tail groups).
## Conventions

- Per-atom tables live in `elements/<Atom>/TABLES.md` (committed,
  log-provenanced); README · Results carries the Ti worked example; cross-atom
  comparisons live in the element files. New atom = new element file + a
  data row, not
  README growth.
- Atomic units, ℏ = 1. Spin-orbitals: alpha = `0..norb-1`, beta = `norb..2norb-1`.
- Determinants: sorted tuples of occupied spin-orbitals.
- Every operator enters the determinant space as (one-body matrix,
  two-body integral tensor) — never creation/annihilation algebra.
  Two documented exceptions: the explicit ladder operators
  (`spin_ladder_matrix`) apply a†/a directly — an inter-sector map has
  no sector-symbol representation; its signs are pinned by the
  twin-construction tests (L²/S² from L±/S±) — and the many-body
  parity (`parity_operator`) is multiplicative — the exterior power
  built as Cauchy–Binet minors over determinant pairs, pinned by the
  hand-derived mixed-window matrix and the involution on closed spans.
- New physics code ships with a test that would catch a sign or phase error
  (spectra, commutators, or independent-construction agreement).

## Tests

`pytest` runs 184 tests in ~13 min unloaded (chkfile caching; ROHF computed
once per
session; the Ti rungs block is the slowest, ~7 min cold; the Sc rungs
~2.5 min; under external load ~230 the same suite took 31 min, 2026-10-28).
Layout: `tests/unit/` = the seconds tier (pure numpy/libcint,
synthetic substrates, and the FROZEN carbon ROHF --
`tests/carbon_reference.py`, repr-round-trip-exact source, loaded by
the root `carbon` fixture: NO SCF runs in the unit tier and the
downstream is bit-stable across runs, certified by double-run
sha256; regenerate only on basis/default changes via
`scripts/generate_fixture_carbon.py`);
`tests/integration/` = the campaign chemistry (TM references, CASSCF/SA,
term tables); conftest is split the same way (shared root, campaign
fixtures in `integration/conftest.py`). Shared DATA and helpers live in
plain modules with UNIQUE basenames (`integration/refs.py`,
`unit/helpers.py`) -- never import the `conftest` basename directly:
pytest claims it per-directory, and direct imports resolve
order-dependently across tiers (the 2026-10-28 suite catch). Tiers:
`pytest -m fast` = the
machine-light selection (the unit tier plus the frozen-reference
integration files, e.g. test_stages); `pytest -m slow` = the
SCF-bearing chemistry (the
fast/slow MARKS are the machine axis, the directories the human one;
test_cells mixes tiers per function as stated policy -- 1:1 with src,
its module docstring says so; every other file marks at module level).
Run the full suite
before committing changes to
physics-bearing modules (operators/stages/cells/casscf/scans); targeted
files or `-m fast` are fine for everything else. Suite flakes: the rungs
fixtures are session-scoped — a chaotic SCF draw can fail whichever tight
band grazes it (rotating failures across runs, each isolation-green on
redraw); widen only a grazed band, and only with the observed fresh-draw
value in a comment (precedents: sa_pur 2.2e-2 2026-10-27; the Ti escape
gauge 6.40 and the Sc partial-swap 0.0736, 2026-10-29). FREE-MACHINE
AMPLIFICATION (2026-10-29): with the machine unloaded (all BLAS threads
live), full-suite sessions intermittently crash the C/Si purified rung
inside PySCF's Davidson (NaN CI overlaps; one fixture death cascades to
~12 errors) while every isolated rerun is green — load-dependent chaos,
NOT a regression; a red suite alone is never the verdict, isolation
reruns are (measured: 2 red + 1 near-green suites on three identical
consecutive runs). Test files map
1:1 to src modules (`test_substrate_chemistry.py` carries the quartet's
SCF-bearing chemistry over reference/windows/core/casscf, imported
OWNER-DIRECT since the phase-2 walk -- the facade retired in phase 5;
`test_reference.py` + `test_windows.py` are the
unit tiers of reference and windows -- the latter: the selection
layer, the window-columns generalization, the
`Window.columns` one-treatment collapse, and the AO-rows ≡ AO-columns
adapter, on synthetic + frozen-carbon substrates, plus the STO-3G
degenerate cases:
null core, the whole AO space as window -- certificates exact by
construction); tests are named for the
invariant they pin. The critical invariants: `test_cells.py`
(the cell layer: address↔key round-trip + the pur-fold
canonicalization; THE CI-ordering bridge certified by the
eigen-equation against an independently built Slater–Condon H on five
sectors, plus the identity-permutation pin of the carbon coincidence
-- fact 24; the store round-trip incl. the SA `_roots` preference;
run_cell ≡ the direct registry construction, analyze_cell(stored) ≡
in-session -- the store certificate -- and post-hoc gauge-safe
observables from the stored CI, all on carbon; the mo×pur key valid +
the declaration groups: `window_grid` = the one render order, the
spine a subset, `diagnostics`/`planned` the filters, `lattice_keys`
the d²-lattice derivation),
`test_hamiltonian_matches_pyscf_fci`
(Slater–Condon vs PySCF FCI), `test_elementwise_identical_to_slater_condon`
(second quantization ≡ first quantization, elementwise) +
`test_parts_match_slater_condon_on_random_g` (the name-twinned parts on
random complex non-particle-symmetric G — the string contraction terms
physical integrals are blind to) +
`test_rhf_closed_shell_bridge_be_sto3g` (closed-shell RHF bridge, S=0) +
`test_vectorized_matches_slow_loop_multichunk` (the chunk-offset path:
the production memory budget yields ONE chunk at test sizes, so the
budget is monkeypatched),
`test_full_space_vs_pyscf_casci_core_explicit` (core-explicit CAS(6,5) vs
PySCF ncore=0), `test_full_multiplet_ms0` (textbook ¹S/¹D/³P),
`test_commutators_vanish_in_shell_space` (exact symmetry),
`test_analytic_matches_libcint` (two independent L constructions),
`test_orbital_ladders_twin_construct_l2` +
`test_spin_ladders_twin_construct_s2_across_sectors` (the ladder
operators built explicitly: L² = Lz² + ½(L₊L₋ + L₋L₊) and S² = Sz² +
½(S₋S₊ + S₊S₋) elementwise on every p² sector, S± as inter-sector
maps),
`test_atom_term_tables[*]` (generic pipeline per atom; new atoms = new rows),
`test_treatments.py` (the treatment ladder through CELL ADDRESSES: per-rung fingerprints, parameterized C/Si; Sc TM block:
exact closure + free-ion tables, relax honesty, S-basin escape, the bimodal
SA basins, radial ceiling vs escape, complete-window still-escapes, and the
purified-substrate block: machine-exact closure on C/Si, Sc's
`sa_relax_pur`/`relax_pur` at 1e-13 with unbroken Lz/L², the flat-gauge
floor, the surviving swap attractor — branched on window weights; the Ti
block: TI_EXACT_TERMS, the gated-basin relax honesty, pur machine-exactness,
the escape + swap basins, sa_radial's orthogonalized-substrate gain),
`test_purified_core_orbitals_exact_closed_span` +
`test_purified_substrate_closes_the_input_windows` (fact 20: the closure
certificate on the core span and both input windows, m-scalar density,
contamination/density-shift scales — the gate-(a)/(b) pins),
`test_ti_rohf_reference_needs_character_forced_occupations` (fact 21a:
the runaway + the hook + the basin gate's reference numbers),
`test_radial_rung_full_gain_without_l_leakage` (fact 16a: radial == free
gain, zero l-leakage), `test_divergence_guard_fires_on_broken_orthonormality`
(the ticket-09 guard's negative control; in test_stages, the fast tier),
`test_casscf_input_window_orthogonalization_defect` (fact 17),
`test_radial_basis_labels_and_partners` + `test_lm_sectors_partition_fused_shells`
(the labeled substrate), `test_avas_cross_validates_the_shell_window`
(independent ecosystem window selection), `test_window_closure_certificate_and_blind_spots`
and `test_window_closure_sc_shell_and_defective` (the one-particle meter:
exact certificates, the fact-18c blind spots closed by [P,L⊥]),
`test_window_avas_reproduces_the_shell_cascade` (fact 19: AVAS inherits
the mean-field relaxation — C = the full θ* gain near-symmetric, Sc = +80 mEh
and O(0.5) closure defects), `test_stages.py` (the facts-module
contract), `test_scans.py` (frozen-scan anatomy: full radial gain at exact
symmetry, exactly-flat single-member tilts, the Lz blind spots, the Sc
inversion, and the Ti d² verdict — tilts descend, radial rises, the
Lz-blind ±pair pattern — deterministic, no optimizer),
`test_half_integer_spin_labels` (doublet
labeling), `test_natural.py` (the γ-ladder claim: m-scalar densities yield whole m-sets, single roots do not; the 1-RDM pinned against a brute-force pair sum; real-AO certificates), `test_substrate_chemistry.py` (build_full_mo block properties — the
ticket-09 regression guard — the frac-occ ordering; the
selection/meter pins incl. shell_weights partition live in the unit
`test_windows.py`),
`test_analysis.py` (configuration labels, the PySCF-order bridge pin,
the complex-H guard), `test_system.py` (validation).

### The module walk

The review proceeds module by module: read the module and its test
file together; restructure per Refactor; certify per Verification;
collapse any facts whose homes the module now owns to index lines;
update the layout row in the same commit.

## Refactor

The rules every restructuring of src/ follows (formulated 2026-10-28;
refined as the review walk proceeds; certifying the result is the
Verification section's job):

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
5. **Illegal states become unrepresentable.** Delete the fake degree
   of freedom rather than guard it; guard what the interface can
   still express, failing fast at the seam that owns the invariant,
   the message stating why.
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

The rules every change to physics-bearing code is verified by
(formulated 2026-10-28; refined as the review walk proceeds). The
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
   (deterministic constructions only); elementwise pins for algorithm
   restructures; AST comparison with docstrings neutralized for
   docstring-only passes (`scripts/check_docstring_only.py`); the
   full suite before physics-bearing commits.
   Mathematically idempotent steps still count -- they perturb
   degenerate gauges and break bit-exactness.
7. **Measured pins.** Bands cite the observed fresh value in a
   comment; docstring examples run once before being written; hand
   derivations defer to the machine certification.

A verification pass ends the same way: rule by rule, applied or
explicitly declined.

## Planning

`wayfinder/map.md` holds the effort map (destination, decisions, fog) and
lists the live tickets. Update it when a decision is made; keep tickets
one-session sized. The atom campaign (C, Si, Sc, Ti) is COMPLETE — tickets
01–16 all RESOLVED (facts 18–21; Ti delivered the d² tables incl. ¹G, the
character-forced/basin-gated ROHF reference, the JT verdict for ³F, and the
machine-exact pur rungs on a second TM). Next wave: ECP (design pre-verified
in map.md "Post-wave survey"); beyond it the user's candidate list:
spin–orbit, transition moments, g-factors; each ticket is self-contained in
`wayfinder/tickets/` (no conversation context needed).

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
  (candidates, declinations, and why) in map.md under "Post-wave
  survey" with the date. Widen a grazed test band only with the
  observed fresh-draw value in a comment (2026-10-15 precedent).
