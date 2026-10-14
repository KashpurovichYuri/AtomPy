# AGENTS.md — AtomTerm

Working agreement for coding agents on this repository.

## What this project is

**AtomTerm** classifies electronic states of atoms by their atomic term
symbols (`^(2S+1)L`) using *exact* angular momentum operators. Pipeline:
ROHF → CASCI → CASSCF → SA-CASSCF (PySCF), then a first-quantized analysis:
`Lz`, `L^2`, `S^2`, `N` built as plain matrices in the Slater-determinant
basis (no second quantization anywhere), symmetry quantified by
`||[H, Lz]||_F` and `||[H, L^2]||_F`.

## Environment

```bash
source ~/miniconda3/etc/profile.d/conda.sh && conda activate pyscf
# or call ~/miniconda3/envs/pyscf/bin/python directly
pytest                      # from repo root; pythonpath=src is preconfigured
python scripts/run_carbon.py [workdir]
```

Python 3.12, PySCF 2.13, numpy 2.5, scipy 1.18 (note: `scipy.special.sph_harm`
is gone; use `sph_harm_y(n, m, polar, azim)`).

## Layout (deep modules, small interfaces)

| Module | Interface | Notes |
|---|---|---|
| `src/atomterm/system.py` | `AtomConfig`, `build_mol` | validates: single atom at origin, spherical basis |
| `src/atomterm/angular.py` | `l_matrices`, `l2_matrix`, `l_matrices_analytic` | one-electron `Lx/Ly/Lz` AO matrices. Primary: libcint `int1e_irp` (decoded layout in the docstring). Cross-check: exact tesseral algebra via overlaps. |
| `src/atomterm/detbasis.py` | `determinant_basis`, `slater_condon_matrix`, `commutator` | **pure numpy, no pyscf** — Slater–Condon rules over explicit determinants |
| `src/atomterm/secondquant.py` | `annihilate`, `create`, `secondquant_matrix` | **pure numpy, no pyscf** — same matrices via explicit a+/a application on occupation vectors (correct operator application, sign tracking) |
| `src/atomterm/operators.py` | `build_operators`, `full_space_integrals` | assembles H, Lz, L^2, S^2, N as full `(ndet, ndet)` matrices; core folded into h1 (mcscf convention, spherically averaged); `method="slater_condon"|"secondquant"` switches the many-body evaluation |
| `src/atomterm/drivers.py` | `run_rohf`, `run_frac_occ_scf`, `run_casci`, `run_casscf`, `run_sa_casscf`, `run_casscf_core_relax`, `run_casscf_pyscf_start`, `select_valence_ao`, `select_shell_ao`, `shell_weights`, `purified_s_core` | all checkpointed (`.chk`); active space = raw shell AO columns placed at `[ncore:ncore+ncas]`; core-explicit tests use `purified_s_core` |
| `src/atomterm/analysis.py` | `analyze`, `format_report` | term symbols from `<L^2>`, `<S^2>`; commutator norms |
| `src/atomterm/pipeline.py` | `run_atom`, `run_carbon` | generic end-to-end orchestration for any atom (Ms sectors derived from spin + active-electron count); `run_carbon` = carbon wrapper; writes `results/` |
| `src/atomterm/ladder.py` | `run_carbon_ladder` | ticket-08 rung comparison (raw/avg/frac-occ core, core relax, CASSCF starts, SA), all-vs-all matrix |

## Hard-won facts (do not re-derive)

1. **PySCF AO labels are truthful** (`2px` really is px); shells may be fused
   (`nctr > 1`): cc-pVDZ C shell 0 is one contracted s shell with two AOs.
   Do not use `mol.search_ao_label` — patterns are regex-flavored (`"2p*"`
   matches `"2s"`); prefix-match the labels instead.
2. **`int1e_irp`** (libcint, 9 components) is `T[3a+b] = <mu| r_a d_b |nu>`;
   `Lz = -i(T1 - T3)`, `Lx = -i(T5 - T7)`, `Ly = -i(T6 - T2)`. Pinned by tests.
3. **Two integral conventions collide**: PySCF `int2e`/`ao2mo` are **chemist**
   `(pq|rs) = int p(1) r(1) q(2) s(2)`; `detbasis` uses **physicist**
   `<pq|rs> = int p(1) q(2) r(1) s(2)` = chemist `(pr|qs)`.
   `active_space_integrals` returns physicist (transpose(0,2,1,3) applied).
   `fci.direct_spin1.kernel` 4-fold input is chemist; pass `g.transpose(0,2,1,3)`.
   The textbook hand-check `E(3P det) = 2 h_2p + J - K` is the arbiter.
4. **Matrix products are exact only in an orthonormal basis** (`l2_matrix`
   needs the `S^{±1/2}` metric dance). Shell AOs of one contracted shell are
   self-orthonormal, so the active space is exactly rotation-closed.
5. **Core density must be spherically averaged** (`_spherical_average_density`):
   the SCF core carries ~1e-5 numerical p-contamination that would otherwise
   split the ³P multiplet by ~4e-5 Eh. Physical (l=0) component, not a hack.
6. **mcscf objects mutate**: `mc.kernel()` canonicalizes `mc.mo_coeff`;
   `get_h1eff()/get_h2eff()` re-process `mc.mo_coeff` — pass `mo` explicitly
   and/or use a fresh object for integrals.
7. **CASCI/CASSCF take the FULL MO matrix** and slice `[ncore:ncore+ncas]`;
   a truncated matrix fails `nvir >= 0` in `check_sanity`.
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
   `drivers.purified_s_core`: canonical core MOs projected onto the s-AO
   subspace (dropping the SCF p/d tail) and Löwdin-reorthonormalized —
   close to the SCF solution yet exactly l = 0. Raw s AOs work equally well;
   unpurified canonical MOs break orthonormality and the [H, Lz] commutator
   at the 1e-4 level (measured contamination: 1s 1.9e-7, 2s 6.3e-6). The
   MO s-weight lives in the AO metric (`c†Sc = 1`), not the Euclidean norm.
   Folded (active-space) H + explicit core energy == the fixed-core block
   of the full-space H (tier 3, machine precision).
10. **From-scratch runs are not bit-reproducible**: threaded BLAS makes the
   SCF converge to last-digit-different MOs, and orbital-optimized stages
   (CASSCF/SA) amplify that chaotically (CASSCF root splittings wobble
   ~1e-2 Eh between identical reruns). Behavior-equivalence checks MUST pin
   the same ROHF chkfile as input — CASCI-stage reports are then
   byte-identical (residual: signed zeros `0.0`/-`0.0` on noise-level Lz
   eigenvalues). Never judge a refactor by diffing two from-scratch runs.
11. **Generic pipeline** (`pipeline.run_atom`): Ms sectors derive from
   `cfg.spin = 2S` + active-electron parity — primary `(na, nb) =
   ((n+2S)/2, (n−2S)/2)`, plus `(n/2, n/2)` when n is even and S > 0.
   Sector chkfile tags label M_S = (na−nb)/2 (`casci_ms1`, `casci_ms0`,
   `casci_ms0p5` for a doublet) — na−nb itself would collide. New atoms =
   a new row in `tests/test_atoms.py` `CASES`, no new code.
12. **Carbon ladder facts** (ticket 08, `ladder.py`): raw (un-averaged) core
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
13. **CASSCF multi-solution is not TM-only**: for Si 3s²3p² the single-root
   CASSCF optimizer escapes the 3p shell (pulls 3s in, lands in an
   ³S-like well) and the ³P-weighted SA collapses too — so orbital
   optimization beyond carbon needs a stability strategy (CASCI-seeded
   start or frozen inter-shell rotations); Si/TM rows run CASCI only.
   Also: folded-core vs PySCF max|dE| scales with core size (Si 3.1e-4 Eh
   vs C ~1e-7) — same spherical-average cause as fact 5, harmless.

## Conventions

- Atomic units, ℏ = 1. Spin-orbitals: alpha = `0..norb-1`, beta = `norb..2norb-1`.
- Determinants: sorted tuples of occupied spin-orbitals.
- Every operator enters the determinant space as (one-body matrix,
  two-body integral tensor) — never creation/annihilation algebra.
- New physics code ships with a test that would catch a sign or phase error
  (spectra, commutators, or independent-construction agreement).

## Tests

`pytest` runs 40 tests in ~15 s (chkfile caching; ROHF computed once per
session). The critical invariants: `test_hamiltonian_matches_pyscf_fci`
(Slater–Condon vs PySCF FCI), `test_elementwise_identical_to_slater_condon`
(second quantization ≡ first quantization, elementwise),
`test_full_space_vs_pyscf_casci_core_explicit` (core-explicit CAS(6,5) vs
PySCF ncore=0), `test_full_multiplet_ms0` (textbook ¹S/¹D/³P),
`test_commutators_vanish_in_shell_space` (exact symmetry),
`test_analytic_matches_libcint` (two independent L constructions),
`test_atom_term_tables[*]` (generic pipeline per atom; new atoms = new rows),
`test_ladder.py` (per-rung symmetry fingerprints of ticket 08).

## Planning

`wayfinder/map.md` holds the effort map (destination, decisions, fog) and
lists the live tickets. Update it when a decision is made; keep tickets
one-session sized. The live wave (08: carbon ladder, 09: Sc, 10: Ti, then 11 optional, 12
review protocol) is fully specified in `wayfinder/tickets/` — each ticket
is self-contained (no conversation context needed). Review protocol:
`/code-review` (two axes, parallel sub-agents) after each ticket wave, then
`/improve-codebase-architecture` survey; see ticket 12.
