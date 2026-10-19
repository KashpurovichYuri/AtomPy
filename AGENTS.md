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
python scripts/run_c.py [workdir]
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
| `src/atomterm/operators.py` | `build_operators`, `full_space_integrals`, `window_closure` | assembles H, Lz, L^2, S^2, N as full `(ndet, ndet)` matrices; core folded into h1 (mcscf convention, spherically averaged); `method="slater_condon"|"secondquant"` switches the many-body evaluation; `window_closure` = one-particle closure certificate ‖[P,Lz]‖/‖[P,L⊥]‖ of any window span (fact 19) |
| `src/atomterm/drivers.py` | `run_rohf`, `run_frac_occ_scf`, `run_casci`, `constrained_casscf`, `core_virtual_only_mask`, `no_core_window_exchange_mask`, `lm_sectors`, `radial_basis`, `radial_mask`, `run_sa_casscf`, `run_casscf_*`, `select_valence_ao`, `select_shell_ao`, `shell_weights`, `purified_s_core`, `purified_core_orbitals`, `avas_window`, `span_overlap` | all checkpointed (`.chk`); `constrained_casscf(mo_full=, mask=, weights=)` is the one deep CASSCF driver — a mask is any `uniq_var_indices`-shaped callable, and the `run_casscf_*` names are one-line adapters over it; `lm_sectors`/`radial_basis(c_core=)`/`radial_mask` build the (l, member)-labeled substrate for same-sector (radial) constraints (facts 16e/17); `purified_core_orbitals` = the exactly L-closed core span (fact 20: shared radial per closed shell, E ⊗ members closure theorem); active space = raw shell AO columns placed at `[ncore:ncore+ncas]`; core-explicit tests use `purified_s_core`; `avas_window`/`span_overlap` = the ecosystem window probe (fact 19) |
| `src/atomterm/stages.py` | `StageFacts`, `stage_facts` | ONE analysis pass over a finished stage: fold (with the stage's own core) + analyze + energies vs PySCF + window purity + ortho divergence guard; the single owner of the bookkeeping pipeline and ladder share — callers only render (report text / flat json payload) |
| `src/atomterm/scans.py` | `scan_window`, `scan_state`, `radial_scan`, `member_tilt_map` | frozen-orbital scans, optimizer-free: explicit m-symmetric radial curves (exactly closed windows) + member-tilt maps; deterministic given the ROHF reference (fact 18) |
| `src/atomterm/analysis.py` | `analyze`, `configuration_summary`, `format_report` | term symbols from `<L^2>`, `<S^2>`; commutator norms; approximate configuration labels per root |
| `src/atomterm/pipeline.py` | `run_atom`, `run_carbon` | generic end-to-end orchestration for any atom (Ms sectors derived from spin + active-electron count); `run_carbon` = carbon wrapper; writes `results/` |
| `src/atomterm/ladder.py` | `Rung`, `run_ladder`, `run_carbon_ladder` | orbital-treatment rungs as a DATA registry (`_registry` returns `list[Rung]`): new rung = new entry, report iterates the registry (cannot drop a rung); `rung_tags=` selects a subset (window studies); `purified=True` appends the `*_pur` rungs (same treatments on the exactly-closed purified core substrate, fact 20 — incl. `sa_relax_pur`, the exact-symmetry SA rung); all-vs-all matrix + per-rung comparison columns (win_l purity, core anisotropy), ortho_err divergence guard |

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
13. **The Si "multi-solution well" was a construction bug** (ticket 09;
    retires the ticket-07 interpretation): Si's CASSCF/SA did NOT escape to
    an ³S-like basin — `build_full_mo` silently returned a NON-orthonormal
    MO matrix (raw 3p AOs overlap the canonical 2p core MOs by ~0.3, so the
    complement projector `I − P` was inexact and extra junk columns
    appeared) and the optimizer diverged on it (energies below the
    variational bound, C†SC − I ~ 1). With the block-safe construction Si
    behaves exactly like carbon: CASSCF ³P with −0.127 Eh correlation, SA
    restores ⟨L²⟩ = 2.000000. Lesson: check C†SC − I of anything fed to an
    optimizer — every ladder rung carries `ortho_err`.
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
    average-of-configuration SCF must select its fractional block by
    CHARACTER (AO projection), not energy order: on Sc the aufbau block
    runs away to {4s, 4p, ~2d} — a converged but p-containing ensemble,
    not 3d1 4s2. The projection weight must be the S-projector
    (Sc)_o† S_oo^{-1} (Sc)_o (a bare quadratic form c_o†S_oo c_o is
    unbounded for diffuse bases). (d) TM ladder fingerprints: the raw
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
    `radial_mask`): rotations allowed only between same-(l, member) slots
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
17. **The ~1e-3 CASSCF-rung commutator class is the INPUT window's
    core-orthogonalization defect** (ticket 14; resolves the relax-rung
    attribution left open in fact 12/15): the raw shell window is exactly
    closed, but every CASSCF must START from its core-orthogonalized form
    (build_full_mo), and that object already carries the defect BEFORE
    any optimization (Sc: 6.8e-4 on the input window == the relax rung's
    converged value; Si ~1e-3; C ~1e-9). Mechanism: projecting window AOs
    against the canonical core span subtracts directions that are not
    rotation-closed (canonical MOs carry m-asymmetric tails); the leak is
    O(⟨window|core⟩ × core nonclosure) — Sc's 4s↔3p-MO overlap ~0.3 makes
    it 1e-3. NOT core polarization (refolding the relax rung's optimized
    core with the raw window commutes at 1e-13) and NOT the optimizer
    (sa_radial, frozen canonical core, shows the same 1e-3). Removing the
    baseline needs a rotation-closed projector — i.e. an l-purified core
    span, which is exactly ticket 13. Read every Si/Sc CASSCF-family
    commutator as baseline + optimizer gauge.
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
    Numbers/curves: README §7.4; anatomy: §3.6.
19. **The one-particle closure meter is the universal certificate**
    (ticket 16 P2; `operators.window_closure`, `avas_window`,
    `span_overlap`): ‖[P,Lz]‖ and ‖[P,L⊥]‖ of the window SPAN projector
    (S-metric) close the fact-18c blind spots — [P,L⊥] sees EVERY
    member asymmetry; both ≈ 0 certifies closure completely. Carried by
    every ladder row (|[P,Lz]|, |[P,Lp]|) and StageFacts payloads.
    Span overlap is the IDENTITY axis, independent of closure (Sc's
    AVAS window: 0.9855 overlap yet O(0.5) closure defects — overlap
    says nothing about closure). AVAS verdict: right size,
    mean-field-inherited span; raw shells strictly dominate. SA
    deviations decompose 100% window / 0% core-fold (witness: refold
    the converged core with the raw shell window; artifacts
    results/*_sa_decomposition.md). Verdict tables: README §7.3.
20. **The closed core span** (ticket 13; `purified_core_orbitals`, the
    `*_pur` rungs — fact 17's fix, verified C/Si/Sc/Ti): per l, take the
    member-resolved radial blocks R of the core MOs; E_l = top-k span of
    the S-metric COLUMN Gram RᵀSR (eigenvalues = S-weight × member
    multiplicity; the tail level is both the selection and the
    certificate); replicate every radial on ALL members; one global
    Löwdin (S = 0 between different l keeps l-purity). TRAPS: the
    S-generalized eigenproblem of RRᵀ returns S⁻¹-warped radials — the
    column Gram is the correct object; canonical core MOs are NOT
    shell-confined (C 2s: 66% of S-weight off its shell) — per-shell
    projection is WRONG. Result: both input windows (build_full_mo and
    radial_basis) close at machine precision; the purified core density
    is exactly m-scalar (avg fold == raw fold); frozen-window pur rungs
    (relax_pur, sa_relax_pur) deliver CASSCF/SA-CASSCF with unbroken
    Lz/L² on every atom incl. Ti (comm ~1e-13, max|dE| ≤ 1e-7). Free SA
    is machine-exact on p² on this substrate; on Sc/Ti the member gauge
    is exactly FLAT, so sa_radial_pur stops at its floor (~1e-6..1e-4
    window tilt) that no convergence tightening removes, and the 4s↔3p
    swap attractor SURVIVES purification (branch on window weights).
    Energy bookkeeping: the pur family carries the window-definition
    offset vs raw-fold CASCI rows (Sc +88.7 mEh — the projected window
    is the honest common definition WITHIN the pur family);
    relax-family max|dE| ~1e-3 vs unpolarized CASCI is fold bookkeeping
    (bound, don't pin). Narrative: README §5.3, §7.3.

21. **Ti: d², the aufbau runaway, and the basin gate** (ticket 12):
    (a) THE REFERENCE IS TWO-LAYERED CHAOS. Default ROHF on Ti/ccpvdz
    (spin=2) is an aufbau runaway (doubly-occupied 3d inside the
    canonical core, +123.07 mEh above the physical 3d2 4s2 ³F);
    energy-ordered get_occ falls back from EVERY start. Fix:
    `run_rohf(occ_closed=, occ_open=)` — character-forced occupations
    via the fact-15c S-projector metric; lands −848.40648875 with
    pure-3d open MOs. Second layer: even with forced occupations the
    SCF is MULTI-BASIN (orbitally degenerate open shell): a metastable
    +11.2 mEh basin exists (same occupations; its polarized core smears
    the build_full_mo input window at 1.6e-2 vs ~1e-8 in the good
    basin; default minao/vsap guesses reach it, atom/huckel reliably do
    not) — add the `occ_e_ref` energy gate (retried over cycled guesses
    until |E−ref| ≤ 1e-4, else RuntimeError; cached chkfiles are
    trusted). The gate also fixes the relax rung's substrate (good-basin
    input defect ~1e-8, so casscf_relax sits at 2.3e-8 — NOT Sc's
    fact-17 class; that class is basin/substrate-specific, not
    atom-intrinsic). (b) Tables (pinned in tests): option A CAS(2,5)/{3d}
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
    README §5.5, §7.2–7.3, §7.4.

## Conventions

- Atomic units, ℏ = 1. Spin-orbitals: alpha = `0..norb-1`, beta = `norb..2norb-1`.
- Determinants: sorted tuples of occupied spin-orbitals.
- Every operator enters the determinant space as (one-body matrix,
  two-body integral tensor) — never creation/annihilation algebra.
- New physics code ships with a test that would catch a sign or phase error
  (spectra, commutators, or independent-construction agreement).

## Tests

`pytest` runs 91 tests in ~13 min (chkfile caching; ROHF computed once per
session; the Ti ladder is the slowest block, ~7 min cold; the Sc ladders
~2.5 min). Test files map
1:1 to src modules; tests are named for the
invariant they pin. The critical invariants: `test_hamiltonian_matches_pyscf_fci`
(Slater–Condon vs PySCF FCI), `test_elementwise_identical_to_slater_condon`
(second quantization ≡ first quantization, elementwise),
`test_full_space_vs_pyscf_casci_core_explicit` (core-explicit CAS(6,5) vs
PySCF ncore=0), `test_full_multiplet_ms0` (textbook ¹S/¹D/³P),
`test_commutators_vanish_in_shell_space` (exact symmetry),
`test_analytic_matches_libcint` (two independent L constructions),
`test_atom_term_tables[*]` (generic pipeline per atom; new atoms = new rows),
`test_ladder.py` (per-rung fingerprints, parameterized C/Si; Sc TM block:
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
(the ticket-09 guard's negative control),
`test_casscf_input_window_orthogonalization_defect` (fact 17),
`test_radial_basis_labels_and_partners` + `test_lm_sectors_partition_fused_shells`
(the labeled substrate), `test_avas_cross_validates_the_shell_window`
(independent ecosystem window selection), `test_window_closure_certificate_and_blind_spots`
and `test_window_closure_sc_shell_and_defective` (the one-particle meter:
exact certificates, the fact-18c blind spots closed by [P,L⊥]),
`test_avas_window_reproduces_the_shell_cascade` (fact 19: AVAS inherits
the mean-field relaxation — C = the full θ* gain near-symmetric, Sc = +80 mEh
and O(0.5) closure defects), `test_stages.py` (the facts-module
contract), `test_scans.py` (frozen-scan anatomy: full radial gain at exact
symmetry, exactly-flat single-member tilts, the Lz blind spots, the Sc
inversion, and the Ti d² verdict — tilts descend, radial rises, the
Lz-blind ±pair pattern — deterministic, no optimizer),
`test_half_integer_spin_labels` (doublet
labeling), `test_drivers.py` (build_full_mo block properties — the
ticket-09 regression guard — shell_weights partition, frac-occ ordering),
`test_analysis.py` (configuration labels), `test_system.py` (validation).

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
