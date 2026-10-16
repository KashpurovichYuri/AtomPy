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
| `src/atomterm/operators.py` | `build_operators`, `full_space_integrals`, `window_closure` | assembles H, Lz, L^2, S^2, N as full `(ndet, ndet)` matrices; core folded into h1 (mcscf convention, spherically averaged); `method="slater_condon"|"secondquant"` switches the many-body evaluation; `window_closure` = one-particle closure certificate ‖[P,Lz]‖/‖[P,L⊥]‖ of any window span (fact 19a) |
| `src/atomterm/drivers.py` | `run_rohf`, `run_frac_occ_scf`, `run_casci`, `constrained_casscf`, `core_virtual_only_mask`, `no_core_window_exchange_mask`, `lm_sectors`, `radial_basis`, `radial_mask`, `run_sa_casscf`, `run_casscf_*`, `select_valence_ao`, `select_shell_ao`, `shell_weights`, `purified_s_core`, `purified_core_orbitals`, `avas_window`, `span_overlap` | all checkpointed (`.chk`); `constrained_casscf(mo_full=, mask=, weights=)` is the one deep CASSCF driver — a mask is any `uniq_var_indices`-shaped callable, and the `run_casscf_*` names are one-line adapters over it; `lm_sectors`/`radial_basis(c_core=)`/`radial_mask` build the (l, member)-labeled substrate for same-sector (radial) constraints (facts 16e/17); `purified_core_orbitals` = the exactly L-closed core span (fact 20: shared radial per closed shell, E ⊗ members closure theorem); active space = raw shell AO columns placed at `[ncore:ncore+ncas]`; core-explicit tests use `purified_s_core`; `avas_window`/`span_overlap` = the ecosystem window probe (fact 19c/d) |
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
18. **The frozen-scan anatomy of the radial direction** (ticket 16,
    `scans.py`, optimizer-free): (a) on p² atoms the ENTIRE CASSCF gain
    is a ONE-PARAMETER m-symmetric radial rotation at EXACT symmetry —
    C: θ* = 0.61 rad over 2p→3p, gain −0.2508 Eh (= the full free/radial
    CASSCF gain); Si: θ* = 0.85 over 3p→4p, −0.1267 Eh; [H, Lz] = 0 at
    every point of the curve. Mechanism visible in the one-/two-body
    split: active-active repulsion falls faster than the (nucleus+core)
    attraction rises. (b) Single-member tilts cost EXACTLY 0 Eh (to
    1e-12): each M_L component of the ³P multiplet avoids one member, so
    the CI finds an untouched component — the term's rotational
    degeneracy IS the flat direction of fact 16; tilting ≥2 members costs
    real energy (C: +1.6 mEh/pair at δ=0.1). This is why free single-root
    CASSCF breaks m-symmetry for free, and why SA (which averages all M_L
    components, every one of which feels a single-member tilt) has the
    symmetric window as its stationary point. (c) METER BLIND SPOTS:
    [H, Lz] cannot see a tilt of the m=0 member (Lz pz = 0) nor EQUAL
    tilts of a full ±m pair (Lz mixes only within ±pairs) — such
    deformations leak only through L± / L², which is trivially zero in
    maximal-Ms single-irrep sectors. Sector commutators are necessary,
    not sufficient, closure witnesses. (d) Sc INVERTS the structure: the
    m-symmetric 3d→4d rotation gains NOTHING (the curve rises — radial
    diffusion of a one-electron d shell only loses nuclear attraction;
    4s→5s likewise harmful), while MEMBER-ASYMMETRIC tilts DESCEND
    (−0.3 mEh for one tilted member at δ=0.1, −0.9 mEh for four, with
    comm ~ 0.2–0.45): a genuine Jahn–Teller-like instability along the
    degenerate ²D manifold — the drive behind Sc's chaotically broken
    single-root radial solutions; SA over the quintet averages the
    instability away (sa_radial stays symmetric). Scan numbers are
    DETERMINISTIC (no optimizer, no basin chaos).
19. **One-particle closure meter, SA-deviation decomposition, AVAS
    verdict** (ticket 16 Phase 2; `window_closure`, `avas_window`,
    `span_overlap`, `scripts/sa_decomposition.py`): (a) ‖[P,Lz]‖ and
    ‖[P,L⊥]‖ = ‖[P,Lx]‖²+‖[P,Ly]‖‖¹ᐟ² of the window SPAN projector
    (S-metric) are the complete closure certificate with NO blind spots:
    [P,Lz] reproduces the many-body blind spots exactly (m=0-member,
    equal-±pair tilts invisible), [P,L⊥] sees EVERY member asymmetry.
    Carried by every ladder row (|[P,Lz]|, |[P,Lp]| columns) and by
    StageFacts payloads. It exposes what sector commutators cannot:
    C's converged casscf_radial window is member-tilted at [P,Lp] ≈ 1.1
    (the flat gauge of fact 18b the single-root optimizer wandered into)
    while its many-body commutators sit at 5e-7; sa/sa_radial anchor the
    symmetric gauge ([P,Lp] ~ 1e-4). (b) Decomposition verdicts
    (artifacts results/*_sa_decomposition.md): the core-fold witness
    (refold the rung's converged core with the RAW shell window) is
    ≤ 1e-13 for every SA rung on all three atoms — SA deviations are
    100% window (identity + member structure), 0% core/fold. Si refines
    fact 16: FREE sa/sa_pyscf heal below the 3.4e-3 input baseline
    (→ 4e-6/7e-8) — only the constrained sa_radial is pinned at it.
    Sc (this run): sa = the 4s↔3p swap basin (sv 0.00, comm 0.83/3.8);
    sa_radial = exactly the input baseline; sa_pyscf = collapse;
    canonical-START baselines are huge (Si 0.24/0.32, Sc 2.8/3.7 — the
    SCF window block is member-polarized; that is what sa_pyscf eats).
    (c) Span overlap ≠ closure — independent axes: Sc's AVAS window
    overlaps the shell span at 0.9855 yet carries closure defects
    0.5/0.4 (≈ 0.17-rad member tilt) and CASCI commutators 0.59/2.6;
    overlap-based "window quality" says nothing about rotational
    closure. (d) AVAS verdict (verified; refutes the same-span guess):
    right SIZE, but the span is the mean-field-inherited relaxation —
    on C/Si exactly the fact-18 m-symmetric optimum (−250.8/−126.7 mEh
    = the full θ* gains, near m-symmetric), on Sc the diffuse canonical
    virtual d-set (+80.8 mEh WORSE than raw shells, <L²> 5.90
    masquerading as ²D). Raw shells strictly dominate for our purpose;
    sv(AVAS) ≈ sv(CASSCF family) on C (0.917) and Si (0.883) — three
    independent constructions agree on the relaxed span. (e) Sc
    tilt-follow (deterministic, <L²>-tracked): the fact-18d descent
    deepens monotonically to −2.7 mEh (subset {0,1,2}) / −3.6 mEh
    ({0,1,3,4}) at δ=0.35 with ground-state L drifting 2.00 → 1.87/1.84
    — the anatomy of the optimizer's ²D→²P relabeling; {0,1,3,4} is
    Lz-BLIND ([H,Lz] ~ 1e-14) while [H,L²] ≈ 1–2 fires (this Ms=1/2
    multi-term sector catches it; maximal-Ms single-irrep sectors would
    catch NOTHING). m-symmetric reference rises steeply (+51.6 mEh at
    δ=0.30) with <L²> = 6.000 exact. The m-scalar-constraint rung is
    NOT built (decision): the member gauge is exactly flat on p²
    (fact 18b) and sa_radial already realizes the m-symmetric optimum
    reproducibly — there is nothing to gain.
20. **The closed core span, realized** (ticket 13; `purified_core_orbitals`,
    the `*_pur` ladder rungs; fact 17's fix, verified on C/Si/Sc):
    (a) PREMISE CORRECTION: canonical core MOs are NOT confined to their
    own (n,l) AO shells (C's 2s MO has 66% of its S-weight on other
    radials — same-l radial mixing is physical), so "project each core
    MO onto its own shell" is WRONG (gate (a) failed at density 0.77).
    The correct closure theorem is sharper: for one-center bases the
    l-chain rows decompose into MEMBER sectors with identical radial
    chains (S = I_m ⊗ S_rad exactly), and L acts ONLY on the member
    index — so E ⊗ members is L-invariant for ANY radial subspace E.
    Shared radial comes free; full AO shell subspaces are unnecessary.
    Construction: per l, the member-resolved radial blocks R = [r_{i,m}]
    of the core MOs; E_l = top-k span of the S-metric column Gram RᵀSR
    (eigenvalues count S-weight × m-multiplicity — ~2l+1 per occupied
    radial, tail-level beyond: the selection AND the certificate);
    replicate every radial on ALL members; one global Löwdin (S = 0
    between different l keeps l-purity). Beware the metric trap: the
    S-generalized eigenproblem of RRᵀ returns S⁻¹-warped radials — the
    COLUMN Gram is the correct object. (b) Verified gates: input windows
    (build_full_mo AND radial_basis) close at machine precision —
    Si 3.4e-3/4.5e-3 → 1.1e-15/1.8e-15, Sc 5.7e-3/2.1e-2 →
    6.6e-15/8.7e-15; the purified core's density is m-scalar EXACTLY
    (equal member-diagonal blocks, zero cross-member coherences,
    anisotropy 1e-16, raw fold == averaged fold → casci_avg_pur's
    max|dE| = 0.0) — it is the orbital realization of the scalar core
    the fold already uses. Reference shift = the discarded core
    polarization only (spherically averaged density change C 4.2e-6,
    Si 3.6e-6, Sc 1.05e-4; energy C +3.0 µEh, Si +128 µEh, Sc
    +674 µEh — C's equals the fact-12 avg-vs-raw shift exactly).
    (c) THE DELIVERABLE: Sc CASSCF/SA-CASSCF with unbroken Lz and L².
    Frozen-window rungs close at machine precision on every atom:
    Sc relax_pur 1.3e-13, sa_relax_pur (SA + core-relax mask, window
    frozen) 6.7e-14/3.6e-13 with the full 90-state term table and
    max|dE| 5e-8. FREE SA closes exactly where the optimizer has a
    restoring force: C sa_pur/sa_radial_pur = 0.0 exactly, Si 3.3e-16
    (vs 4e-6 healing on the canonical substrate) with the full θ*
    gains. On Sc the member gauge is exactly FLAT (fact 18d), so
    sa_radial_pur stops at a ~1e-6 window tilt (comm 5.7e-7/1.8e-5)
    that NO convergence tightening removes (conv_tol 1e-7→1e-12: [P,Lp]
    3.8e-6→1.0e-6, energy identical) — an optimizer floor with no
    restoring force, not a baseline. The free-sa 4s↔3p swap attractor
    SURVIVES the purified substrate (window p = 1/6 exactly again) —
    it is an optimizer basin, closing the fact-17 attribution question.
    Single-root casscf_radial_pur keeps the flat member-gauge chaos of
    fact 16 (C [P,Lp] 1.1 with [H,Lz] 0.0 — the blind-spot demo again).
    (d) Phase 4 (post-hoc window repair) NOT needed — free SA is
    machine-exact on p² atoms on this substrate. Bookkeeping note:
    pur-rung energies carry the same fact-14 window-definition offset
    vs raw-fold CASCI rows (Sc +88.7 mEh; the projected window is the
    honest common definition *within* the pur family); relax-family
    energy SIGNS vs the unpolarized CASCI are fold bookkeeping (our
    averaged fold discards the polarization the optimizer bought:
    max|dE| ~ 1e-3 is the detector), not physics — bound, don't pin.

## Conventions

- Atomic units, ℏ = 1. Spin-orbitals: alpha = `0..norb-1`, beta = `norb..2norb-1`.
- Determinants: sorted tuples of occupied spin-orbitals.
- Every operator enters the determinant space as (one-body matrix,
  two-body integral tensor) — never creation/annihilation algebra.
- New physics code ships with a test that would catch a sign or phase error
  (spectra, commutators, or independent-construction agreement).

## Tests

`pytest` runs 83 tests in ~4.5 min (chkfile caching; ROHF computed once per
session; the Sc ladders are the slowest block, ~2.5 min cold). Test files map
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
floor, the surviving swap attractor — branched on window weights),
`test_purified_core_orbitals_exact_closed_span` +
`test_purified_substrate_closes_the_input_windows` (fact 20: the closure
certificate on the core span and both input windows, m-scalar density,
contamination/density-shift scales — the gate-(a)/(b) pins),
`test_radial_rung_full_gain_without_l_leakage` (fact 16a: radial == free
gain, zero l-leakage), `test_divergence_guard_fires_on_broken_orthonormality`
(the ticket-09 guard's negative control),
`test_casscf_input_window_orthogonalization_defect` (fact 17),
`test_radial_basis_labels_and_partners` + `test_lm_sectors_partition_fused_shells`
(the labeled substrate), `test_avas_cross_validates_the_shell_window`
(independent ecosystem window selection), `test_window_closure_certificate_and_blind_spots`
and `test_window_closure_sc_shell_and_defective` (the one-particle meter:
exact certificates, the fact-18c blind spots closed by [P,L⊥]),
`test_avas_window_reproduces_the_shell_cascade` (fact 19d: AVAS inherits
the mean-field relaxation — C = the full θ* gain near-symmetric, Sc = +80 mEh
and O(0.5) closure defects), `test_stages.py` (the facts-module
contract), `test_scans.py` (frozen-scan anatomy: full radial gain at exact
symmetry, exactly-flat single-member tilts, the Lz blind spots, the Sc
inversion — deterministic, no optimizer),
`test_half_integer_spin_labels` (doublet
labeling), `test_drivers.py` (build_full_mo block properties — the
ticket-09 regression guard — shell_weights partition, frac-occ ordering),
`test_analysis.py` (configuration labels), `test_system.py` (validation).

## Planning

`wayfinder/map.md` holds the effort map (destination, decisions, fog) and
lists the live tickets. Update it when a decision is made; keep tickets
one-session sized. The live wave (12 NEXT: Ti — options A then B in one
session; tickets 13/16 RESOLVED, facts 18/19/20; the ticket carries a
pre-flight-verified ROHF occupation fix — default ROHF on Ti is an aufbau
runaway; PHYSICS.md narrative deferred by the user)
is fully specified
in `wayfinder/tickets/` — each ticket
is self-contained (no conversation context needed). Review protocol:
`/code-review` (two axes, parallel sub-agents) after each ticket wave, then
`/improve-codebase-architecture` survey; see ticket 13.
