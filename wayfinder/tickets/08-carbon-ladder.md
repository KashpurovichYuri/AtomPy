# Ticket 08: Carbon ladder — verify the commutator/symmetry story end to end

Status: resolved

## Question

Every symmetry claim of the project (exact terms in a shell CASCI, spherical
averaging, SA restoration) has so far been demonstrated only in its final
form. Can we verify the *whole ladder* of orbital treatments on carbon —
including deliberately un-averaged and freely optimized rungs — and compare
all rungs against each other?

## Resolution (to implement)

New `ladder.py`: one analysis pass per rung, then an all-vs-all comparison
matrix (energies, term tables, <L^2>, both commutator norms, window
composition, max|dE| vs PySCF). Carbon only.

Rungs (7):

1. **CASCI, raw core** — `spherical_core=False` through
   `build_operators`. Expect: 3P split ~4e-5 Eh in the {M_L=+-1}/{0} pattern
   (quadrupole k=2, m=0 selection rule: dM_L = 0, +-2), commutators ~1e-5,
   term letters still right. Report ||rho - rho_sph||.
2. **CASCI, averaged core** (current default) — exact zeros (gates).
3. **CASCI, fractional-occupation core** — average-of-configuration SCF
   (open shell occupied fractionally: 2/3 per 2p AO), fold with
   `spherical_core=False`: Fock is exactly spherical, so commutators ~ 0 by
   construction; energies must match rung 2 to ~1e-10 (O(eps^2) test of the
   "averaging = ensemble projection" theory).
4. **CASCI + core relaxation** — window span frozen, only inactive<->virtual
   rotations allowed (gradient-block zeroing in PySCF CASSCF if a clean hook
   exists; otherwise defer as follow-up). Expect E <= rung 2, exact
   analysis commutators, max|dE| revealing the relaxed core's anisotropy.
5. **CASSCF, shell start** (current) — <L^2> ~ 1.67, comm ~ 2e-2 (fact 8).
6. **CASSCF, PySCF-default start** — same class; start-dependence measured.
7. **SA(3), uniform** — <L^2> = 2.000000, comm ~1e-11.

Machinery:

- folding fix: optimized stages fold with the stage's OWN inactive columns
  (`cass.mo_coeff[:, :ncore]`), not the ROHF ones; CASCI stages unchanged
  (identical columns by construction).
- window composition: `shell_weights(mol, c, s)` — S-metric projection
  norms `w_i(n,l) = c[:,i]^dag S_shell c[:,i]` per column (rows sum to 1);
  stored per stage, printed in reports.
- docstrings: commutator semantics ([H,Lz] vs [H,L^2] vs full group);
  `H_eff = H0 + lambda W` formalization; parity -> quadrupole (m=0) ->
  {+-1}/{0} split chain; spherical averaging = ensemble average over the
  degenerate ROHF manifold, agrees with fractional-occupation core at
  O(eps^2); core-relaxation peculiarities (relaxed core anisotropy measured
  via max|dE|, never hidden).

Tests (`tests/test_ladder.py`): pattern assertions — raw-core split > 0 but
< 1e-3 with bounded commutators; averaged/frac-occ exactness gates; frac-occ
vs averaged energies ~1e-10; CASSCF <L^2> ~ 1.67 +- tolerance; SA <L^2> =
2.000000. All 34 existing tests green.

## Acceptance

- Ladder runs end-to-end on carbon; comparison matrix in report + JSON.
- Pattern assertions green; existing suite green.
- AGENTS.md gains the raw-core fingerprint + fractional-occupation
  equivalence facts. Si/Sc/Ti tickets then reuse the machinery unchanged.

## Resolution (closed)

All 7 rungs implemented and compared (40 tests green). Measured matrix
(carbon, E_ground / |[H,Lz]| / max|dE| vs PySCF / window p-weight):

- casci_raw: −37.43161, 5.6e-6, 0.0, p=1.0 — the internal field is real
- casci_avg: −37.43158, 0.0, 2.9e-5, p=1.0 — exact invariance
- casci_frac: −37.42500, 1.2e-20, 0.0, p=1.0 — invariant by construction
- casscf_relax: −37.44037, 0.0, 2.4e-5, p=1.0 — −8.8 mEh core relaxation,
  zero escape (uniq_var_indices mask hook worked as designed)
- casscf_shell/pyscf: −37.68236 (same basin, 2e-9 apart), comm 1e-2..1e-1
  (chaotic gauge), p>0.95
- sa3: −37.68239, 3.3e-12, <L^2> = 2.000000

Findings that CORRECTED predictions (the point of the ladder):

1. **[H, L²] = 0 identically in the 3-det Ms=1 sector**: it is a complete
   L=1 irrep, L² = 2·Id — the raw-core fingerprints there are [H, Lz] and
   the SPLITTING, which shows the predicted {M_L=0} vs {±1} pattern
   (4.4e-5 Eh; roots [−37.4316107, −37.4315670, −37.4315670]).
2. **Fractional-occupation core ≠ averaged core at 1e-10** (ticket
   prediction wrong): they differ by 6.6 mEh. Average-of-configuration is
   a DIFFERENT reference (open-shell exchange with smeared vs localized
   valence), not the ε→0 limit of core polarization. Both invariant;
   difference physical. AGENTS.md fact 12.
3. **(n,l)-grouped window weights are ill-defined** (2p/3p AO subspaces
   overlap — weights double-count). Replaced by l-grouped weights, which
   partition exactly. Consequence: the free-CASSCF window is l-PURE —
   fact 8 refined: the breaking is m-incompleteness within l=1, not s/d
   leakage. l-purity necessary, not sufficient; commutators remain the
   closure authority.
4. **PySCF SA kernels return only the ensemble energy** (scalar e_cas) —
   run_sa_casscf now extracts per-root energies via a CASCI pass in the
   converged window (the pipeline's old "roots" print was a single number).
5. Folding fix applied everywhere: optimized stages fold with their OWN
   inactive columns; for casci_raw/frac the max|dE| = 0.0 confirms our
   fold == PySCF's fold exactly, and casci_avg's 2.9e-5 isolates the pure
   averaging effect.

Tests: 6 ladder pattern tests + all previous 34 green. Si/Sc/Ti reuse
`shell_weights`, the fold fix, and the rung machinery unchanged.
