# U / Stuttgart-RLC-ECP78 + cut-Dyall — CASCI campaign report (u_14e)

Reference: ROHF E = −50.38185233 Eh, 5f³6d¹7s² (Newton + forced whole-l
spans "l0","l1" closed / "l2","l3" open, occ_e_ref-gated; 4/4 guesses
identical). ncore = 5 ({6s,6p×3,7s} = the occupied closed block —
canonical order is MEANINGLESS on this basis: the closed block sits at
canonical indices 22/61/88/120/143). 7s² frozen in the core: 7s² is
¹S-closed, so the ground term of 5f³6d¹7s² is the term of 5f³6d¹.

## Window selection (measured, not labeled)

PySCF numbers uncontracted shells from the TIGHT end starting at the
valence name ("5f" = the exp-36.8 primitive — a trap). The 5f open-MO
weight lives on mid-tight primitives: kin 1.325 (0.145), 3.032 (0.351),
6.478 (0.503), 13.16 (0.203); the diffuse end is empty. Windows are
selected by these weights via select_diffuse-style kin cuts; primitive
m-set windows close at 2e-13. Multi-primitive same-l windows are NOT
S-orthonormal: analysis must use the build_full_mo-orthonormalized
window (fact 4; raw columns break product operators — L² eigenvalues
172.98, non-quantized, measured).

## Substrate: pur is mandatory on U

Raw canonical core: input-window defect 0.27, [H,L²] = 2.2 (the fact-17
mechanism at U scale — the 6p core is directly adjacent to 5f).
purified_core_orbitals (ncore=5, contamination 1.3e-2): window closes
at 5.1e-13, commutators 1.3e-12, exact multiplets. All tables below on
the pur substrate.

## The window-convergence series (sector (4,0); block grounds in Eh)

| window | bottom configuration | f³d¹ relative gap |
|---|---|---|
| 1f + 1d (12 AOs, 495 dets) | f¹d³ | +306 mEh |
| 2f + 1d (19 AOs, 3,876 dets) | f²d² (even ⁵L ×17) | +52 mEh |
| 3f + 1d (26 AOs, 14,950 dets) | **⁵L° odd, L=8, S=2, ×17 exactly degenerate** | — |

Block grounds, 2f window: f⁴ −3.1776, f¹d³ −3.3346, f²d² −3.5034,
f³d¹ −3.4513. 3f window ground: −3.745754, ⟨L²⟩ = 72.00 exactly,
odd, ⟨N_f⟩ = 2.33 / ⟨N_d⟩ = 1.67 — an f³d¹-dominated mix with f¹d³
(pair transfer f³d¹↔f¹d³ conserves parity; the even f²d² manifold is
EXACTLY excluded by parity). Second multiplet: ⁵K° (L² = 56) at
−3.740461.

**Verdict (part 1)**: the physical ground term of the 14-electron U valence
space in this basis/ECP is **⁵L° (5f³6d¹-parented, L=8, S=2, odd)** —
the LS parent of the experimental ⁵L°₆ ground level. The configuration
ordering is window-sensitive (1f: f¹d³; 2f: f²d²; 3f: f³d¹): radial
flexibility of the f shell decides it; window convergence is an
experiment axis for actinides. A 4f window point (adds kin 13.16,
40,920 dets) is queued as the stability confirmation.

*[Part-2 amendment]*: confirmed on the f side (4f stable), transiently
retracted on the d side (2d flip at 3f), and CONFIRMED at the symmetric
corner: 4f+2d restores ⁵L° f³d¹ by > 37 mEh and crosses below the ROHF
(−50.3920 < −50.38185). The ⟨N_f⟩ = 2.33 mix reading is also corrected
(configuration-pure f³d¹; see part 2).

## Engine notes (this campaign's enabling work)

- slater_condon_matrix VECTORIZED (elementwise-identical to the kept
  loop reference, <1e-13 reassociation noise; 18,411 dets: hours →
  21.7 s). Real operators now build real matrices (2× memory traffic,
  ~2.5× faster eigh). Expectation sums over eigenvector COLUMNS are
  truncated to the roots actually read (full-spectrum einsums are the
  O(n³) trap).
- BLAS thread policy: 128 threads by default (user rule, 2026-10-21);
  255-thread auto-pools LIVELOCK on our matrix sizes (the "silent
  deaths": 36 GB RSS, no crash, no output — spin-wait contention).
- Per-block PySCF cross-checks: our fold vs PySCF g at 6.7e-16, h1 raw
  5.3e-15 / averaged 5.5e-3 (the standard class); max|dE| per sector
  2.5–5e-3 (tier 1a).

## Part 2 (2026-10-21): exponent bands, standard tables, the two-axis race

### Window protocol, made mechanical

`select_diffuse` (kinetic-diagonal proxy) is RETIRED: `primitive_chain`
+ `select_primitives` walk the basis `lm_sectors`-style and select by
EXPONENT — bands are the contiguous sets around the measured open-MO
Mulliken peaks (printed by the driver at start; never hand-picked):

- f chain: 1.44 (1.12) > 0.674 (0.91) > **2.92 (0.53)** > 0.294 (0.41) > 0.113 (0.06)
- d chain: 0.172 (0.41) > **0.372 (0.30)** > 0.076 (0.22)

Bands: 1f=(1.0,2.0), 2f=(0.5,2.0), 3f=(0.2,2.0), 4f=(0.2,4.0); the 2d
variant widens d (0.1,0.5). Contracted shells raise a clear error
(they are pre-converged radials — use label/MO selection there).

### The standard report (per-configuration term tables)

`scripts/run_u.py` emits term × parity × ⟨N_f⟩/⟨N_d⟩ × dominant
configuration × E per root, block grounds + per-block max|dE| vs PySCF
(`population_operator` = the exact S-metric row projector; integer
per-determinant N_f gives the exact config weights). Full tables:
`results/u_14e/u_run_part2.log` (10 sectors, 2h17m).

**A part-1 correction**: the ⁵L° ground is configuration-PURE f³d¹
(⟨N_f⟩ = 3.000, w = 1.00). Part 1's ⟨N_f⟩ = 2.33 "f³d¹/f¹d³ mix" was an
artifact — L = 8 is unreachable from f¹d³ (max L = 5): angular-momentum
coupling forbids the mixing the old number implied.

### The window lattice (Ms=2 block grounds, Eh, total scale)

| window | f³d¹ odd (⁵L° ×17) | f²d² even (⁵L ×17) | other blocks in 48 roots | ground |
|---|---|---|---|---|
| 1f+1d | — | — | f¹d³ −49.4587 (⁵I° ×13) | f¹d³ |
| 2f+1d | — | −49.7506 | — | **f²d²** |
| 3f+1d | −49.9928 | (above roots) | — | **f³d¹ ⁵L°** |
| 4f+1d | −50.3320 | (above roots) | f⁴ −50.3070 (+25 mEh) | **f³d¹ ⁵L°** |
| 2f+2d | — | −49.8943 | f¹d³ −49.8763 | f²d² |
| 3f+2d | −50.0617 | **−50.0809** | — | **f²d² (⁵L° +19.1 mEh)** |
| 4f+2d | **−50.3920** | > −50.3547 | — | **f³d¹ ⁵L° ×17, margin > 37 mEh** |

Verdicts: (a) the f-side tight neighbor (4f) does NOT reorder — ⁵L°
stable, +339 mEh absolute correlation, f⁴ block surfaces at +25 mEh;
(b) the d-side neighbor (2d) alone FLIPS the 3f ordering to even ⁵L by
19 mEh — but (c) at the SYMMETRIC corner (both axes grown once) the
flip un-flips: ⁵L° f³d¹ wins by > 37 mEh (f²d²'s ground above the 48th
root at −50.3547) AND crosses below the ROHF (−50.3920 < −50.38185,
the window finally represents the mean field; measured-gradient
extrapolation predicted −50.40). The 3f2d flip was a transient of the
asymmetric lattice point. **The ground term of the U 14-electron
valence space is ⁵L° (5f³6d¹-parented, L=8, S=2, odd), stable at every
saturated corner.** Last open neighbor: the d-chain's 3rd primitive
(exp 0.076, weight 0.216 — f side saturated at 98.2% weight; a
flip-back needs f²d² to gain > 37 mEh MORE than f³d¹ from it).

### Engine notes (part 2)

- **PySCF's Davidson straggles inside degenerate multiplets**: up to
  38 mEh within exact 13–17-fold clusters at nroots=48 (2f2d, measured)
  — single-root tolerances pass while cluster interiors drift. Fix:
  `run_casci(fci_opts={"pspace_size": ndet})` forces their exact branch;
  both-sides-exact max|dE| = 7–8e-13/1e-15 (our Lanczos vs their dense).
  Above EXACT_MAX=12k dets the PySCF column is Davidson-class (labeled).
- **Dense LAPACK eig is the wall**: memory-bound ~40–80 GFLOPs
  effective regardless of threads (370 s at 14,950 part 1; >45 min
  incomplete at 18,411 — killed). A subset `dsyevr` still pays the full
  tridiagonalization. The truncated-roots path now uses ARPACK Lanczos
  (`low_roots`, k=nroots+16 for the degenerate clusters, fixed v0):
  40,920 dets ≈ 35 min/sector all-inclusive; validated against the
  both-sides-exact 2f2d spectrum (8.6e-13).
- Commutators stay machine class everywhere (4e-11 at 40,920 dets);
  pur substrate throughout (raw core would read 2.2).

## Pending (next waves)

4f+2d resolution of the two-axis race (running; 73,815 dets); then the
SA wave (ticket 23: 17-root SA over the ⁵L° M_L set — does orbital
relaxation replace window growth on either axis); window methodologies
(ticket 22: replicated SCF radials / NO selection as compression of
this lattice); NIST comparison (deferred).
