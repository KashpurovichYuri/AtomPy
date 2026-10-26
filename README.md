# AtomTerm

Atomic term-symbol analysis on top of [PySCF](https://pyscf.org): the pipeline
runs ROHF → CASCI/CASSCF/SA-CASSCF, then classifies every state by
**exact** angular-momentum operators L̂z, L̂², Ŝz, Ŝ² built as plain matrices in
the Slater-determinant basis, and certifies rotational-symmetry preservation
by ‖[Ĥ, L̂z]‖_F and ‖[Ĥ, L̂²]‖_F (Ŝz certifies nothing — spin-free Ĥ commutes
with it identically; it is the exact M_S label of a sector). Atoms covered:
**C** (2p²), **Si** (3p²), **Sc** (3d¹4s²), **Ti** (3d²4s²), **U** (f³d¹;
ECP + cut-Dyall).

## Quick start

```bash
conda activate pyscf          # PySCF 2.13 environment (any platform; in a bare
                              # shell: source <conda-root>/etc/profile.d/conda.sh first)
pytest                        # 128 tests, ~13 min (chkfile caching); -m fast = sub-second tier
python scripts/run_c.py       # writes elements/C/ (report.md/.json + chkfiles)
```

Run from the repo root — `pyproject.toml` preconfigures `pythonpath = src`
and the test paths; no install step.

## What you get

The results, at a glance — committed in `elements/<Atom>/` (Results ·
Files maps every file to its regeneration command).

**Term tables** — every state classified, per Ms sector. Titanium
(CAS(2,5) over {3d}, Ms = +1 sector — the d² triplet manifold; rows
compressed by multiplicity):

```
sector: Ms=+1 (na=2, nb=0)   (10 determinants)

 #           E (Eh)      <L^2>  L      <S^2>    S  term
 0      -1.46866251  12.000000  3   2.000000    1    3F   (x7)
 7      -1.41396932   2.000000  1   2.000000    1    3P   (x3)

Lz eigenvalues (M_L content): [-3.0, -2.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 2.0, 3.0]
Sz eigenvalues (M_S content): [+1] x10   (exactly M_S·1 in a fixed sector)
||[H, Lz]||_F  = 1.386e-14   (relative 3.018e-15)
||[H, L^2]||_F = 5.819e-14   (relative 1.267e-14)
```

The textbook d² multiplet comes out exactly (Hund: ³F ground, ³P
excited) — commutators at round-off, because the {3d} window is exactly
rotation-closed. The Ŝz line documents the sector's M_S; it certifies
nothing, since [Ĥ, Ŝz] ≡ 0 identically (Motivation). Details: Results
(Term table, Multiplet energies); regeneration: Results (Files).

**Treatment rungs** — every orbital treatment priced against the ROHF
anchor (a rung = one treatment row of this comparison), the free variant
beside the symmetric one, each classified by its diagnostics: ΔE vs the
ROHF anchor, ⟨L²⟩, commutators. Titanium's headline rows (ROHF reference
−848.406489 Eh):

| treatment                        | E, Eh       | ΔE (mEh) | ⟨L²⟩    | ‖[H,Lz]‖ | verdict         |
| ---                              | ---         | ---      | ---     | ---      | ---             |
| CASCI, averaged fold             | −848.3959   | +10.5    | 12.0    | 4.2e-09  | exact           |
| CASSCF, free                     | −848.4377   | −31.3    | 0.0503  | 4.3e-01  | escape — broken |
| CASSCF, window frozen (core-relax) | −848.3961 | +10.4    | 12.0    | 4.2e-09  | baseline        |
| SA-CASSCF, free                  | −849.2676   | −861.1   | 12.003  | 1.1e+00  | swap — broken   |
| SA-CASSCF, frozen, purified core | −848.3948   | +11.7    | 12.0    | 1.0e-13  | exact           |

Verdicts: **exact** = commutators and span certificates at the 1e-13
class; **escape / swap** = window-identity basins (Motivation);
**baseline** = a named non-closure mechanism at its signature scale.
Methods compared (Results) is the reading guide; the full per-atom
tables (span certificates, cross-check column, configurations) live in
`elements/<Atom>/TABLES.md`.

---

## Motivation

The exact Hamiltonian of an atom is SO(3) × SU(2) invariant, and by
Schur's lemma every term (2S+1)L is exactly (2L+1)(2S+1)-fold
degenerate, sharp in L and S.

**The restriction.** A CAS run takes its core from a prior SCF
reference, confines the configurations to a small active set — the
**window** — and represents the omitted core electrons by a folded mean
field. The effective Hamiltonian is

  Ĥ = Σ_pq h_pq a†_p a_q + ½ Σ_pqrs (pq|rs) a†_p a†_q a_s a_r.

Splitting the orbitals into core (i — doubly occupied, frozen), active
(t — the CAS window) and virtual (a — empty), the partition gives

  Ĥ_CAS = E_c + Σ_tu f_tu Ê_tu + ½ Σ_tuvw (tu|vw) ê_tuvw,

- **E_c** — the frozen core's own energy (constant).
- **Ê_tu = Σ_σ a†_tσ a_uσ** — spin-summed excitation ("move an active
  electron from u to t"); ê_tuvw likewise for pair excitations —
  together they generate every configuration inside the active space.
- **f_tu = h_tu + Σ_{i∈core} [2(tu|ii) − (ti|iu)]** — the folded
  one-electron part: the bare h (nucleus + kinetic) plus the averaged
  effect of the core electrons on an active electron (the local Coulomb
  potential of the core density, plus nonlocal exchange). This is the
  mean-field approximation: core–active interaction kept static,
  core–active correlation dropped.
- **(tu|vw)** — four-center two-electron integrals among active
  orbitals.

CASCI is the exact diagonalization of this effective Hamiltonian; its
only approximation is the frozen mean field of the inactive space. The
partition is **exact iff the window is orthogonal to the core
orbitals**: a window vector leaning into the core, w = w⊥ + w∥ with
w∥ ∈ span(core), admits configurations with extra density in *occupied*
core orbitals that the antisymmetrized problem forbids — the CI
silently harvests deep-core stabilization the Pauli principle would
have priced. Every window is therefore built core-orthogonal.

**Injection point 1 — the fold.** The active–active part (tu|vw) is a
rotational scalar; **the folded part f carries whatever symmetry the
core density ρ_c has**. ROHF's ρ_c is axial → f is axial → Ĥ_CAS is
axial: the first entry of symmetry breaking. The same is seen through
orbital rotations: the optimized core orbital sits in its polarized
shape φ_c = 1s + λ·d₀ (λ = the polarization amplitude — the core, not
the window, is the problematic object here); rotating it by 90° about x
remixes Y₂₀ into dz² and dx²−y² combinations the frozen space
{φ_c, 2s} does not contain. Since the frozen core sits inside Ĥ_CAS
through f_tu, the leak is inherited:

  U(R)⁻¹ Ĥ_CAS U(R) = Ĥ_CAS[core on a new axis] ≠ Ĥ_CAS unless R is
  about z.

The symmetry group of Ĥ_CAS is the stabilizer of the frozen core
density — rotations about the torus axis, when such an axis exists.
The cure is a symmetric ρ_c; three core treatments are compared in the
project:

- **raw** — ρ_c as the SCF left it: ε unscreened, the detector.
- **averaged** — ρ_c spherically averaged (its l = 0 part kept) before
  the fold: one post-processing step, the ROHF total unchanged to
  sub-mEh, ε removed; the per-atom ΔE anchors: Results, `elements/`.
- **fractional** — a separate SCF with the open shell smeared uniformly
  over the shell (n_open/(2l+1) per orbital): the anisotropy never
  develops; a different functional, not a variational neighbor of ROHF
  (its ΔE: Results, `elements/`).

**What survives: the stabilizer.** For every generator of
Stab(ρ_c) = {R ∈ SO(3) : R ρ_c R† = ρ_c}, [Ĥ_CAS, generator] = 0 — up
to the numerical axiality of the converged ρ_c; [Ĥ_CAS, Ŝz] = [Ĥ_CAS,
Ŝ²] = 0 identically (the field is spin-free); ±M_L partners stay glued
(real Hamiltonian). SO(3) descends to the stabilizer, and the term
branches into its components — never "splits arbitrarily".

**The worked instance — carbon 2p²** (an illustration, not the general
theory). ROHF places the two ³P electrons in {px, py}: the open-shell
density is a torus, axial about z. The ³P multiplet branches by |M_L|
into Π (M_L = ±1) and Σ (M_L = 0) — exactly the σ/π branching of a
linear molecule or a ligand field, which makes the axial-field response
testable against known crystal-field splittings (the verification
program, Results). Other references have other, generally smaller
stabilizers — hence measured raw-fold commutators from the 1e-5 class
(C, an almost perfect axis) to percent class (Ti).

**Injection point 2 — the window.** With P the projector onto the
active configuration space, the exact commutator [Ĥ, L̂] = 0 gives

  [PĤP, PL̂P] = PL̂(1−P)ĤP − PĤ(1−P)L̂P:

if the space and its complement are rotation-closed ([P, L̂] = 0), the
projected generators still close the algebra and Schur survives. This
is what the code checks: ‖[P, L̂z]‖ and ‖[P, L̂⊥]‖ (`window_closure`).
Carbon shows why the window must span the whole valence shell
{px, py, pz} — rotating py about x gives cos θ·py + sin θ·pz, which
stays inside the p span; a truncated set leaks. The active set — the
**window** — is therefore always chosen shell-complete and built
core-orthogonal: closure is imposed by design, never converged to.
Windows are used in three kinds throughout, compared at matched
determinant cost (Results): **m-complete AO shells and bands** —
rotation-closed by design; **raw-ROHF MO selections** — the open
reference orbitals plus valence virtuals, never rotation-closed (the
SCF's open orbitals are m-delocalized); and **natural-orbital (NO)
occupation groups** — NOs of the term-averaged density, which by
Schur's lemma arrive as complete m-sets: rotation-closed by
construction.

**The chain, measured.** ε — the core anisotropy — against the
splitting it produces (raw fold):

| Atom | ε (relative) | splitting (Eh) | splitting/ε (Eh) | ⟨L̂²⟩ of the ground term |
|---|---|---|---|---|
| C | 2.405e-3 | 4.373e-5 (³P) | 0.02 | 2.0000 |
| Si | 2.109e-3 | 2.431e-4 (³P) | 0.12 | 2.0000 |
| Ti | 7.754e-3 | 4.144e-3 (³F septet) | 0.53 | 12.0000 |
| Sc | 1.259e-2 | 1.857e-3 (²D quintet) | 0.15 | 6.0001 |
| U | — (not measured; U runs on the pur substrate) | | | |

Two readings. The splitting is finite while ⟨L̂²⟩ stays exact: the fold
is spin-free, so terms of different S never mix (¹D admixes ¹S; ³P can
never admix ¹D) — the degeneracy dies (the branching), the high-spin
label is SU(2)-protected. And the breaking is a frozen *response*: in
the exact theory the core follows each M_L component and the
degeneracy survives; freezing the core at the reference's polarization
cuts that loop. Hence the splitting scales as (valence quadrupole) ×
(core quadrupolar polarizability): no polarization functions — closed
channel, accidental degeneracy; every shell added reopens it (the
verification program, Results). The ratio splitting/ε is the effective
coupling — an atomic property. Commutators and splittings are the
detectors; a numerically pretty ⟨L̂²⟩ is not.

Throughout, a treatment is examined through three questions:

1. does every diagonalized root carry a sharp term label (⟨L̂²⟩ =
   L(L+1), ⟨Ŝ²⟩ = S(S+1))?
2. how large are the commutators ‖[Ĥ, L̂z]‖, ‖[Ĥ, L̂²]‖, evaluated as
   exact matrices over the determinant basis?
3. does the term table (state counts per term) equal the combinatorial
   free-ion table?

And for every remedy we report how it achieves this, and the energetic
cost relative to the broken-symmetry solutions.

### Notation and operators

We work with spin-orbitals (α = 0..n−1, β = n..2n−1) and
**determinants** — sorted occupation tuples. A **sector** fixes the
numbers (n_α, n_β) of α/β electrons, hence M_S = (n_α−n_β)/2 and N. A
**member** is one m-function of a shell (px, py, pz; the d and f sets
below). The operators N̂, L̂z, L̂², Ŝz, Ŝ² all preserve sectors (each
commutes with Ŝz and N̂), so every matrix lives inside one sector of
fixed M_S.

Every operator enters the determinant space as a **one-body part** and,
where it has one, a **two-body part**,

  O₁ = Σ_pq o_pq a†_p a_q,   O₂ = ½ Σ_pqrs g_pqrs a†_p a†_q a_s a_r,

with these symbols (spin-orbitals ordered (α; β); L̂^a = C† L̂^a_AO C,
a = x, y, z — the AO angular-momentum matrices transformed into the
window):

```
N̂  :  o = 𝟙                                          (counts each electron once)
L̂z :  o = L̂^z
L̂² :  o = Σ_a (L̂^a)²
      g = 2 Σ_a (L̂^a)_pr (L̂^a)_qs
Ŝz :  o = diag(½·𝟙, −½·𝟙)
Ŝ² :  o = ¾·𝟙
      g = 2·[(ŝ_z)_pr (ŝ_z)_qs + ½(ŝ₊)_pr (ŝ₋)_qs + ½(ŝ₋)_pr (ŝ₊)_qs]
L̂± :  o = L̂^x ± iL̂^y                                (orbital ladders; sector operators)
Ŝ± :  o = ŝ±   (ŝ₊ = |α⟩⟨β|, the β→α flip)           (inter-sector maps: M_S → M_S ± 1)
Ĥ  :  o = h_eff  (the folded core field, above)
      g = (pq|rs)   (the active-orbital Coulomb integrals)
```

The ladder algebra is pinned elementwise against the delivered
operators: L̂² = L̂z² + ½(L̂₊L̂₋ + L̂₋L̂₊), Ŝ² = Ŝz² + ½(S₋S₊ + S₊S₋),
[L̂z, L̂±] = ±L̂±. Ŝz carries no certificate — [Ĥ, Ŝz] ≡ 0 identically;
it only labels the sector (its matrix is exactly M_S·**1**).

**ε — the core anisotropy** — is the relative weight of the non-scalar
part of the frozen core density, ‖ρ_c − ρ̄_c‖_F / ‖ρ̄_c‖_F (ρ̄_c = the
spherical average; measured directly from the SCF density by
`core_anisotropy()`, no CI involved). The meters, and what each is
blind to:

| Diagnostic | What it measures | Sensitive to | Blind to |
|---|---|---|---|
| ‖[Ĥ, L̂z]‖_F, ‖[Ĥ, L̂²]‖_F | the assembled treatment (window + core fold) | everything, order ε | vacuous in single-irrep sectors (L̂²); ±pair structure (L̂z) |
| ‖[P, L̂z]‖, ‖[P, L̂⊥]‖ | the window span alone | member asymmetry, missing members; no optimizer noise | — (as a pair) |
| ⟨L̂²⟩, ⟨Ŝ²⟩ per root | term identity | order ε² | small defects |
| multiplet splitting | fold anisotropy | order ε | — |
| orthonormality error | validity of the stage's MO matrix (C†SC = 1) — a guard, not a symmetry diagnostic | construction bugs | — |

### The optimizer side

Over m-labeled determinants Ĥ_CAS is block-diagonal in M_L: the fold is
m-diagonal, and the Coulomb term conserves total M_L (Δm₁ + Δm₂ = 0).
Hence, for a window deformation (a **member tilt** = one member
replaced by a radial admixture, cos θ·(n l m) + sin θ·(n′ l m)): in
one-dimensional M_L blocks the tilted member drops out of every lowest
determinant, so the energy is *flat* (a **flat direction** — a
deformation the energy cannot see) — breaking m is free, the endpoint
a *gauge* (energetically indistinguishable descriptions); in
multi-dimensional blocks the tilt splits the degenerate block
symmetrically about its trace, so the lowest level *descends*
(Jahn–Teller-like) — the optimizer is driven to break M; and the state
average is stationary with zero curvature — no restoring force. The
frozen-scan verdicts (deterministic, no optimizer; curves in
`elements/*/NOTES.md`): on p² the entire CASSCF gain is one m-symmetric
radial rotation at exact closure (C −250.8, Si −126.7 mEh) and tilts
are exactly flat — breaking never buys energy; on d atoms radial
rises, tilts descend — and the large free-run gains are
window-identity changes instead: the **escape** (window rebuilt from
orbitals the configuration never implied; states describe no term) and
the **swap** (a valence 4s for a core 3p — a different core/active
partition; totals only). Every CASSCF start also carries the
**start-window defect**, ‖[P_w̃, L̂]‖ ≲ ⟨window|core⟩·‖[P_core, L̂]‖,
from core-orthogonalizing against a non-closed core — zeroed by the
purified substrate (Methods).

### The failure inventory

| Level | Symptom (order of magnitude) | Mechanism | Fix | Cost |
|---|---|---|---|---|
| 1. Basis bookkeeping | labels mislead; shells fuse | wrong inputs → wrong constructions | label + metric discipline | none |
| 2. Sector witnesses | ‖[Ĥ,L̂²]‖ = 0 with defects present | single-irrep sectors make L̂² = L(L+1)·**1**; ±pair blind spots | projector pair; multi-term sectors | sector size |
| 3. Core field | splittings and commutators 1e-5…1e-1 | the polarization chain (above) | fold cures / purification (Methods) | ΔE per Results |
| 4. Orbital optimization | tilted windows, floors, escape/swap | flat and descending directions; identity change | restrictions + substrate (Methods) | 0…~30 mEh above the deepest minimum |
| 5. SCF reference | ~1e2 mEh wrong configuration; ~10 mEh metastable basin | energy-ordered filling; degenerate-open-shell basins | the Methods toolbox | per-atom reference data |

**Resolutions.** Every cure removes ε at its source: (i) state
averaging over rotation-complete state sets — the polarizations of all
M_L components average to a spherical core, and every state of the
average shares one scalar Ĥ_CAS; (ii) term-adapted multi-determinant
references; (iii) promoting the core into the active space; (iv)
shell-complete active windows (m-complete by construction); (v)
verification — ⟨L̂²⟩, the ±M_L pairing, the L̂₊ tests and the stabilizer
test (the verification program, Results). Class mapping: *exact* —
raw-window CASCI, frozen-window rungs on the isotropic substrate;
*baseline* — the start-window defect and the flat-direction floor;
*broken* — the escape and the swap.

---

## Methods

**ROHF — the reference that defines the core.** Everything downstream
folds it. Where energy-ordered filling fails (a doubly-occupied-d
runaway ~1e2 mEh wrong; a metastable polarized basin ~10 mEh up) it is
kept honest by character-forced occupations (orbital type measured by
S-projector weight, not energy order), an energy gate over cycled
starting guesses, and cached trusted references — E_ref and the AO
spans are the atom's data. It is also the source of ε; the fold
variants are raw, averaged, fractional. Records:
`elements/*/NOTES.md`.

**CASCI — at fixed core and window.** Diagonalizes one sector; the
window is the raw shell AOs (rotation-closed by construction). The core
enters through the fold (raw, averaged, fractional).

**CASSCF — letting the orbitals move.** Minimizes one root over orbital
rotations — the relaxation CASCI cannot express. Run as a grid of
*starts × restrictions*: the **substrate** (the starting MO matrix:
block-safe, the window core-orthogonalized — the leak guard of the
partition exactness, Motivation) and **masks** on rotation pairs:

1. window frozen (correlation through the core alone),
2. radial restriction (same-(l, m) chains only — the
   symmetry-preserving channel),
3. free.

A **basin** is the region of attraction of one local minimum; the free
landscape has several. Start and converged windows are compared by the
**span overlap** (smallest principal cosine between the two spans —
did the window become a different space?).

**Precaution (both objectives).** The free mask can rotate a core
orbital into the window — the partition then changes, totals stop
being comparable, and our fold's re-evaluation of such an endpoint is
not an energy of the original problem (the bookkeeping double-counts
the moved core orbital). The run's own total is quoted and the row is
compared with care — flagged by the fold shift.

Naming bridge for the tables: the frozen-window mask is what the tables
call **core-relax** — the window span is untouched (closure preserved
whatever the core), only core↔virtual rotations re-optimize the core's
description. Its purified twin (**core-relax, pur**) runs the same mask
on the purified core substrate, so the run inherits a scalar fold — the
certified class. Single-root CASSCF stays gauge-chaotic even on the
purified substrate (its breaking directions are its own), which is why
the tables report the plain core-relax row for CASSCF and the purified
one for SA-CASSCF — the state average restores what a single root
cannot.

**SA-CASSCF — the multiplet average.** Replaces the single root by the
average over the ground term's M_L components — the rotationally scalar
object. It anchors the term identity but not the window. The certified
variants:

1. frozen window,
2. frozen window on the **purified core substrate** — the core span
   rebuilt as complete m-sets on shared radials, exactly isotropic.

**Scans — the frozen landscape.** A diagnostic tool, not a deliverable:
freeze everything, deform one window coordinate, re-diagonalize the
CASCI — deterministic, no optimizer. This shows *where the optimizer's
gain lives*: on the m-symmetric radial channel (symmetry-preserving) or
on member tilts (symmetry-breaking), and whether tilts are flat
(breaking is free) or descend (breaking is driven). Per-atom curves:
`elements/*/NOTES.md`.

**Running.** Every stage checkpointed (chkfiles cached and trusted);
tests in two tiers (`-m fast`, sub-second; the full suite ~13 min); a
new atom is per-atom data, not code (the pipeline, below).

| Standard approach | Expected result | Explanation | Alternative |
|---|---|---|---|
| CASCI over canonical virtual MOs | closure defects at the 1e-1 level; term labels unreliable | canonical MOs diagonalize a Fock built from an anisotropic density; quasi-degenerate same-l members mix arbitrarily | raw shell windows, averaged fold |
| free CASSCF, default start | p²: full gain with a silently member-tilted window; d-shells: escape/swap basins, ⟨L̂²⟩ off by units | flat directions (endpoint = noise) and descending directions; window identity free to change | restrictions + substrate control |
| free SA-CASSCF | ⟨L̂²⟩ ≈ L(L+1) restored, closure not (Sc: 6.002 with ‖[H,L̂z]‖ ≈ 0.8) | the average is scalar to first order only; no restoring force on flat directions; other basins still reachable | frozen-window SA on the purified substrate |

### The pipeline

A new atom is a declared case, executed top-down — per-atom *data*
(closed/open AO spans, the reference energy E_ref, the combinatorial
term tables), not per-atom code:

```
CaseSpec — the declaration
│  atom, basis/ECP, spin, nelec_act, ncore, shell_patterns
│  ref: ReferenceSpec (occ spans l0..l3 | label prefixes,
│        conv_accelerator ∈ {"diis","newton","diis+newton"}, occ_e_ref)
│  want: {terms, casci-variants, casscf, sa-casscf, treatments, scans}
│
└─ run_case(spec)
   ├─ 1. ROHF reference (reference.py)
   │      a cached checkpoint is reused as is; a fresh run applies
   │      character-forced filling, the declared acceleration, and
   │      retries guesses until the gate |E − E_ref| is satisfied.
   ├─ 2. CASCI per sector — fold ∈ {averaged, raw, fractional}
   │      term tables produced and compared to the combinatorial
   │      contract; deviations are reported.
   ├─ 3. if "casscf" → single-root CASSCF
   ├─ 4. if "sa-casscf" → SA-CASSCF over the ground term's M_L set
   ├─ 5. if "treatments" → the declared window × treatment product
   │      (study.py, on the treatments.py registry): the windows are
   │      walked for frozen CASCI as well as for the optimizer
   │      treatments.
   └─ 6. if "scans" → radial / tilt scans (scans.py)
```

The reference gate's energy `occ_e_ref` is self-calibrated atom data:
several starting guesses are tried, their solutions classified by
character (the S-projector weight), and the passing basin's energy
recorded.

### The protocol card

**Run all five deliverables, report all five** — the certified treatment
and the direct values are complements, not alternatives. Each deliverable
is a named configuration in the Methods space; each optimizer
deliverable is reported as the pair *free variant / symmetric variant*.

| Deliverable | Configuration | Purpose | Acceptance |
|---|---|---|---|
| ROHF reference | forced filling + gate | the anchor every ΔE hangs on | occupations & core purity by projection; E = E_ref |
| CASCI anchor | raw shell window, averaged fold | the exact term tables | tables = combinatorial; commutators exact |
| Certified treatment | SA-CASSCF + window frozen + purified core | exact-class symmetry with correlation | commutators & window spans exact; tables intact |
| Direct CASSCF | free, shell substrate | the energy-seeking value + basin anatomy | classified (span overlap, window weights); totals quoted |
| Direct SA-CASSCF | free SA-CASSCF, shell substrate | best direct value; term-readable even in the swap basin | term label; basin classified; totals quoted |

For a new atom: (1) run the reference (spans and E_ref are the atom's
data); (2) the CASCI anchor — term tables must come out exact (Results);
(3) the certified treatment (expect the exact class); (4) both direct
runs, classified; (5) optionally, frozen scans for the radial/tilt
anatomy.

---

## Results

All raw-shell windows and the averaged fold unless stated; every number
is this project's measurement. **Every other atom lives in
`elements/<Atom>/TABLES.md`** (+ `NOTES.md` for the dated specifics;
C and Si mirror each other throughout, Sc and U carry the
transition-metal and f-scale stories); the artifacts are under
`elements/Ti/` (regeneration: Files, below). Ti: ROHF reference
−848.406489 Eh (3d²4s², the gated reference), M_S = 1, ground term ³F.
On symmetry-broken basins the term table itself scrambles — independent
evidence beyond the commutators.

### Term table — the contract each sector must meet

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,5) / {3d} | Ms=1 | 10 | ³F: 7, ³P: 3 |
| CAS(2,5) / {3d} | Ms=0 | 25 | ³F: 7, ¹D: 5, ³P: 3, ¹G: 9, ¹S: 1 |
| CAS(4,6) / {3d,4s} | Ms=1 | 120 | ³D: 15, ³F: 35, ³G: 18, ³H: 22, ³P: 15, ⁵D: 5, ⁵F: 7, ⁵P: 3 |

The combinatorial free-ion content — the exact stage reproduces it
state-for-state (the Ms=1 block comes out exact in the showcase at the
top of this note).

### Multiplet energies

CAS(2,5), averaged fold; exact within each term to ≤1e-9; the M_L
column from the sector Lz spectra. The ordering matches the free-ion
(NIST) tables; the gaps are CAS-level values, not experimental ones.

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ³F | 7 | −1.468663 | — | −3…3 |
| ¹D | 5 | −1.423747 | +0.0449 | −2…2 |
| ³P | 3 | −1.413969 | +0.0547 | −1…1 |
| ¹G | 9 | −1.398223 | +0.0704 | −4…4 |
| ¹S | 1 | −1.295050 | +0.1736 | 0 |

### Methods compared

Every row: E — the run's variational total; ΔE vs the ROHF anchor;
**cross-check** = |E(analysis fold) − E(as run)|, max over roots (the
mEh-class averaging effect; the fold convention is machine-exact,
1e-15 vs PySCF's get_h1eff). Verdicts (the classes above): **exact** —
commutators and span certificates at the 1e-13 class; **baseline** — a
named non-closure mechanism at its signature scale (the raw fold's
percent-class commutators at unchanged energy make it the
polarized-core detector); **escape / swap** — the window-identity
basins, compared in totals only. The certified rungs sit mEh-tight
around the frozen value (CASCI +10.5 ≈ core-relax +10.4 ≈ SA-pur
+11.7): orbital optimization buys ≈nothing on the as-run window — the
negative free numbers are window-identity changes, classified by the
span certificates; the SA swap dives −861 mEh in its own fold — the
851 mEh cross-check is exactly that convention break. Full diagnostics
(span certificates, configurations): `elements/Ti/TABLES.md`; the
cross-atom summary: the element files.

| treatment | E, Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ | cross-check (mEh) | verdict |
|---|---|---|---|---|---|---|---|
| ROHF (anchor) | −848.4065 | 0 | — | — | — | 0 | — |
| CASCI, averaged | −848.3959 | +10.5 | 12.0000 | 4.2e-09 | 2.3e-08 | 2.8 | exact |
| CASCI, raw | −848.3959 | +10.5 | 12.0000 | 4.9e-02 | 2.6e-01 | ~0 | baseline (detector) |
| CASCI, fractional | −848.3585 | +48.0 | 12.0000 | 1.1e-09 | 7.6e-09 | ~0 | exact |
| CASSCF, free | −848.4377 | −31.3 | 0.0503 | 4.3e-01 | 1.1e+00 | 2.8 | escape — broken |
| CASSCF, core-relax | −848.3961 | +10.4 | 12.0000 | 4.2e-09 | 2.3e-08 | 2.7 | baseline |
| SA-CASSCF, free | −849.2676 | −861.1 | 12.0029 | 1.1e+00 | 6.0e+00 | 851 | swap — broken |
| SA-CASSCF, core-relax, pur | −848.3948 | +11.7 | 12.0000 | 1.0e-13 | 4.0e-13 | 0.0 | exact |

### The window lattice — window kinds at matched determinant cost

Three window kinds (Motivation), same mean field, compared at **equal
determinant count** (equal cost). Every row carries the span
certificate ‖[P,Lz]‖ / ‖[P,L⊥]‖: the AO and NO windows close at
~1e-13; the raw-MO spans never close (O(1)–9 at every size measured,
4–100 orbitals — the SCF's open orbitals are m-delocalized), so the
certificate disqualifies them without running any CI.

| window | CAS | dets | AO, Eh | ⟨L²⟩ | ROHF MO, Eh | ⟨L²⟩ | NO |
|---|---|---|---|---|---|---|---|
| 3d | (4,5) | 50 | −848.054 | 6.000 | — | — | — |
| 3d4s | (4,6) | 120 | −848.395 | 12.000 | −848.405 | 0.000 | −848.399 |
| 3d4d | (4,10) | 1,200 | −848.232 | 6.000 | — | — | — |
| 3d4s4p | (4,9) | 756 | −848.429 | 12.000 | — | — | −848.432 |
| 3d4s4d | (4,11) | 1,815 | −848.409 | 12.000 | −848.414 | 0.010 | — |
| 3d4s4p4d | (4,14) | 5,096 | **−848.438** | 12.000 | −848.420 | 0.201 | **−848.438** |
| 3d4s4p4d5d | (4,19) | 18,411 | −848.439 | 12.000 | — | — | — |
| 3d4s4p4d4f | (4,21) | 27,930 | **−848.446** | 12.000 | — | — | — |

What it says, in one pass: the AO lattice descends monotonically below
ROHF at exact ⟨L²⟩ (the d-atom window already correlates); the MO
column wins only the smallest pair — its window contains the reference
determinant (representability) — then saturates with ⟨L²⟩ never
quantizing; the NO column matches the AO value at equal dets and
exceeds it at the top. Both transition metals cross below their
escaped CASSCF basins at exact symmetry — the escape ceilings are
window statements, not method limits. U: the same at f-scale, with a
reference determinant that is a numerically pure ⁵S° state
(`elements/U/`).

### Files and regeneration

**The module map.**

| Module | Owns |
|---|---|
| `system` | system validation — one atom at the origin, spherical basis (L̂ defined about r = 0) |
| `angular` | one-electron L̂x/L̂y/L̂z (+ L̂± ladder) AO matrices — two independent constructions pinning each other (libcint; exact tesseral algebra) |
| `detbasis` / `secondquant` | many-body matrices over determinants — Slater–Condon and explicit a†/a, elementwise-identical twins |
| `operators` | operator assembly — `build_operators`: the symbols and the ladder operators, the window-closure certificates (`window_closure`) |
| `reference` / `casscf` | the engine — SCF reference with gates; the one deep optimizer (`constrained_casscf`) with its masks and adapters |
| `windows` / `core` | window construction, selection, the kind registry (`Window.columns`); the purified core substrate and fold energy |
| `natural` | NO mathematics — densities, occupation groups, NO-window selection (the NO columns of the lattice) |
| `stages` | the one uniform analysis pass per finished stage — the source of every number in this section |
| `scans` | frozen-landscape measurements — deterministic |
| `analysis` | term labels, configuration summaries, report rendering |
| `cases` / `treatments` / `study` | orchestration — the declared per-atom campaign (`cases`: CaseSpec → run_case); rungs as data (`treatments`); the declared window × treatment study runner (`study`) |

**The files.** Everything human-readable lives committed in
`elements/<Atom>/`: `TABLES.md` (term content, spectrum, methods
compared, window lattice, study inserts), `report.md/json` (the term
tables), `ladder.json` (the full rung registry — Methods compared),
`scan.md`, `sa_decomposition.md`, `basin_map.md`, `NOTES.md` (dated
specifics). New atom = new element file, not README growth; the
regenerators write straight into these files (`scripts/regen_ladders.py
Ti` -> `elements/Ti/ladder.json`; chkfiles land next to them and are
gitignored). Regeneration: term tables — `scripts/run_c.py` and the
per-atom `run_atom` calls (recorded in each report header); rungs — the
per-atom `run_ladder` calls; scans — `scripts/radial_scan.py`; the
lattices — `scripts/run_tm_lattice.py`, `run_tm_methb.py`, `run_u.py`,
`run_u_methb.py`, `run_u_sa.py`. Executable documentation: the test
suite (128 tests) — term tables, exact-class closure pins, basin
branches. Plain-language companion: `RUNGS.md`.

### Status and directions

**Status.** C, Si, Sc, Ti delivered end-to-end; U campaign complete
through ticket 23 (the relaxation-vs-growth verdict: the deep
correlation is radial relaxation, held at exact ⁵L by the ensemble,
destroyed by free single-root). Per-atom records: `elements/<Atom>/`;
cross-atom story: the element files.

**Open physics.** The ⁵S° mechanism (U reference + free-run attractor;
Ti ⟨L²⟩ = 0 MO windows — same class). Inner core-pair valence (the
m-complete design; the accidental 6p probe moved ~−0.5 Eh). U
f³d¹-vs-f²d² ordering inside one NO window. NIST comparison.

**Method gaps (ticket 20).** Axes-product runners behind a
`Window.columns(mol, mf)` resolver; per-cell budgets for long runs;
`--fresh` / chkfile key versioning; the intent manifest; the
full-matrix slice guard; the needle linear-dependence helper (fired
once); promote `closed_shell_core_energy` to the library.

**Forward directions.** Spin–orbit (splits the certified multiplets;
extends in place: Ĵ², Ĵz as determinant-basis matrices). Transition
moments / g-factors (between the certified eigenstates). Scalar
relativity (h₁ only; re-measure the reference basins once). Further
atoms = element files + data rows.

**Verification program (planned).** The Motivation narrative is
checkable end to end; the ladder operators make part of it free. In
order: (1) *toy axial p²* — the 15 p² determinants with exact shell
Coulomb plus a one-body axial field δ_m ∝ R₂(2−3m²)/5: the Λ
branching, the 3R₂/5 Σ–Π gap, ±M_L pairing, ¹S shifting only at R₂²,
¹D–¹S second-order repulsion, ‖[Ĥ,L̂z]‖ ≈ 0 with ‖[Ĥ,L²]‖ ∝ |R₂| (pure
numpy, no SCF); (2) *operator-level* — L̂₊ tests (pure ³P roots give
the exact ladder amplitudes; ¹D purity via its deficit), the variance
⟨L̂⁴⟩−⟨L̂²⟩², and the stabilizer test (rotate the frozen core, rebuild
the fold: the spectrum is invariant exactly for rotations in
Stab(ρ_c)); (3) *real-code headline* — the core quadrupole and the
per-determinant Σ/Π split tracked as polarization shells enter the
basis one at a time (the accidental-degeneracy corollary), SA-15
restoration, state-specific re-breaking from a spherical start;
(4) *response toy* — a core orbital free to mix λ·d₀: λ ∝
α_core·Q_val, the shielding sign, λ = 0 for spherical valence
occupations.
