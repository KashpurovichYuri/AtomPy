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
pytest                        # 184 tests, ~13 min (chkfile caching); -m fast = the seconds tier
python scripts/run.py C       # writes elements/C/ (report.md/.json + chkfiles)
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
rotation-closed. Every atom delivers this per Ms sector; any atom's
tables: `elements/<Atom>/TABLES.md`. The Ŝz line documents the
sector's M_S; it certifies nothing, since [Ĥ, Ŝz] ≡ 0 identically
(Motivation). Details: Results (Term table, Multiplet energies);
regeneration: Results (Files).

**Treatment rungs** — every orbital treatment priced against the
atom's ROHF anchor, the free variant beside the symmetric one, each
classified by its diagnostics. The manifest: Results (Methods
compared); the per-atom tables — `elements/<Atom>/TABLES.md`.

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
  potential of the core density, plus nonlocal exchange).
- **(tu|vw)** — four-center two-electron integrals among active
  orbitals.

CASCI is the exact diagonalization of this effective Hamiltonian; its
only approximation is the frozen mean field of the inactive space, so core–active interaction is kept static, core–active correlation is dropped. The
partition is **exact iff the window is orthogonal to the core
orbitals**: a window vector leaning into the core, w = w⊥ + w∥ with
w∥ ∈ span(core), admits configurations with extra density in *occupied*
core orbitals that the antisymmetrized problem forbids — the CI
silently harvests deep-core stabilization the Pauli principle would
have priced. Every window is therefore built core-orthogonal.

**Injection point 1 — the fold.** The active–active part (tu|vw) is a
rotational scalar; **the folded part f carries whatever symmetry the
core density ρ_c has**. ROHF's ρ_c is axial (e.g., for Carbon atom) → 
f is axial → Ĥ_CAS is axial: the first entry of symmetry breaking. The same
is seen through orbital rotations: the valence open shell exerts a
quadrupolar (rank-2, m = 0 about the axis) field on the core, and
first-order perturbation theory mixes the s core orbital with what the
selection rule allows — l → l ± 2, m fixed — the d shell's m = 0
member d₀ (= dz²). The optimized core orbital sits in its polarized
shape φ_c ≈ 1s + λ·d₀ (λ = the polarization amplitude); rotating by
90° about x turns the needed admixture into a dz²/dx²−y² combination —
a d direction the frozen space {φ_c, 2s} (an s span plus that one d₀
direction) does not contain. Since the frozen core sits inside Ĥ_CAS
through f_tu, the leak is inherited:

    U(R)⁻¹ Ĥ_CAS U(R) = Ĥ_CAS[core on a new axis] ≠ Ĥ_CAS unless R is
    about z (in case of Carbon atom).

The symmetry group of Ĥ_CAS is the stabilizer of the frozen core
density (e.g., rotations about the torus axis, when such an axis exists as for
Carbon atom). The cure is a symmetric ρ_c; three core treatments are compared
in the project:

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

**The worked instance — carbon 2p²** (illustration). ROHF places
the two ³P electrons in {px, py}: the open-shell
density is a torus, axial about z. The ³P multiplet branches by |M_L|
into Π (M_L = ±1) and Σ (M_L = 0) — exactly the σ/π branching of a
linear molecule or a ligand field, measured on carbon (the ³P Σ/Π gap,
the chain table below); a systematic cross-check against known
crystal-field splittings belongs to the verification program (Status
and directions). Other references have other, generally smaller
stabilizers — hence measured raw-fold commutators from the 1e-5 class
(C, an almost perfect axis) to percent class (Ti).

**Injection point 2 — the window.** With P the projector onto the
active configuration space, the exact commutator [Ĥ, L̂] = 0 gives

    [PĤP, PL̂P] = PL̂(1−P)ĤP − PĤ(1−P)L̂P,

and if the space and its complement are rotation-closed ([P, L̂] = 0), the
projected generators still close the algebra and Schur survives. This
is what the code checks: ‖[P, L̂z]‖ and ‖[P, L̂⊥]‖ (`window_closure`).
Carbon shows why the window must span the whole valence shell
{px, py, pz} — rotating py about x gives cos θ·py + sin θ·pz, which
stays inside the p span; a truncated set leaks. The active set — the
**window** — is therefore always chosen shell-complete and built
core-orthogonal: closure is imposed by design, never converged to.
Starts come in three kinds throughout, compared at matched determinant
cost (Results); their spans are the windows: **m-complete AO shells
and bands** — rotation-closed by design; **raw-ROHF MO picks** —
energy order slices through the degenerate open shell, whose members
are m-delocalized, so these never close (occupancy-guided picks from a
symmetric reference — the fractional SCF, whose exactly degenerate
shells are flagged by occupation — take whole blocks and close as
naturally as shell spans); and **natural-orbital (NO) occupation
groups** — NOs of the term-averaged density, which by Schur's lemma
arrive as complete m-sets: rotation-closed by construction.

**The chain, measured.** ε (the core anisotropy; defined with the
meters, below) against the splitting it produces (raw fold):

| Atom | ε (relative) | splitting (Eh) | splitting/ε (Eh) | ⟨L̂²⟩ of the ground term |
|---|---|---|---|---|
| C | 2.405e-3 | 4.373e-5 (³P) | 0.02 | 2.0000 |
| Si | 2.109e-3 | 2.431e-4 (³P) | 0.12 | 2.0000 |
| Ti | 7.754e-3 | 4.144e-3 (³F septet) | 0.53 | 12.0000 |
| Sc | 1.259e-2 | 1.857e-3 (²D quintet) | 0.15 | 6.0001 |
| U | — (not measured; U runs on the pur substrate) | | | |

Two readings of one measured number — the splitting:

- *The labels survive; the degeneracy dies.* The fold is spin-free, so
  different-S terms never mix (¹D ↔ ¹S only, never ³P — the spin-free
  commutation, above): ⟨L̂²⟩ stays exact while the M_L degeneracy
  branches. The splitting, not the label, is the signal.
- *A frozen response.* In the exact solution the core repolarizes for
  each M_L component and the degeneracy survives; freezing ρ_c at the
  reference's shape cuts that loop. The splitting scales as valence
  quadrupole × core quadrupolar polarizability — zero without
  polarization functions (closed channel, accidental degeneracy),
  reopened by every polarization shell (the verification program,
  Status and directions). Splitting/ε is the effective coupling, an
  atomic property.

The detectors are the commutators and the splittings; a numerically
pretty ⟨L̂²⟩ is not.

Every remedy reports how it achieves this, and the energetic cost
relative to the broken-symmetry solutions.

### Notation and operators

We work with spin-orbitals (α = 0..n−1, β = n..2n−1) and
**determinants** — sorted occupation tuples. A **sector** fixes the
numbers (n_α, n_β) of α/β electrons, hence M_S = (n_α−n_β)/2 and N. A
**member** is one m-function of a shell (px, py, pz; the d and f sets
below). The operators N̂, L̂z, L̂z², L̂², Ŝz, Ŝ² all preserve sectors (each
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
L̂z²:  o = (L̂^z)²
      g = 2·(L̂^z)_pr (L̂^z)_qs
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
| ⟨L̂z²⟩ per root | the root's M_L² content | broken folds: each branch's composition | isotropic folds: degenerate-cluster gauge (starred) — not a per-root observable |
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
average (SA-CASSCF, Methods) is stationary with zero curvature — no
restoring force. The frozen-scan verdicts (deterministic, no optimizer;
curves where measured, `elements/<Atom>/scan.md`): on p² the entire
CASSCF gain is one m-symmetric radial rotation at exact closure and
tilts are exactly flat — breaking never buys energy; on d atoms radial
rises, tilts descend — the drive behind their broken basins. The large
free-run gains are window-identity changes instead: the **escape** —
the optimizer rebuilds the window from orbitals the configuration
never implied, and the converged states match no term; the **swap** —
a valence orbital exchanged for a core one, a different core/active
partition (totals only, Partition bookkeeping). Which basin a free run
finds is start-dependent; instances and verdicts: Results (Methods
compared) and the element files. Every CASSCF start also carries the
**start-window defect** (closure lost at core-orthogonalization
against a non-closed core) — zeroed by the purified substrate
(Methods).

### The failure inventory

Two injection points carry every symmetry failure — the fold (IP 1)
and the window (IP 2); the SCF reference sits upstream of both, and
the instrument (labels, sectors, orthonormality) can lie on its own.

| Where | Failure | Symptom | Mechanism | Cure |
|---|---|---|---|---|
| the instrument | basis bookkeeping | labels mislead; shells fuse | wrong inputs → wrong constructions | label + metric discipline |
| the instrument | sector witnesses | ‖[Ĥ,L̂²]‖ = 0 with defects present | single-irrep sectors make L̂² = L(L+1)·**1**; ±pair blind spots | projector pair; multi-term sectors |
| upstream of both | SCF reference | ~1e2 mEh wrong configuration; ~10 mEh metastable basin | energy-ordered filling; degenerate-open-shell basins | character-forced occupations + energy gate (Methods) |
| **IP 1 — the fold** | core field | splittings and commutators 1e-5…1e-1 | the polarization chain (above) | averaged fold / SA / purification (Methods) |
| **IP 2 — the window** | orbital optimization | tilted windows, floors, escape/swap | flat and descending directions; identity change | shell-complete starts + masks (radial / frozen — imposed on CASSCF and SA alike) + substrate (Methods) |

**Resolutions.** Every cure removes ε at its source, at its injection
point: (i) state averaging over rotation-complete state sets (IP 1) —
the polarizations of all M_L components average to a spherical core,
and every state of the average shares one scalar Ĥ_CAS; (ii)
term-adapted multi-determinant references (upstream); (iii) promoting
the core into the active space (dissolves IP 1); (iv) shell-complete
active windows, m-complete by construction (IP 2); (v) verification —
⟨L̂²⟩, the ±M_L pairing, the L̂₊ tests and the stabilizer test (the
instrument; the verification program, Status and directions). Class
mapping: *exact* — raw-window CASCI, frozen-window rungs on the
isotropic substrate; *baseline* — the start-window defect and the
flat-direction floor; *broken* — the escape and the swap.

---

## Methods

**ROHF — the reference that defines the core.** Everything downstream
folds it. Where energy-ordered filling fails (a runaway to a wrong
closed-shell configuration, ~1e2 mEh too high; a metastable polarized
basin, ~10 mEh up) it is
kept honest by character-forced occupations (orbital type measured by
S-projector weight, not energy order), an energy gate over cycled
starting guesses, and cached trusted references — E_ref and the AO
spans are the atom's data. It is also the source of ε; the fold
variants are the three core treatments above.

**CASCI — at fixed core and window.** Diagonalizes one sector; the
window is a stated *start*'s span — AOs, MOs, or NOs (the start axis,
below); closure is a property of the span. The core enters through the
fold (the three treatments above).

**CASSCF — letting the orbitals move.** Minimizes one root over orbital
rotations — the relaxation CASCI cannot express. Run as a grid over
*starts × masks × substrate*: the **substrate** is the starting MO
matrix (block-safe, the window core-orthogonalized — the leak guard of
the partition exactness, Motivation); **masks** on rotation pairs:

1. window frozen (correlation through the core alone),
2. radial restriction (same-(l, m) chains only — the
   symmetry-preserving channel),
3. free.

A **basin** is the region of attraction of one local minimum; the free
landscape has several. Start and converged windows are compared by the
**span overlap** (smallest principal cosine between the two spans —
did the window become a different space?).

**The start axis.** The columns that define the window span — the same
axis CASCI's windows and the optimizer grid use: *AOs* (shell AOs in
the optimizer-runnable orthogonalized form), *MOs* (canonical MOs with
a stated picking rule; the stock default is PySCF's low-energy slice,
other picks admissible per slot; closure is a property of the pick —
Motivation), *NOs* (occupation groups of the term-averaged CI — the
default start, rotation-closed by construction); *virtual-mix* appears
only as a registry diagnostic. The manifest tables read mechanically:
a blank factor cell inherits from the row above, and a pending row is
a planned run.

**Partition bookkeeping (the free-mask risk, made precise).** The
quantity at risk is the energy itself: a folded total is an energy *of
a stated partition*. For a partition (C, W) — core columns C, window
W —

    E(C, W) = E_core[C] + min_Ψ∈W ⟨Ψ| Ĥ[h₁(ρ_C)] |Ψ⟩,

exactly (Slater–Condon: core-index integrals fold into
h₁(ρ_C) = h + J[ρ_C] − K[ρ_C] — our fold, `build_operators`; pairs in
C into E_core[C]; machine-exact against PySCF's own fold `get_h1eff`,
1e-15 — the pin lives in test_operators). Our analysis splits
`mo[:, :ncore] | mo[:, ncore:]` (`stage_facts`) — valid iff that split
equals the run's partition. A free mask may exchange core and window
orbitals; the kernel's final canonicalization then reorders columns by
energy, and split ≠ partition: the moved orbital φ is folded into ρ_C
*and* correlated inside W — double-counted. Such a number is an energy
of no partition and reads *below* the kernel's own total. The **fold
shift** catches it: ~0 consistent; mEh the averaging effect;
large-and-below partition changed (compare kernel totals only).

**SA-CASSCF — the multiplet average.** Replaces the single root by the
average over the ground term's M_L components — the rotationally scalar
object. It anchors the term identity but not the window. The certified
variants are similar to CASSCF theory.

**Scans — an optional, experimental diagnostic.** A tool, not a
deliverable: freeze everything, deform one window coordinate,
re-diagonalize the CASCI — deterministic, no optimizer; the anatomy of
the optimizer's gain (Motivation, above). Curves where measured:
`elements/<Atom>/scan.md`.

**Running.** Every stage checkpointed (chkfiles cached and trusted);
tests in two tiers (`-m fast`, sub-second; the full suite ~13 min); a
new atom is per-atom data, not code (the pipeline, below).

**The treatment grid.** Not a menu of standard approaches — one
declared grid, priced and classified. Holding the atom's window fixed,
the axes are:

- the **fold** — raw / averaged / fractional (the core field; IP 1),
- the **start** — AOs / MOs / NOs (the columns defining the window
  span; IP 2),
- the **mask** — free / radial / frozen (which rotations the optimizer
  may use; IP 2, imposed on CASSCF and SA alike),
- the **objective** — the single root or the ground term's M_L average,
- the **substrate** — polarized / purified (the core span under the
  fold; touches both injection points).

Every cell is priced against the ROHF anchor and classified by the
diagnostics; the grid localizes each failure by which axis flips it,
and the certified cells prove restoration at a stated cost. Any run
declares a **subset** of cells — the spine below is the default.

**Which cells are needed.**

- *The spine (default subset, every atom):*
  - the three folds;
  - frozen × purified under both objectives — the certified class;
  - one free row per substrate — the basin detectors.
- *Diagnostics, run as needed:*
  - radial × purified × SA — the gauge anchor;
  - the MOs start — representability, at the smallest windows;
  - virtual-mix — a registry check.
- *Redundant by design (optional):*
  - extra starts under free — demonstrate basin identity;
  - the single-root twin of a certified cell — reads the gauge.

**Gates (what the code enforces).**

1. the reference gate — character-forced occupations plus an energy
   gate over cycled starting guesses (`run_rohf(ref=…)`);
2. checkpoint trust — cached chkfiles reused as-is;
3. the orthonormality guard — C†SC = 1 checked before any analysis
   (the divergence guard);
4. construction gates — block-safe core-orthogonalization;
   near-linear dependence raises a named error rather than NaN-ing;
5. basin branches — where SCF chaos is bimodal, verdicts branch on
   window weights, never on noise.

**Direct runs — analyzed, not avoided.** Gating is half the design; the
other half runs the *ungated* variants (PySCF used directly:
energy-ordered MO picks, free optimization, the raw fold) through the
same analysis. The direct class is identifiable in the grid — stock MO
starts under the free mask, raw fold — and every broken row of the
manifest is one of these, classified rather than discarded.

### The pipeline

A new atom is a declared case, executed top-down — per-atom *data*
(closed/open AO spans, the reference energy E_ref, the combinatorial
term tables), not per-atom code:

```
CaseSpec — the declaration
│  atom, basis/ECP, spin, nelec_act, ncore, shell_patterns
│  ref: ReferenceSpec (occ spans l0..l3 | label prefixes,
│        conv_accelerator ∈ {"diis","newton","diis+newton"}, occ_e_ref)
│  run_casscf: the single-root stage flag; planned_windows: the
│        not-yet-run window labels (the cell grid derives from them)
│
└─ run_case(spec) — the classic staged pipeline (cases.py)
   ├─ 1. ROHF reference (reference.py)
   │      a cached checkpoint is reused as is; a fresh run applies
   │      character-forced filling, the declared acceleration, and
   │      retries guesses until the gate |E − E_ref| is satisfied.
   ├─ 2. CASCI per sector — fold ∈ {averaged, raw, fractional}
   │      term tables produced and compared to the combinatorial
   │      contract; deviations are reported.
   ├─ 3. if run_casscf → single-root CASSCF
   └─ 4. SA-CASSCF over the ground term's M_L set

The methods grid is the OTHER deliverable path (cells.py): every
treatment × start × mask × substrate × fold combination is one
addressed CellKey derived from the same declaration — the rungs-era
modules (treatments.py, study.py) retired 2026-10-29, their pins
migrated to cell addresses. The frozen scans (scans.py) are the
optimizer-free anatomy layer, run from scripts.
```

The reference gate's energy `occ_e_ref` is self-calibrated atom data:
several starting guesses are tried, their solutions classified by
character (the S-projector weight), and the passing basin's energy
recorded.

### The variety axes

Every point of the campaign's variety is **data** — an axis value in a
`CellKey` address, never a new script. The declaration (`CaseSpec` rows
in `cases.CAMPAIGN`) says which slices are interesting; the runners
execute the declared set:

| Axis | Values | Declared in | Run by |
|---|---|---|---|
| atom × basis × ECP | any `AtomConfig` (all-electron or `ecp=`) | `CAMPAIGN` rows | both runners |
| window | any shell set the basis carries (`3d4s`, …); `planned_windows` = declared, not yet run | `CaseSpec` | both runners |
| sector (M_S) | `(na, nb)`, from spin + electron parity | derived from the case | `run.py` (all sectors) |
| method | casci / casscf / sa-casscf | the grid grammar | `regen_cells.py` |
| mask | free / frozen / radial / virtmix | the grid grammar | `regen_cells.py` |
| start | ao / mo / no | the grid grammar | `regen_cells.py` |
| substrate | pol (canonical core) / pur (closed span) | the grid grammar | `regen_cells.py` |
| fold | raw / avg / frac | the grid grammar | `regen_cells.py` |
| reference | default / fractional-occupation / character-forced gated | `CaseSpec.ref` | both runners |

Two runners, two deliverables: `run.py <atoms>` runs the fixed staged
pipeline (the term tables — the stages *are* the deliverable, no method
filters); `regen_cells.py <atoms>` runs the methods grid and selects
within it by key identity:

```bash
python scripts/regen_cells.py Ti                        # one atom, full grid
python scripts/regen_cells.py Ti --method sa-casscf     # one methodology
python scripts/regen_cells.py Ti --window 3d4s --mask radial
python scripts/regen_cells.py --group diagnostics       # a named declared group
```

Filtered runs **merge** into the existing `cells.json` (untouched cells
keep their stored payloads); only the unfiltered run replaces the
store. Anything undeclared remains explorable directly —
`run_atom(cfg, patterns, nelec_act, …)` for a one-off staged look,
`run_cell(CellKey(...), mf, ...)` for a single cell: the declaration
records what is *interesting*; it does not gate the machinery.

---

## Results

All raw-shell windows and the averaged fold unless stated; every number
is this project's measurement. **Every atom lives in
`elements/<Atom>/TABLES.md`**, each file repeating this structure —
term table, multiplet energies, the manifest, the lattice, others (C
and Si mirror each other throughout, Sc and U carry the
transition-metal and f-scale stories). The worked instance is Ti — the
artifacts under `elements/Ti/` (regeneration: Files, below); its ROHF
reference −848.406489 Eh (3d²4s², the gated reference), M_S = 1,
ground term ³F. On symmetry-broken basins the term table itself
scrambles — independent evidence beyond the commutators.

### Term table — the contract each sector must meet

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,5) / {3d} | Ms=1 | 10 | ³F: 7, ³P: 3 |
| CAS(2,5) / {3d} | Ms=0 | 25 | ³F: 7, ¹D: 5, ³P: 3, ¹G: 9, ¹S: 1 |
| CAS(4,6) / {3d,4s} | Ms=1 | 120 | ³D: 15, ³F: 35, ³G: 18, ³H: 22, ³P: 15, ⁵D: 5, ⁵F: 7, ⁵P: 3 |

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

Every row: E — the run's variational total (Eh and eV); ΔE vs the ROHF
anchor (also stated in each caption); **Span (Lowdin)** — the
converged window's shell composition (S-metric weights; fractions of
the window's columns per l — the *span*, not the electron
configuration: the swap reads s→0, the escapes show their f content);
**Elec (Lowdin)** — the ground root's electron count per l (partitions
the active electrons: what the electrons *do* inside the span).
⟨L̂z²⟩ is a **gauge diagnostic**: in an exactly closed fold the ground
is a degenerate M_L multiplet and the solver returns an arbitrary
member mix (fractional values, draw-dependent); it quantizes to M²
only when a field lifts the degeneracy — the raw fold's Λ-branching
picking |M_L| = 2 on Ti. The term labels come from ⟨L²⟩/⟨S²⟩; the
M_L-resolved reading (L̂z diagonalized inside the ground eigenspace:
⟨L̂z²⟩ = 9, 4, 1, 0, 1, 4, 9 for ³F) is the multiplet demonstration,
`elements/Ti/TABLES.md`. What is exact without any canonicalization
is the sector's **L̂z spectrum** — the M_L content of the determinant
space itself (carbon's triplet sector: [+1, 0, −1]), reported with
each sector's term analysis (`report.md`): a property of the space,
not of the solver's vector — gauge-free and fold-free. The fold-shift
bookkeeping (|E(analysis fold) − E(as run)|) is not a table column;
its reading: Partition bookkeeping, Methods.

Verdicts (the classes of the
failure inventory, Motivation): **exact** — commutators and span
certificates at the 1e-13 class; **baseline** — a named non-closure
mechanism at its signature scale (the raw fold's percent-class
commutators at unchanged energy make it the polarized-core detector);
**escape / swap** — the window-identity basins, compared in totals
only; **NOT-CONVERGED** — the optimizer failed, the row waits for a
re-run.

The certified rungs sit mEh-tight around frozen CASCI (+10.5):
orbital optimization buys ≈ nothing on the as-run window. Every other
row of the manifest is a *classified* window-identity change — escape,
swap, gauge — the class the fold-shift bookkeeping flags (Partition
bookkeeping, Methods). Values, certificates: `elements/Ti/TABLES.md`;
other atoms: their element files.

| treatment | E, Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ | verdict |
|---|---|---|---|---|---|---|
| ROHF (anchor) | −848.4065 | 0 | — | — | — | — |
| CASCI, averaged | −848.3959 | +10.5 | 12.0000 | 4.2e-09 | 2.3e-08 | exact |
| CASCI, raw | −848.3959 | +10.5 | 12.0000 | 4.9e-02 | 2.6e-01 | baseline (detector) |
| CASCI, fractional | −848.3585 | +48.0 | 12.0000 | 1.1e-09 | 7.6e-09 | exact |
| CASSCF, free | −848.4377 | −31.3 | 0.0503 | 4.3e-01 | 1.1e+00 | escape |
| CASSCF, frozen window | −848.3961 | +10.4 | 12.0000 | 4.2e-09 | 2.3e-08 | baseline |
| SA-CASSCF, free | −848.4166 | −10.1 | 12.0029 | 1.1e+00 | 6.0e+00 | swap |
| SA-CASSCF, frozen, pur | −848.3948 | +11.7 | 12.0000 | 1.0e-13 | 4.0e-13 | exact |

### The window lattice — the span question, at exact symmetry

The lattice holds the method fixed (CASCI, symmetric folds) and grows
the **window** — one nominal span per row (3d, 3d4s, 3d4s4p, …). Each
row is realized through different **starts** (the construction
columns): the AO span *is* the nominal window; MO and NO spans are
different windows of the same size, compared at **equal determinant
count**. The lattice asks: can the symmetric solution reach the broken
ones' energies by growing the span instead of breaking symmetry?
(Measured: yes — both transition metals cross below their escaped
CASSCF basins; the escape ceilings are window statements, not method
limits.) The treatment grid (Methods compared) asks the complementary
question at fixed window: which optimizations break, and what restores
them.

Every row carries the span certificate ‖[P,Lz]‖ / ‖[P,L⊥]‖: the AO and
NO spans close at ~1e-13; the raw-MO spans never close (O(1)–9 at
every size measured, 4–100 orbitals — the SCF's open orbitals are
m-delocalized), so the certificate disqualifies them without running
any CI. Per-series tables (3d² and 3d²4s², construction columns,
pending slots): `elements/Ti/TABLES.md`. The AO lattice descends
monotonically below ROHF at exact ⟨L²⟩ (the d-atom window already
correlates); the MO column wins only the smallest pair — its window
contains the reference determinant (representability) — then saturates
with ⟨L²⟩ never quantizing; the NO column matches the AO value at
equal dets and exceeds it at the top. U: the same at f-scale, with a
reference determinant that is a numerically pure ⁵S° state
(`elements/U/`).

Two construction facts, measured as-run on the 3d² series. First:
**CASCI over the ROHF MO window equals the ROHF exactly** — singles
are Brillouin-dead, and doubles are dead unless same-l (the
same-center different-l exchange integral vanishes; parity kills
d¹p¹↔d² as well), so the first live channel is the same-l double
shell: +0.00 µEh on {3d} and {3d,4p}, −39.66 µEh once 4d enters,
while the AO rows are unbounded (the raw span does not contain the
reference determinant). Second, the **NO** source: the ROHF density
is idempotent — its natural orbitals *are* the MOs with occupations
exactly {0, 1, 2}, so fractional-occupancy selection from ROHF
returns the MO window; the meaningful source is the term-averaged
CASCI density on the ROHF-orthogonalized window, whose occupation
groups are l-pure m-sets plus correlating tails.

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
tables), `cells.json/md` (the full treatment cell set, key-addressed —
Methods compared), `scan.md`, plus the frozen artifacts
`sa_decomposition.md` / `basin_map.md` (their one-shot generators
retired with the rungs era; committed data, not regenerable).
New atom = new element file, not README growth; the
regenerators write straight into these files (`scripts/regen_cells.py Ti`
-> `elements/Ti/cells.{json,md}`; chkfiles land under `rungs/` and are
gitignored). Regeneration: term tables — `scripts/run.py` (the declared
campaign; the per-atom `run_atom` calls are recorded in each report
header); cells — the
declared key-set of `scripts/regen_cells.py`; the methods tables —
`scripts/build_methods_tables.py` (renders the factor grid from
`cells.json` + the audit verdicts); the cell-sourced lattice sections —
`scripts/build_lattice_tables.py`; scans —
`scripts/radial_scan.py`; the U lattices —
`scripts/run_u.py`, `run_u_methb.py`, `run_u_sa.py`. Executable documentation: the test
suite (184 tests) — term tables, exact-class closure pins, basin
branches.

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
