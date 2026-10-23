# AtomTerm

Atomic term-symbol analysis on top of [PySCF](https://pyscf.org): the pipeline
runs ROHF → CASCI → CASSCF → SA-CASSCF, then classifies every state by
**exact** angular-momentum operators L̂z, L̂², Ŝ² built as plain matrices in
the Slater-determinant basis, and certifies rotational-symmetry preservation
by ‖[Ĥ, L̂z]‖_F and ‖[Ĥ, L̂²]‖_F. Atoms covered: **C** (2p²), **Si** (3p²),
**Sc** (3d¹4s²), **Ti** (3d²4s²).

## Quick start

```bash
conda activate pyscf          # PySCF 2.13 environment
pytest                        # 91 tests, ~13 min (chkfile caching)
python scripts/run_c.py       # writes results/ (chkfiles + report.md/.json)
```

## What you get

Three result types, all written under `results/` (§8 maps every file to
its regeneration command).

**Term tables** — every state classified, per Ms sector. Carbon
(CAS(2,3), Ms=0 sector):

```
 #           E (Eh)      <L^2>  L      <S^2>  S  term
 0      -1.10060875   2.000000  1   2.000000  1    3P
 1      -1.10060875   2.000000  1   2.000000  1    3P
 2      -1.10060875   2.000000  1   2.000000  1    3P
 3      -1.01251492   6.000000  2   0.000000  0    1D      (x5)
 8      -0.88037417   0.000000  0   0.000000  0    1S
Lz eigenvalues (M_L content): [-2, -1, -1, 0, 0, 0, 1, 1, 2]
||[H, Lz]||_F  = 1.5e-15      ||[H, L^2]||_F = 3.2e-15
```

The textbook 2p² multiplet (Hund: ³P < ¹D < ¹S) comes out exactly —
commutators at round-off, because this window is exactly rotation-closed
(§2.3). §7.1–7.2; regenerate: `python scripts/run_c.py`.

**The ladder** — every orbital treatment priced against the ROHF
reference, broken and symmetric variants side by side. Titanium's
summary row set (ROHF reference −848.406489 Eh; ΔE vs ROHF, mEh):

| deliverable | free variant | symmetric variant |
|---|---|---|
| CASCI, raw shell | — | +10.6 (exact) |
| CASSCF | **−31.3** (escape — broken) | +10.4 frozen (baseline) / −2.1 radial (single root, gauge-chaotic) |
| SA-CASSCF | −10.1 (swap — broken) | +1.0 radial purified (floor) / +11.7 frozen purified (exact) |

§7.3; regenerate: the `run_ladder` calls recorded in each ladder header.

**Scan curves** — the deterministic energy landscape: freeze everything,
deform one window coordinate, re-diagonalize. §7.4; regenerate:
`python scripts/radial_scan.py`.

---

## Introduction

The exact states of an isolated atom carry sharp term symbols (³P, ¹D, …):
the Hamiltonian commutes with all rotations — the SO(3) × SU(2)
invariance of coordinates and spins — so every state has definite L and
S, and the states of one term are exactly degenerate. CASCI and CASSCF
treatments *can* lose this structure. What breaks is not quantum
mechanics but the *restriction*: the calculation is confined to a small
set of orbitals, and two properties of that restriction can fail.

First, CASCI inherits an anisotropic field from its reference (§2.4).
Carbon: the two ³P electrons occupy two of the three 2p orbitals — the
choice is arbitrary, all pairs are degenerate — say px and py; their
charge is then unevenly distributed (heavier along x and y than z), the
closed-shell core polarizes in response, and the resulting field acts
unevenly on the valence states, splitting the ³P term. The pz-pair case
is the same story in another orientation; no pair is spherical.

Second, CASSCF additionally lets the active orbitals themselves relax
(§3.3, §2.4). The active orbital set — the *window* (§2.3) — can then
deform in three ways with different consequences: some deformations
leave the energy unchanged, and the optimizer stops at an arbitrary
point of them; some *lower* the energy by tilting orbitals (on d-shells;
the Jahn–Teller mechanism, §2.4); and some change *which* orbitals are
active at all — trading a valence 4s for a core 3p (§5).

State-averaged CASSCF exists because a single optimized state represents
just one component of the multiplet; the average over all components is
the rotationally symmetric object (§3.4). It restores sharp L values,
but it cannot alone guarantee the full symmetry, and on heavier atoms it
can still trade orbitals.

Heavier atoms add failure modes (§5): the energy-lowering tilts, the
orbital trades, and the reference calculation itself can converge to a
wrong or metastable solution (§5).

Is the minimum correct? Where the energy is insensitive to a
deformation, the energy is still well defined — the orbitals are not,
and that is fine (§2.4). Across genuinely different solutions, every
free run must be inspected for which orbitals it actually uses before it
is believed (§7.3).

---

*How to read.* §1 states the goal. §2 fixes the objects, the language
and the diagnostics — later sections use nothing that is not defined
there. §3 walks the method hierarchy (why each exists, what it computes,
how we run it, and what textbook usage gets wrong). §4 maps the methods
onto the code. §5 is the core: where symmetry breaks, level by level.
§6 is the operational protocol card. §7 holds the results. §8 is the
file map; §9 the roadmap. Plain-language companions: `RUNGS.md` and the
glossary appendix below; agent instructions: `AGENTS.md`.

---

## 1. Goal and scope

**The goal: atom term assignments from CAS treatments, with certified
L̂² and L̂z symmetry, for any atom.** "Certified" means the three
questions of §2.1 answered against the classes of §2.5 — sharp term
labels, exact-class commutators, exact term tables — not approximate
agreement. The treatment delivered is CASCI/CASSCF-class,
non-relativistic throughout this note; the roadmap (§9) extends to
scalar-relativistic Hamiltonians (which change h₁ only — a scalar is
rotationally invariant, so the entire machinery survives) and to
ECP-based full relativity. The evidence base is four atoms (C, Si, Sc,
Ti; cc-pVDZ); the machinery is atom-generic, with per-atom *data*
(closed/open AO spans, the reference energy E_ref, the combinatorial
term tables) rather than per-atom code.

---

## 2. The object and the diagnostics

### 2.1 The symmetry an atom demands — and the question of this project

An isolated atom is invariant under rotations. Its Hamiltonian therefore
commutes with L̂², L̂z, Ŝ², Ŝz, and every exact eigenstate is
|γ; L, M_L, S, M_S⟩ — it carries a **term symbol** (2S+1)L. Because of
L̂² and Ŝ² symmetry, the (2L+1)(2S+1) components of one term are *exactly
degenerate*: ³P is one energy, nine times.

The goal of this project is to understand **why practical CAS stages
(CASCI over an SCF reference, CASSCF, SA-CASSCF) lose this symmetry**, and
how to follow it after all — without giving up the correlation these
methods buy. A treatment is examined through three questions:

1. does every diagonalized root carry a sharp term label
   (⟨L̂²⟩ = L(L+1), ⟨Ŝ²⟩ = S(S+1))?
2. how large are the commutators ‖[H, L̂z]‖, ‖[H, L̂²]‖ — evaluated as
   exact matrices over the determinant basis — compared against the
   numerical classes of §2.5?
3. does the term table (state counts per term) equal the combinatorial
   free-ion table?

### 2.2 The operators, over determinants

We work with spin-orbitals (α = 0..n−1, β = n..2n−1) and
**determinants** — sorted occupation tuples. A **sector** fixes the
numbers (n_α, n_β) of α/β electrons, hence M_S = (n_α−n_β)/2 and N. The
operators N̂, L̂z, L̂², Ŝz, Ŝ² all preserve sectors (each commutes with Ŝz
and N̂), so every matrix lives inside one sector of fixed M_S.

A CAS treats a chosen set of active orbitals — the **window** (analyzed
in §2.3). Every operator is brought in as **one-body symbol + two-body
symbol**: the symbol is the operator's matrix in the active-orbital
basis, built from its AO-basis matrix via the orbital coefficients; both
implementations below consume the same symbols. Evaluation is by
Slater–Condon rules over determinants (`detbasis`), or — elementwise
identically, pinned by test — by explicit creation/annihilation on
occupation vectors (`secondquant`). The two implementations validate each
other.

| Operator | One-body symbol | Two-body symbol |
|---|---|---|
| N̂ | identity | — |
| L̂z | Lz matrix in the active-orbital basis | — |
| L̂² | Σᵢ L̂ᵢ² (l² AO matrices, transformed) | 2 Σᵢ<ⱼ L̂ᵢ·L̂ⱼ (composite integrals) |
| Ŝ² | — | closed combinatorial form (pair exchange) |
| Ĥ | folded h₁ (nucleus + core field) | active-orbital two-electron integrals |

### 2.3 The closure theorem

**We call a space rotation-closed if L̂W ⊆ W, and an operator a scalar if
it commutes with all rotations.** The members of an l-shell are its 2l+1
m-functions (px, py, pz; dxy, dyz, dz², dxz, dx²−y² — where we refer to
member *sets*, the indices 0–4 follow this order).

**Exact symmetry imposes two conditions on the treatment.**

1. **The window is rotation-closed.** Otherwise L̂ maps window functions
   outside the window, and the CAS Hamiltonian inherits the leak:
   [H, L̂] ≠ 0 even though the physical atom is symmetric. Violations
   arise in practice because canonical-MO windows and optimizer output
   carry member-dependent radial mixtures, because truncation-by-weight
   can drop members, and because even core-orthogonalization of a raw
   shell window can inject non-closed directions (the mechanisms, with
   numbers, are §5's subject).
2. **The core fold is a scalar.** The inactive core enters the active
   problem as a one-body field inside h₁; if that field is angularly
   anisotropic, it breaks L̂z directly. This happens because the SCF core
   polarizes in the anisotropic field of the open shell (§2.4).

**Rotation generators never touch radial coordinates.**

> **Theorem (closure).** For φ(r)·Y_lm one has L̂(φ·Y) = φ·(L̂Y) — the
> radial factor is untouched. Therefore a finite one-particle space W is
> rotation-closed **iff**, for each l present, W contains the complete
> set of 2l+1 members, all sharing the *same* radial subspace E:
> W = ⊕ₗ (E_l ⊗ {all members of l}). The radial subspace E itself is
> arbitrary.

**Exactly two failure modes exist: members missing, and radials not
shared across members.**

**Window operators never leave the window Fock space.** Constructing
"L̂²" from a†_p a_q (p, q ∈ W) therefore gives the *projected* L̂²; it
equals the true L̂² on W's determinants iff W is rotation-closed.

**Two projector norms check the window alone.** With no CI and no
optimizer involved, ‖[P, L̂z]‖ and ‖[P, L̂⊥]‖ — P the orthogonal
projector onto the window span in the AO overlap metric — catch exactly
the two failure modes above.

**A member tilt is the elementary window deformation.** One member is
replaced by a radial admixture, cos θ·(n l m) + sin θ·(n′ l m) — same
l, same m, neighboring radial. In the real AO basis L̂z mixes only the
two functions of each ±m pair and annihilates m = 0 members, so
‖[P, L̂z]‖ has **blind spots** (an m = 0 tilt, or equal tilts of both
partners of a pair, are invisible to it); L̂ᵔ rotates any member into
any other, so ‖[P, L̂ᵔ]‖ sees every member asymmetry — together, no
blind spots. Both directions of the asymmetry are measured (per-atom
records: `elements/*/NOTES.md`).

**A flat direction is a deformation along which the energy is constant
to numerical precision.** The optimizer's endpoint on one is decided by
numerical noise alone — also called a *gauge choice*: energetically
indistinguishable alternatives, distinguishable only through the
closure diagnostics (the working census: `elements/C/NOTES.md`).

**Sector content decides what its commutators can testify.** On a
sector carrying a single irreducible L-representation, L̂² =
L(L+1)·**1** identically, so ‖[H, L̂²]‖ vanishes *whatever is wrong* —
multi-term sectors are the honest courts. The projector pair, being
sector-independent, is the universal certificate.

### 2.4 The diagnostics and their sensitivity

**ε measures the size of the symmetry breaking.** Throughout, it is the
weight of the angularly non-scalar part of the core field's density (a
scalar core field contributes nothing to any commutator), reported as
core anisotropy; ‖A‖_F = sqrt(Σ|A_ij|²) is the Frobenius norm.

| Diagnostic | What it measures | Sensitive to | Blind to |
|---|---|---|---|
| ‖[H, L̂z]‖_F, ‖[H, L̂²]‖_F | the assembled treatment (window + core fold) | everything, order ε | vacuous in single-irrep sectors (L̂²); ±pair structure (L̂z) |
| ‖[P, L̂z]‖, ‖[P, L̂⊥]‖ | the window span alone (members §2.3) | member asymmetry, missing members; no optimizer noise | — (as a pair) |
| ⟨L̂²⟩, ⟨Ŝ²⟩ per root | term identity | order ε² | small defects |
| multiplet splitting | fold anisotropy | order ε | — |
| orthonormality error | validity of the stage's MO matrix (C†SC = 1) — a guard, not a symmetry diagnostic | construction bugs | — |

**All five diagnostics respond to one object**: H = H₀ + εV, with
ε the *measured* size of the breaking (core anisotropy, 1e-3…1e-2
relative across the atoms) and V its shape. Commutators and splittings
track ε to first order; expectation values shift at ε². The
bookkeeping is verified across the atoms.

**The breaking must act twice to move an expectation value, but only once
to split a degenerate set.** An expectation value like ⟨L̂²⟩ changes only
because the breaking mixes a term with other L-terms, with amplitude of
order ε; the eigenstate is contaminated at order ε and its expectation
shifts by the squared amplitude. Within one degenerate multiplet there is
no gap to cross: the anisotropic field acts directly on the degenerate
set and separates its components at order ε, in close analogy with Stark
splitting. We call these component energy differences the **splitting**
and their mean the **centroid**; the centroid, averaged over the
multiplet, again moves only at order ε².

**Core anisotropy has a physical origin and three cures.** The open
shell's quadrupole polarizes the core during the SCF (§2.3, condition
2); the cures — spherical averaging of the core density, the
fractional-occupation reference, core purification — are compared by
cost in §7.

**Orthonormality is a validity guard, not a symmetry diagnostic.** If
the stage's MO matrix is not orthonormal in the AO overlap metric the
construction is invalid (a non-orthonormal input once made an optimizer
return sub-variational energies; `elements/Si/NOTES.md`). Every stage
reports its orthonormality error.

**The chain from quadrupole to splitting closes quantitatively.**
Inversion symmetry forbids a dipole, so the leading anisotropic term is
the quadrupole — whose signature (the {M_L = 0} vs {M_L = ±1} pattern)
identified it. The measured ε is the induced quadrupole weight; the
splitting follows at order ε with an atom-dependent coupling:
| Atom | ε (core anisotropy, relative) | splitting (Eh) | splitting/ε (Eh) | ⟨L̂²⟩ of the ground term (raw fold) |
|---|---|---|---|---|
| C | 2.405e-3 | 4.373e-5 (³P) | 0.02 | 2.0000000000 |
| Si | 2.109e-3 | 2.431e-4 (³P) | 0.12 | 2.0000000000 |
| Ti | 7.754e-3 | 4.144e-3 (³F septet) | 0.53 | 12.0000000000 |
| Sc | 1.259e-2 | 1.857e-3 (²D quintet) | 0.15 | 6.0001 |

The ratio splitting/ε is the effective coupling strength — an atomic
property. The order relations hold throughout: commutators and
splittings are the detectors; a numerically pretty ⟨L̂²⟩ is not.
### 2.5 Numerical conventions — what “exact” costs

| Class | Typical scale | Meaning | Always resolvable? |
|---|---|---|---|
| exact | ~1e-13 | round-off of the construction | — |
| baseline | 1e-8…1e-2 | a *named* mechanism: core-orthogonalization defect, core polarization, optimizer floor | each has a signature (checkable); all resolvable except the optimizer floor — certifiable and boundable, but with no restoring force to remove it (§5) |
| broken | 1e-1…10 | window identity change (swap/escape, §5) | yes — recognized and rejected by the diagnostics |

§5 maps every rung to a class.

**The optimizer's deformations (derived).** Over m-labeled determinants
H is block-diagonal in M_L: the one-body fold is m-diagonal, and the
two-body Coulomb term conserves total M_L (Δm₁ + Δm₂ = 0). Three
consequences:

- *Tilt flatness.* Where the ground multiplet's M_L blocks are
  one-dimensional, H is diagonal there and each determinant's energy
  depends only on its own members. Tilting one member (its angle
  δ the parameter) leaves the determinants avoiding that member
  untouched; the blocks are degenerate, so the ground energy is
  constant in δ — the optimizer may break m at zero cost, and its
  endpoint is noise.
- *Tilt descent (Jahn–Teller-like).* In multi-dimensional M_L blocks a
  tilt also moves the couplings. On the degenerate multiplet the
  deformation's block is traceless — its trace is the derivative of
  the equal-weight ensemble energy, zero at a stationary point — so
  the eigenvalues split symmetrically about zero and the lowest
  descends linearly: E₀(δ) = E₀ − |V(δ)| + O(δ²). A degenerate
  state, a symmetry-lowering coordinate, a guaranteed descent.
- *The state average.* Because the block is traceless, Σₐ ∂Eₐ/∂δ = 0:
  the ensemble is stationary at the m-symmetric window — stationary
  with zero curvature, however, which is not a minimum with restoring
  force (§3.4).


---

## 3. The method hierarchy

Each method gets the same triple: why it exists, what it computes, how we
run it. Two vocabulary items are defined here and used everywhere after:
the *substrate* (§3.3) and the *basin* (§3.3).

### 3.1 ROHF — the reference that defines "core"

**The reference calculation fixes what the core is: everything downstream
folds it.** Two of its aspects matter here: how an orbital's *type* is measured,
and how the reference is kept honest where energy-ordered filling
fails.
**An orbital's type is measured by projection.** For an AO span with
index set o, the S-metric projector is P_o = S[:,o] S_oo⁻¹ S[o,:], and an
orbital c carries weight w = c†P_oc ∈ [0,1] on that span (bounded — a
true projector weight; the same metric appears in §3.2's window
construction). Occupations can then be assigned by *character* rather
than energy order.

**Energy-ordered filling can converge to the wrong answer — the
well-known improper-convergence problem, here at SCF level.** Two
measured faces: it fills the wrong *configuration* (a doubly-occupied d
inside the core, ~1e2 mEh above the physical one; every energy-ordered
start tested lands there), and with occupations fixed the degenerate
open shell still leaves a *metastable basin* (~10 mEh up, same
occupations) whose polarized core degrades everything downstream. The
toolbox, in order: **character-forced filling** (occupy each cycle by
projector weight); the **S-projector character metric** for
near-degenerate bases; an **energy gate** with cycled starting guesses
until the reference basin is reached; **cached trusted references**
(E_ref and the AO spans are the atom's data). Records:
`elements/*/NOTES.md`.
**The reference is also the source of the core anisotropy ε.** The open
shell's quadrupole polarizes the core during the SCF (§2.4); the fold
variants of §3.2 are three ways of dealing with what this calculation
leaves behind.

### 3.2 CASCI — diagonalization at fixed core and window

**CASCI diagonalizes the CI problem inside one sector, at fixed core and
window.** The window is the raw shell AOs (rotation-closed by
construction, §2.3); the core enters through the one-body fold, in three
constructions that differ exactly in how they treat ε:

- **averaged** — one-shot post-processing of the *same* ROHF core: the
  core density is spherically averaged (its l = 0 component kept),
  ρ → ρ̄, before the fold field is built. The anisotropic field is
  discarded; the core energy and bookkeeping stay those of the ROHF
  reference.
- **fractional-occupation** — a *different SCF run*: the open-shell
  electrons are smeared uniformly over the shell (occupation
  n_open/(2l+1) per AO), and the whole reference is re-converged in this
  ensemble. The anisotropy never develops self-consistently; the reference is a
different functional of the smearing, not a variational neighbor of
ROHF.
- **raw** — the core field as the SCF left it: ε unscreened, the
  multiplet splits by §2.4's mechanism.

With the averaged fold, the CASCI term tables of §7.1–7.2 come out exact
on all four atoms; the averaged–raw and averaged–frac differences are
the ΔE anchors of §7.3.

**The fold, explicitly.** With core orbitals {φᵢ} doubly occupied,
restriction to active-space wavefunctions gives

  H_act = Σ_pq h_eff(p,q) a⁺_p a_q + ½ Σ_pqrs (pq|rs) a⁺_p a⁺_q a_s a_r + E_core,

h_eff(p,q) = h(p,q) + Σᵢ [2(pq|ii) − (pi|iq)] — each core electron's
Coulomb–exchange field folded into the one-body part; E_core is the
core's own energy. That construction is *the fold*. Our fold vs
PySCF's raw fold = spherically averaged vs raw core density inside
those J/K terms.


### 3.3 CASSCF — letting the orbitals move

**CASSCF minimizes the chosen root's energy over orbital rotations.** Its
motivation is visible in carbon's radial gain — the valence 2p functions
want a 3p admixture worth −250.8 mEh of window energy (§7.4 — in totals
this is the recovery of the raw-shell confinement cost, §7.3), a
relaxation CASCI cannot
express. How we run it is a grid of *starts × restrictions*, and both
axes matter.

**The substrate is the starting MO matrix the optimizer receives** —
core block, window block, complement. Ours is block-safe, and the window
step is explicit mathematics: a raw window column w is replaced by

  w̃ = (1 − P_core) w / ‖(1 − P_core) w‖_S,  P_core = C (C†SC)⁻¹ C†,

i.e. its component in the core span is subtracted in the AO overlap
metric and the remainder re-normalized. The construction is forced by convergence — a non-orthonormal
substrate makes the optimizer diverge to sub-variational energies
(observed once; `elements/Si/NOTES.md`, the §2.4 guard).


**Restrictions are masks on rotation pairs.** *Window frozen*: only
core↔virtual pairs — correlation through the core alone, closure
untouched. *Radial restriction*: only pairs within a same-(l, m) radial
chain — the symmetry-preserving channel; where it exists (C, Si) it
recovers the full gain at exact closure. *Free*: no mask — the textbook
run. **A basin is the region of attraction of one local minimum of this
optimization**; the free landscape has several (§7.3's chaos paragraph).


**The rotations themselves.** The optimizer parameterizes C(κ) =
C₀·exp(κ) with κ antisymmetric in the S-metric — one variable per
allowed rotation *pair*; the masks above select the pairs. A single
pair (active a, virtual v) is the two-dimensional rotation

  c_a(θ) = cos θ·c_a + sin θ·c_v,   c_v(θ) = −sin θ·c_a + cos θ·c_v.

The orbital-content changes we report are projection bookkeeping:
w_X(k) = c_k[X]† S[X,X] c_k[X], the squared S-projection of column k
onto an AO row-block X (an exponent band, or all rows of one l). Same
total l-weight with redistributed bands — *radial relaxation*: the
orbitals re-project onto the AO lattice where the energy wants them.

### 3.4 SA-CASSCF — averaging over the multiplet

**State-averaging replaces the single root by the multiplet average:**

  E_SA(orbitals) = (1/(2L+1)) Σ_a E_a(orbitals),  a over the M_L
  components,

minimized over rotations (³P: three components; ²D: five; ³F: seven).
Where a single root represents one M_L component — and optimizes for it
alone — the average is the rotationally scalar object.

**The average anchors the term identity but not the window.** ⟨L̂²⟩
is restored, yet along flat directions the average is stationary with
no restoring force (§2.4: the deformation block is traceless), so
convergence tolerance, not physics, decides the endpoint; and other
basins remain reachable (the swap, §5). Exactness combines the tools:
state average + frozen window + isotropic substrate. Whether tilting
helps or hurts the average is atom-dependent — the frozen-scan
verdicts of §7.4 (per-atom values in `elements/`).


### 3.5 What textbook runs give you, and why

| Standard practice | What you get | Why | What we do instead |
|---|---|---|---|
| CASCI over canonical virtual MOs | closure defects at the 1e-1 level; term labels unreliable | canonical MOs diagonalize a Fock built from an anisotropic density; quasi-degenerate same-l members mix arbitrarily | raw shell windows, averaged fold (§3.2) |
| free CASSCF, default start | p²: full gain with a silently member-tilted window; d-shells: escape/swap basins, ⟨L̂²⟩ off by units | flat directions (endpoint = noise) and descending directions (§2.4); window identity free to change | restrictions + substrate control (§3.3) |
| free SA-CASSCF | ⟨L̂²⟩ ≈ L(L+1) restored, closure not (Sc: 6.002 with ‖[H,L̂z]‖ ≈ 0.8) | the average is scalar to first order only; no restoring force on flat directions; other basins still reachable | frozen-window SA on the isotropic substrate (§5) |

## 4. The pipeline tour

| Module | Owns |
|---|---|
| `system` | system validation — one atom at the origin, spherical basis (L̂ defined about r = 0) |
| `angular` | one-electron L̂x/L̂y/L̂z AO matrices — two independent constructions pinning each other (libcint; exact tesseral algebra) |
| `detbasis` / `secondquant` | many-body matrices over determinants — Slater–Condon and explicit a†/a, elementwise-identical twins (§2.2) |
| `operators` | operator assembly — the fold (§3.2), the H/L̂z/L̂²/Ŝ²/N̂ symbols, the window-closure certificates (§2.3) |
| `reference` / `casscf` | the engine — SCF reference with gates; the one deep optimizer (`constrained_casscf`) with its masks and adapters |
| `windows` / `core` | window construction, selection, the kind registry (`Window.columns`); the purified core substrate and fold energy |
| `stages` | the one uniform analysis pass per finished stage — the source of every §7 number |
| `scans` | frozen-landscape measurements — deterministic (§7.4) |
| `analysis` | term labels, configuration summaries, report rendering |
| `cases` / `ladder` / `sweep` | orchestration — the declared per-atom campaign (`cases`: CaseSpec → run_case); rungs as data (`ladder`, §7.3); the window × treatment runner (`sweep`) |

**The drivers capability map.** Every treatment is a point in a
five-axis space:

| Capability | What it is |
|---|---|
| ROHF reference | character-forced filling + energy gate — the anchor (§3.1, §5) |
| CASCI | exact sector diagonalization — fold variants raw / averaged / fractional (§3.2) |
| CASSCF | orbital optimization — freedom by mask, start by substrate (§3.3) |
| SA-CASSCF | multiplet-averaged optimization — weights over the M_L components (§3.4) |
| core construction | the purified, exactly isotropic core span (§5) |

The 17-rung registry tabulated in `results/` enumerates the corners of
this space; §7.3 keeps the rows that differ in mechanism.

---

## 5. Where symmetry breaks

| Level | Symptom (order of magnitude) | Mechanism | Fix | Cost |
|---|---|---|---|---|
| 1. Basis bookkeeping | labels mislead; shells fuse | wrong inputs → wrong constructions | label + metric discipline | none |
| 2. Sector witnesses | ‖[H,L̂²]‖ = 0 with defects present | single-irrep sectors make L̂² = L(L+1)·**1**; ±pair blind spots (§2.3) | projector pair; multi-term sectors | sector size |
| 3. Core field | splittings and commutators 1e-5…1e-1 | polarization chain (§2.4) | fold cures / purification | ΔE per §7 |
| 4. Orbital optimization | tilted windows, floors, escape/swap | flat and descending directions (§2.4); identity change | restrictions + substrate | 0…~30 mEh above the deepest minimum |
| 5. SCF reference | ~1e2 mEh wrong configuration; ~10 mEh metastable basin | energy-ordered filling; degenerate-open-shell basins | toolbox of §3.1 | per-atom reference data |

Class mapping (§2.5): *exact* — raw-window CASCI, frozen-window rungs
on the isotropic substrate; *baseline* — the start-window defect and
the flat-direction floor; *broken* — the escape and the swap (L4).

**L1 — basis bookkeeping.** AO labels are truthful yet treacherous:
pattern matching is regex-flavored; shells fuse (one label, several
contracted functions — unpack through contraction data); same-l shells
of different n overlap (~0.3), so every projection and normalization
runs in the AO overlap metric, never the Euclidean one.

**L2 — sector witnesses.** A single-irrep sector makes ‖[H, L̂²]‖ vacuous
and L̂z has ±pair blind spots (§2.3); the projector pair is the
universal certificate, multi-term sectors the honest courts.

**L3 — the core field.** The reference leaves an anisotropic field
acting at first order on the multiplet (§2.4). Two cures act at the
density level — spherical averaging, the fractional-occupation
reference. The third, core purification, exists because of the
*window*, not the core: the start window is (1−P_core)-projected
(§3.3), so its closure inherits the core's; rebuilding the core span
as E_l ⊑ {all members of l} (§2.3's theorem applied to the core)
zeroes that factor and the fold becomes scalar by construction.

**L4 — orbital optimization.** Three mechanisms. **(a) Start-window
defect**: the orthogonalization subtracts non-closed directions when
the core span is not closed, ‖[P_w̃, L̂]‖ ≲ ⟨window|core⟩ ·
‖[P_core, L̂]‖ — percent-class where the window–core overlap is large;
the purified substrate zeroes the second factor. **(b) Flat
directions**: no restoring force; the floor is certifiable but not
removable by convergence tightening — remove the freedom instead
(radial restriction; frozen window). **(c) Identity basins**, measured
by the span overlap σ = σ_min(C_a†SC_b): the **escape** rebuilds the
window from orbitals the configuration never implied — its states
describe no term; the **swap** trades the window's identity (a valence
4s for a core 3p) — variational and approximately term-readable, but
a different core/active partition: compare in totals only, and wanting
that correlation honestly means enlarging the window deliberately. A
rebuilt window is not broken — an m-incomplete one is; closure is
imposed or certified, never converged to.

**L5 — the SCF reference.** The §3.1 toolbox; E_ref is extracted once
per atom by the diagnostics themselves (occupations by projection,
core purity, window closure) — self-calibrating, per-atom data, not
code.

---

## 6. The protocol card

**Run all five deliverables, report all five** — the certified treatment
and the direct values are complements, not alternatives. Each deliverable
is a named configuration in the §4 capability space; each optimizer
deliverable is reported as the pair *free variant / symmetric variant*.

| Deliverable | Configuration | Purpose | Acceptance |
|---|---|---|---|
| ROHF reference | forced filling + gate | the anchor every ΔE hangs on | occupations & core purity by projection; E = E_ref |
| CASCI anchor | raw shell window, averaged fold | the exact term tables | tables = combinatorial; commutators exact |
| Certified treatment | SA-CASSCF + window frozen + purified core | exact-class symmetry with correlation | commutators & window spans exact; tables intact |
| Direct CASSCF | free, shell substrate | the energy-seeking value + basin anatomy | classified (span overlap, window weights); totals quoted |
| Direct SA-CASSCF | free SA-CASSCF, shell substrate | best direct value; term-readable even in the swap basin | term label; basin classified; totals quoted |

For a new atom: (1) run the reference (spans and E_ref are the atom's
data, §3.1/§5); (2) the CASCI anchor — term tables must come out exact
(§7.1–7.2); (3) the certified treatment (expect the exact class, §7.3's
pur rows); (4) both direct runs, classified; (5) optionally, frozen
scans for the radial/tilt anatomy (§7.4).

---

## 7. Results tables

Five table kinds, each answering one question. **T1 term content**
(§7.1) — the combinatorial contract the exact stage must return;
**T2 spectrum** (§7.2) — the sharp CAS-level term energies; **T3
treatment × window** (§7.3) — what each orbital treatment buys and
costs, per window; **T5 frozen scans** (§7.4) — where the gain lives,
optimizer-free; **the window lattice** (§7.5) — window kinds at matched
determinant cost. Per-atom records live in `elements/<Atom>/TABLES.md`
(+ `NOTES.md` for the dated specifics); **Ti is the worked example**
below. All numbers are this project's measurements (C, Si, Sc, Ti on
cc-pVDZ; U on ECP + cut-Dyall); the artifacts live under `results/`
(regeneration in §8). On symmetry-broken basins the term table itself
scrambles — independent evidence beyond the commutators.

### 7.1 Term tables — what each sector must contain

Every CAS below uses raw shell AOs as the window and the spherically
averaged core fold; the tables are the *combinatorial* content — the
exact stage must reproduce them state-for-state.

| Atom, CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| C — CAS(2,3) / {2p} | Ms=1 | 3 | ³P: 3 |
| | Ms=0 | 9 | ³P: 3, ¹D: 5, ¹S: 1 |
| Si — CAS(2,3) / {3p} | Ms=1 / Ms=0 | 3 / 9 | as carbon |
| Sc — CAS(1,5) / {3d} | Ms=½ | 5 | ²D: 5 |
| Sc — CAS(3,6) / {3d,4s} | Ms=½ | 90 | ²D: 20, ⁴F: 14, ²F: 14, ²G: 18, ⁴P: 6, ²P: 6, ²S: 1, ²H: 11 |
| Ti — CAS(2,5) / {3d} | Ms=1 | 10 | ³F: 7, ³P: 3 |
| | Ms=0 | 25 | ³F: 7, ¹D: 5, ³P: 3, ¹G: 9, ¹S: 1 |
| Ti — CAS(4,6) / {3d,4s} | Ms=1 | 120 | ³D: 15, ³F: 35, ³G: 18, ³H: 22, ³P: 15, ⁵D: 5, ⁵F: 7, ⁵P: 3 |

Sc CAS(1,5) is the Sc⁺-like d¹ space on the [Ar]+4s² closed core — the
analog of Ti's CAS(2,5).

### 7.2 Multiplet energies — the sharp, degenerate spectra

Active-space energies (our averaged fold), exact within each term to
≤1e-9; the M_L column comes from the sector Lz spectra. Silicon is as
carbon throughout. Orderings match the free-ion tables (NIST); the gaps
are CAS-level values, not experimental ones.

| Atom | Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|---|
| C | ³P | 3 | −1.100609 | — | −1…1 |
| | ¹D | 5 | −1.012515 | +0.0881 | −2…2 |
| | ¹S | 1 | −0.880374 | +0.2202 | 0 |
| Sc (CAS(1,5)) | ²D | 5 | −0.342549 | — | −2…2 |
| Ti (CAS(2,5)) | ³F | 7 | −1.468663 | — | −3…3 |
| | ¹D | 5 | −1.423747 | +0.0449 | −2…2 |
| | ³P | 3 | −1.413969 | +0.0547 | −1…1 |
| | ¹G | 9 | −1.398223 | +0.0704 | −4…4 |
| | ¹S | 1 | −1.295050 | +0.1736 | 0 |

### 7.3 T3 — the treatment ladders, one row per mechanism

The full per-atom registry (every treatment, start-variant and
diagnostic) lives in `results/*_ladder.md` and
`elements/<Atom>/TABLES.md`. Here: the cross-atom summary, then the
per-atom mechanism tables with full diagnostics.

**How to read.** Every table carries ONE energy column — **E (as
run)**, the run's own variational total, comparable across rungs; ΔE
is vs the ROHF anchor. The **cross-check** column = |E(our analysis
fold) − E(as run)|, max over roots: our analysis re-evaluates every
rung in the spherically averaged core field (the common ε-free ruler
behind ⟨L²⟩/commutators/term tables); ticket 25 decomposed it exactly
into the mEh-class averaging effect plus the window-definition offset
where the valence shell overlaps the core (the 4s|3s-type overlap
carries ~85–95 mEh on Sc/Ti/Si; exactly 0 on C and the pur rows) —
the fold CONVENTION itself is machine-exact (1e-15 vs PySCF's
get_h1eff). Rows that moved the core/active boundary (the swap) are
compared in totals only

**The summary cards.**

**The reference table** (totals in Eh — the anchor every ΔE hangs on;
the averaged fold leaves the ROHF total unchanged to sub-mEh; the
frac-occ SCF is a different functional):


|---|---|---|---|---|
| ROHF (raw core = averaged fold) | −37.682418 | −288.846437 | −759.736255 | −848.406489 |
| frac-occ SCF (ΔE vs ROHF, mEh) | −37.678694 (+3.7) | −288.867793 (−21.4) | −758.014967 (+1721) | −847.296157 (+1110) |

**The summary table** — every deliverable's variants, ΔE vs ROHF in mEh,
all four atoms side by side (per-runk diagnostics: the mechanism tables
below and `results/*_ladder.md`):

| rung | C | Si | Sc | Ti |
|---|---|---|---|---|
| CASCI, averaged fold | +250.8 (exact) | +221.5 (exact) | +6.7 (exact) | +10.6 (exact) |
| CASCI, raw core fold | +250.8 (comm 9.8e-6) | +221.5 (comm 5.6e-4) | +6.7 (comm 2.8e-2) | +10.6 (comm 4.9e-2) |
| CASCI, frac reference | +242.2 | +217.8 | +36.8 | +48.0 |
| CASSCF, free | +0.05 (member-tilted) | +0.1 (lands symmetric) | **−31.5** (escape) | **−31.3** (escape) |
| CASSCF, window frozen | +242.0 | +217.5 | +6.6 | +10.4 |
| CASSCF, radial restriction | +0.05 (exact) | +0.1 (exact) | −6.6 (single root, chaotic) | −2.1 (single root, chaotic) |
| SA-CASSCF, free | +0.03 (lands symmetric) | +0.04 | −13.3 (swap) | −10.1 (swap) |
| SA-CASSCF, frozen, purified core | +242.0 | +217.7 | +7.3 (exact) | +11.7 (exact) |
| SA-CASSCF, radial, purified core | +0.03 (exact) | +0.04 (exact) | +0.0 (floor) | +1.0 (floor) |

Legend: **exact** = commutators and span certificates at the 1e-13
class; **escape** = the free rung rebuilt the window (describes no
term); **swap** = the window's identity traded (4s ↔ core 3p — a
different partition, totals only); **floor** = a flat-direction
endpoint, certifiable but without restoring force; **member-tilted** =
a flat-direction gauge choice. On p² the CASCI rows are the
confinement cost; on the d atoms the as-run window already correlates;
the raw-fold row is the polarized-core detector (same energy,
commutator jumps). The full story of each atom: its element file.

**Sc — CAS(3,6), mechanism table** (E_ground in our fold; ΔE in
totals vs ROHF, mEh; diagnostics per row):

| rung | E_ground | ΔE | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ | ‖[P,Lz]‖ | ‖[P,Lᵔ]‖ | note |
|---|---|---|---|---|---|---|---|---|
| rohf | −759.736255 | 0 | — | — | — | — | — | 1 |
| frac-scf | −758.014967 | +1721 | — | — | — | — | — | 2 |
| casci_avg | −759.817560 | +6.7 | 6.0 | 5.6e-15 | 6.4e-13 | 7e-15 | 9e-15 | 3 |
| casci_raw | −759.819083 | +6.7 | 6.0004 | 2.8e-02 | 1.4e-01 | 7e-15 | 9e-15 | 4 |
| casci_frac | −759.787295 | +36.8 | 6.0 | 3.4e-14 | 3.5e-13 | 7e-15 | 9e-15 | 2 |
| CASSCF (free) † | −759.765900 | **−31.5** | 0.05 | 3.2e-01 | 7.2e-01 | 2.3 | 2.7 | 5 |
| casscf_relax | −759.728374 | +6.6 | 5.9998 | 2.8e-03 | 4.0e-02 | 5.7e-03 | 2.1e-02 | 6 |
| sa † | −760.548082 | −13.3 | 6.002 | 8.3e-01 | 3.8 | 1.4 | 1.4 | 7 |
| sa_radial | −759.735593 | −0.8 | 5.9997 | 3.3e-03 | 5.4e-02 | 6.0e-03 | 2.2e-02 | 6 |
| pur family | −759.7284…−759.7362 | +6.7 / −5.8 / 0.0 | 6.0 | ≤4.8e-13 | ≤4.8e-13 | 7e-15 | 9e-15 | 8 |

1 anchor (§3.1); raw core anisotropy 1.3e-2.
2 the smeared-shell ensemble — a different functional (§3.2); only its
isotropic core field is used.
3 spherically averaged core — the reference CASCI.
4 raw field — percent-level quadrupole (§2.4); same energy, commutators
jump (the polarized-core detector).
5 escape to the L≈0 basin — deepest minimum, identity-destroying (§5).
6 window frozen, core relaxes; carries the start-window defect (§5-L4a).
7 the 4s→3p swap basin; the fold ΔE is partition bookkeeping — totals
only (§5).
8 purified core substrate (§5-L3): machine-exact closures; the SA-radial
floor 4.5e-07; the swap attractor survives (one draw, 9.2e-02).

**Ti — CAS(4,6), mechanism table** (same columns):

| rung | E_ground | ΔE | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ | ‖[P,Lz]‖ | ‖[P,Lᵔ]‖ | note |
|---|---|---|---|---|---|---|---|---|
| rohf | −848.406489 | 0 | — | — | — | — | — | 9 |
| frac-scf | −847.296157 | +1110 | — | — | — | — | — | 2 |
| casci_avg | −848.481312 | +10.6 | 12.0 | 2.2e-14 | 4.1e-13 | 6e-15 | 8e-15 | 3 |
| casci_raw | −848.484093 | +10.6 | 12.0000 | 4.9e-02 | 2.6e-01 | 6e-15 | 8e-15 | 4 |
| casci_frac | −848.443405 | +48.0 | 12.0 | 1.2e-13 | 7.5e-13 | 6e-15 | 8e-15 | 2 |
| CASSCF (free) † | −848.434973 | **−31.3** | 0.05 | 4.3e-01 | 1.1e+00 | 2.8 | 4.0 | 5 |
| casscf_relax | −848.393394 | +10.4 | 12.0 | 2.3e-08 | 1.7e-07 | 3.2e-08 | 5.9e-08 | 6, 9 |
| sa † | −849.267643 | −10.1 | 12.003 | 1.1e+00 | 6.0e+00 | 1.4 | 1.5 | 10 |
| sa_radial | −848.404105 | −0.4 | 12.0 | 8.7e-04 | 3.5e-03 | 1.2e-03 | 7.4e-04 | 6 |
| pur family | −848.3934…−848.4099 | +10.4 / −3.5 / +1.0 | 12.0 | ≤1.1e-12 | ≤1.1e-12 | 6e-15 | 8e-15 | 11 |

9 anchor — character-forced filling + E_ref basin gate (§3.1); the gated
reference keeps the relax rung's start defect at 2.3e-08.
10 the exact 4s→3p swap (window p = 1/6); totals only (§5).
11 machine-exact closures (1e-13); the SA-radial floor 1.1e-04; sa_pur =
a partial swap (s-slot ⅔ s + ⅓ p, −3.8).

Carbon and silicon: the summary table above is the whole story at this
precision — the p² confinement recovered in full by the symmetric
radial rotation, the pur family at the exact class; the dated censuses
(flat-direction draws, the free-CASSCF tilt anatomy) live in
`elements/C/NOTES.md` and `elements/Si/NOTES.md`.

### 7.4 T5 — the frozen-scan verdicts: where the gain lives

Method: freeze everything, deform the window along one coordinate, and
re-diagonalize the CASCI in it — no optimizer anywhere, so the numbers
are deterministic. Two coordinates, defined once: the **radial
rotation** (m-symmetric: the same rotation applied to every member —
closure-preserving) and the **single-member tilt** (one member only,
angle δ — closure-breaking; §2.3). What each can do to the energy is
derived in §2.4's deformations block.

| Atom | radial rotation (m-symmetric) | single-member tilt |
|---|---|---|
| C | recovers everything (θ* = 0.61, −250.8 mEh) | exactly flat |
| Si | gains everything (θ* = 0.85, −126.7 mEh) | exactly flat |
| Sc | rises | descends (~−0.3 mEh at δ = 0.1; deeper with more members) |
| Ti | rises | descends (~−160 µEh at δ = 0.1; pair patterns vary) |

Reading: *flat* → the optimizer may break m at zero cost (a gauge);
*descending* → it is driven to break m (the Jahn–Teller-like descent,
§2.4); *rises* → the direction buys nothing. On the p² atoms symmetry
is free (the symmetric radial rotation recovers the entire window
gain); on the d atoms the best symmetric treatments sit ~1 mEh from
ROHF while the escape basin reaches ~−30 mEh by redefining the
window — and the ticket-24 amendment holds: those ceilings are window
statements (the grown m-complete lattice crosses below them; §7.5).
Per-atom anatomy (multi-member patterns, blind-spot probes, substrate
channels): `elements/*/NOTES.md`.

### 7.5 The window lattice — exact symmetry vs raw MOs at matched cost

Three window kinds on the same mean field, compared at **equal
determinant count** (equal cost): **m-complete AO windows** —
contracted shells or primitive bands (the granularity axis is measured
weak; which granularity a row ran is a code detail); **raw-ROHF MO
windows** — the open reference MOs plus virtuals ranked by valence-band
S-weight; and **natural-orbital (NO) group windows** — occupation
groups of the term-averaged density (methodology below). Every window
row carries the span certificate ‖[P,Lz]‖ / ‖[P,Lᵔ]‖: the AO and NO
windows close at ~1e-13; the raw-MO spans never close (O(1)–9 at
every size measured, 4–100 orbitals — the SCF's open orbitals are
m-delocalized, each added virtual adds fractional m-directions, so a
big CASCI on raw ROHF cannot return symmetric roots at any practical
size, and the certificate says so without running any CI). Full
per-atom lattices: `elements/<Atom>/TABLES.md`.

**Ti — the worked example** (ROHF reference −848.406 Eh, 3d²4s², the
gated reference; M_S = 1; energies in Eh; equal-determinant pairs
where both kinds ran):

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

What it shows, in one pass: the AO lattice is monotone in radial
content and reaches below ROHF (the d-atom window already correlates);
the raw-MO column wins the smallest pair only — by *representability*,
containing the reference determinant where every AO window of the
lower rows sits above the mean field; the NO column matches the AO
value at equal dets and exceeds it at the top pairs; ⟨L²⟩ is exact
on the AO and NO columns at every size and broken on the MO column at
every size. U tells the same story at f-scale, with a descent race
between configuration blocks and a reference determinant that is a
numerically pure ⁵S° state — a different term, not an m-broken
⁵L° (`elements/U/`).

**Why the source density decides the window (the NO column's
methodology).** NOs are the eigenvectors of a one-particle density
γ — occupations = eigenvalues — so *which* γ you diagonalize decides
what "top occupations" can select. Three levels in the (l, m) basis:
*generic* γ carries m↔m′ coherences; [γ, L_z] = 0 (any single M_L
state, any SCF reference) is m-diagonal but m-UNEQUAL — an occupation
cutoff keeps members and drops partners; [γ, L_a] = 0 for all three
generators (the average over a whole term's M_L set, inside one M_S
sector) makes γ, by Schur's lemma, a scalar on each m-set —
γ = R(radials) ⊑ c·𝟙_{2l+1} — so occupation selection ranks only
RADIALS and the window is rotation-closed by construction.
Operational rule: average the whole degenerate ground-term eigenspace
before diagonalizing (the M_L average at fixed M_S is already fully
invariant: spatial rotations don't touch spin, and γ is spin-summed).
Measured on U and Sc/Ti: NO group windows close at ~1e-13 with exact
⟨L²⟩ at every size, realize full-window quality at a fraction of the
dets, and — like every AO-span window — do not contain the ROHF
determinant (representability). Ticket 23, on the U NO windows,
answered the relaxation-vs-growth question: SA over the 17-root ⁵L
block reaches −50.4238 at exactly quantized ⁵L (42 mEh below the
⁵S° mean field) by purely RADIAL relaxation, where free single-root
optimization is destroyed by the same attractor — the full record:
`elements/U/TABLES.md` and `elements/U/NOTES.md`.

**What the tables say.**

1. *Small cost: the raw-ROHF windows win — by representability, not
correlation.* Their windows contain the reference determinant; the AO
windows do not (AO-span nesting is not determinant nesting — on U the
reference determinant's projection onto the 3f+2d span has norm
0.948), and every AO window below the top rows sits above the mean
field.
2. *The MO descent then saturates, correlation-starved* (mEh-class
total gain; ⟨L²⟩ never quantizes), while *the AO lattice descends
monotonically at exact ⟨L²⟩ and owns the top pair on every atom*.
3. *Both transition metals cross below their escaped CASSCF basins at
exact symmetry* — the radial-vs-escape gap (§7.4) is a WINDOW
statement: at a few thousand determinants the exact-symmetry lattice
reaches (Ti: crosses) the broken value.
4. *The NO column beats the shell lattice at every equal determinant
count measured* — a term-averaged-density NO group is a better radial
than the contraction, while staying exactly rotation-closed; the
raw-ROHF column keeps the smallest pairs by representability.

**Band mechanism, in brief.** A "band" is the exponent interval of
admitted primitives for one l. It is read off a measurement: (1) for
every primitive, compute the Mulliken S-weight of the open reference
MOs on its rows; (2) sort primitives by that weight — the peaks rank
themselves; (3) "kf" = the k heaviest primitives, accepted only if
they form one CONTIGUOUS exponent interval (the printed "(lo, hi)" is
read off the set, not chosen); (4) contiguity is the honesty guard —
a scattered selection would carve an artificial radial hole — and a
conclusion holds only if the neighbor steps (k±1) reorder nothing.

---

## 8. Dictionary and file map

- **Per-atom records** (committed): `elements/<Atom>/` — `TABLES.md`
  (term content, spectrum, window × treatment, scan verdicts) and
  `NOTES.md` (dated specifics). New atom = new element file, not
  README growth.
- **Primary §7 sources** (regenerated, under `results/`):
  `*_report.md/json` — §7.1–7.2; `*_ladder.md/json` — §7.3 (the full
  17-rung registry); `*_scan.md` — §7.4's curves; `*_sa_decomposition.md`
  — the §5-L4 window-vs-core attributions.
- **Regeneration, per table kind**: T1/T2 — `scripts/run_c.py` and the
  per-atom `run_atom` calls (recorded in each report header); T3 — the
  per-atom `run_ladder` calls (recorded in each ladder header); T5 —
  `scripts/radial_scan.py` + the scan calls in `results/*_scan.md`;
  the §7.5 lattices — `scripts/run_tm_lattice.py`, `run_tm_methb.py`,
  `run_u.py`, `run_u_methb.py`, `run_u_sa.py` (ticket 23).
- **Executable documentation**: the test suite (116 tests) — term
  tables, exact-class closure pins, basin branches. Plain-language
  companions: `RUNGS.md`, the glossary appendix.

## 9. Status and directions

**Status.** C, Si, Sc, Ti delivered end-to-end; U campaign complete
through ticket 23 (the relaxation-vs-growth verdict: the deep
correlation is radial relaxation, held at exact ⁵L by the ensemble,
destroyed by free single-root). Per-atom records:
`elements/<Atom>/`; cross-atom story: §7.

**Open physics.** The ⁵S° mechanism (U reference + free-run attractor;
Ti ⟨L²⟩ = 0 MO windows — same class). Inner core-pair valence (the
m-complete design; the accidental 6p probe moved ~−0.5 Eh). U
f³d¹-vs-f²d² ordering inside one NO window. NIST comparison.

**Method gaps (ticket 20).** Axes-product runners behind a
`Window.columns(mol, mf)` resolver; per-cell budgets for long runs;
`--fresh` / chkfile key versioning; the intent manifest (precursor:
the six13 printed header); the full-matrix slice guard; the needle
linear-dependence helper (fired once); promote
`closed_shell_core_energy` to the library.

**Forward directions.** Spin–orbit (splits the certified multiplets;
extends in place: Ĵ², Ĵz as determinant-basis matrices).
Transition moments / g-factors (between the certified eigenstates).
Scalar relativity (h₁ only; re-measure the reference basins once).
Further atoms = element files + data rows.

## Appendix. Plain-language glossary

One or two lines per term; the precise statement lives where the
reference points.

### The objects

- **atom config** — the input card: atom, basis, spin (2S), charge;
  everything else derived.
- **reference** — a solved mean field. ROHF is lopsided (it occupied
  one specific open-shell orbital); fractional-occupation is round by
  construction (§3.2).
- **core** — the orbitals we do not correlate; enters as the fold
  field.
- **fold** — the effective active-space Hamiltonian with the core's
  Coulomb–exchange field in h₁ (§3.2, explicit formula).
- **window** — the active space. Kinds: shell/primitive AO bands
  (m-complete), raw-SCF MOs (never closed), NO occupation groups
  (exact from a term-averaged density; §7.5).
- **sector** — determinants at fixed (n_α, n_β); label Ms =
  (nα−nβ)/2.

### The analysis

- **member** — the m-component of a shell (px, py, pz; §2.3).
- **m-completeness** — all members present, sharing the same radial
  mixture per l; l-purity alone is not enough (§2.3).
- **commutator norms** ‖[H, Lz]‖, ‖[H, L²]‖ — the symmetry meters
  (§2.4).
- **closure meter** ‖[P, Lz]‖, ‖[P, Lᵔ]‖ — the same question for the
  window span alone; as a pair, no blind spots (§2.3).
- **span overlap** — smallest principal cosine between two spans; the
  identity gauge (‖did the window become a different space?‖),
  independent of closure (§5-L4).
- **window purity** — l-content of each window orbital; necessary, not
  sufficient.
- **stage / stage facts** — a finished calculation and its one uniform
  analysis pass (§4).
- **ortho_err / DIVERGED** — max |C†SC − 1|; above 1e-6 every number
  of the rung is meaningless (§2.4).
- **PySCF cross-check** (max|dE|) — max over roots of |E_r(ours) −
  E_r(PySCF)|: zero = identical bookkeeping; jumps = a non-scalar raw
  core (§7.3).

### The experiment axis

- **ladder / rung** — the comparison series over orbital treatments /
  one step of it (`RUNGS.md`).
- **rotation mask** — which orbital pairs the optimizer may rotate;
  how physics questions are imposed (§3.3).
- **radial constraint** — the mask "same (l, member) only": pure
  radial breathing (§3.3).
- **gauge** — a description choice that does not change the physics
  (which member combination each window orbital is).
- **flat direction** — a deformation the energy does not see
  (second-order at a stationary point; §2.4).
- **escape / swap** — the two identity basins: window rebuilt from
  unimplied orbitals (no term) / identity traded 4s↔3p (a different
  partition; totals only) — §5-L4.
- **input-window baseline** — the ~1e-3 closure defect inherited from
  core-orthogonalizing the start window against a non-scalar core
  (§5-L4); fixed by the closed core span.
