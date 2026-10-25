# AtomTerm

Atomic term-symbol analysis on top of [PySCF](https://pyscf.org): the pipeline
runs ROHF → CASCI → CASSCF → SA-CASSCF, then classifies every state by
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
python scripts/run_c.py       # writes results/ (chkfiles + report.md/.json)
```

Run from the repo root — `pyproject.toml` preconfigures `pythonpath = src`
and the test paths; no install step.

## What you get

Two result types, all written under `results/` (§7 maps every file to
its regeneration command).

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
rotation-closed (§1.3). The Ŝz line documents the sector's M_S; it
certifies nothing, since [Ĥ, Ŝz] ≡ 0 identically (§1.2). §6.1–6.2;
regeneration: §7.

**Treatment rungs** — every orbital treatment priced against the ROHF
anchor, the free variant beside the symmetric one, each classified by
its diagnostics (E as run = the run's own variational total; ΔE vs the
ROHF anchor). Titanium's headline rows (ROHF reference −848.406489 Eh):

| treatment | E (as run), Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | verdict |
|---|---|---|---|---|---|
| CASCI, averaged fold | −848.394046 | +12.4 | 12.0 | 2.2e-14 | exact |
| CASSCF, free | −848.437740 | −31.3 | 0.0503 | 4.3e-01 | escape — broken |
| CASSCF, window frozen (core-relax) | −848.395884 | +10.6 | 12.0 | 2.3e-08 | baseline |
| SA-CASSCF, free | −848.416563 | −10.1 | 12.003 | 1.1e+00 | swap — broken |
| SA-CASSCF, frozen, purified core | −848.394753 | +11.7 | 12.0 | 5.2e-13 | exact |

Verdicts: **exact** = commutators and span certificates at the 1e-13
class; **escape / swap** = window-identity basins (§4); **baseline** =
a named non-closure mechanism at its signature scale (§1.5). §6.3 is
the reading guide; the full per-atom tables (span certificates,
cross-check column, configurations) live in
`elements/<Atom>/TABLES.md`.

---

## Introduction

The exact states of an isolated atom carry sharp term symbols (³P, ¹D, …):
the Hamiltonian commutes with all rotations — the SO(3) × SU(2)
invariance of coordinates and spins — so every state has definite L and
S, and the states of one term are exactly degenerate. CASCI and CASSCF
treatments *can* lose this structure. What breaks is not quantum
mechanics but the *restriction*: the calculation is confined to a small
set of orbitals, and two properties of that restriction can fail.

First, CASCI inherits an anisotropic field from its reference (§1.4).
ROHF calculation for Carbon: the two ³P electrons occupy two of the
three 2p orbitals — the choice is arbitrary, all pairs are degenerate —
say px and py; their charge is then unevenly distributed (heavier along
x and y than z), the closed-shell core polarizes in response, and
the resulting field acts unevenly on the valence states, splitting the
³P term.

Second, CASSCF additionally lets the active orbitals themselves relax
(§2.3, §1.4). The active orbital set — the *window* (§1.3) — can then
deform in three ways with different consequences: some deformations
leave the energy unchanged, and the optimizer stops at an arbitrary
point of them; some *lower* the energy by tilting orbitals (§1.4); and
some change *which* orbitals are active at all (§4).

State-averaged CASSCF exists because a single optimized state represents
just one component of the multiplet; the average over all components is
the rotationally symmetric object (§2.4). It restores sharp L values,
but it cannot alone guarantee the full symmetry, and on heavier atoms it
can still trade orbitals.

Heavier atoms add failure modes (§4): the energy-lowering tilts, the
orbital trades, and the reference calculation itself can converge to a
wrong or metastable solution (§4).

Is the minimum correct? Where the energy is insensitive to a
deformation, the energy is still well defined — the orbitals are not,
and that is fine (§1.4). Across genuinely different solutions, every
free run must be inspected for which orbitals it actually uses before it
is believed (§6.3).

The deliverable is threefold: **certified term assignments** (the three
questions of §1.1); **comparisons** — every orbital treatment priced
against one common anchor (E as run, ΔE vs the ROHF reference), across
treatments, windows at matched determinant cost, and atoms (§6); and
the **loss accounting** — what each restriction and each symmetry
breaking costs, and which of those losses are removable (§4, §6). The
treatment class is CASCI/CASSCF, non-relativistic throughout; the
roadmap (§8) extends to scalar relativity and ECP-based full relativity.

---

*How to read.* The header, the showcase and the Introduction state the
goal and tell the story. §1 fixes the objects, the language
and the diagnostics — later sections use nothing that is not defined
there. §2 walks the method hierarchy (why each exists, what it computes,
how we run it, and what textbook usage gets wrong). §3 maps the methods
onto the code. §4 is the core: where symmetry breaks, level by level.
§5 is the operational protocol card. §6 holds the results. §7 is the
file map; §8 the roadmap. Plain-language companions: `RUNGS.md` and the
glossary appendix below; agent instructions: `AGENTS.md`.

---

## 1. The object and the diagnostics

### 1.1 Atom symmetries

An isolated atom is invariant under rotations. Its Hamiltonian therefore
commutes with L̂², L̂z, Ŝ², Ŝz, and every exact eigenstate is
|γ; L, M_L, S, M_S⟩ — it carries a **term symbol** (2S+1)L. Because of
L̂² and Ŝ² symmetry, the (2L+1)(2S+1) components of one term are *exactly
degenerate*.

The goal of this project is to understand **why practical CAS stages
(CASCI over an SCF reference, CASSCF, SA-CASSCF) can lose this symmetry**,
and how to follow it without giving up the correlation these
methods buy. A treatment is examined through three questions:

1. does every diagonalized root carry a sharp term label
   (⟨L̂²⟩ = L(L+1), ⟨Ŝ²⟩ = S(S+1))?
2. how large are the commutators ‖[H, L̂z]‖, ‖[H, L̂²]‖ evaluated as
   exact (§1.5) matrices over the determinant basis?
3. does the term table (state counts per term) equal the combinatorial
   free-ion table?

"Certified" is the word for exactly this, throughout this note: the
three questions answered against the classes of §1.5 — sharp term
labels, exact-class commutators, exact term tables — not approximate
agreement.

### 1.2 The operators, over determinants

We work with spin-orbitals (α = 0..n−1, β = n..2n−1) and
**determinants** — sorted occupation tuples. A **sector** fixes the
numbers (n_α, n_β) of α/β electrons, hence M_S = (n_α−n_β)/2 and N. The
operators N̂, L̂z, L̂², Ŝz, Ŝ² all preserve sectors (each commutes with Ŝz
and N̂), so every matrix lives inside one sector of fixed M_S.

A CAS treats a chosen set of active orbitals — the **window** (analyzed
in §1.3). Every operator is brought in as **one-body symbol + two-body
symbol**: the symbol is the operator's matrix in the active-orbital
basis, built from its AO-basis matrix via the orbital coefficients; both
implementations below consume the same symbols. Evaluation is by
Slater–Condon rules over determinants (`detbasis`), or — elementwise
identically, pinned by test — by explicit creation/annihilation on
occupation vectors (`secondquant`). The two implementations validate each
other.

Every operator enters the determinant space as a **one-body part** and,
where it has one, a **two-body part**,

  O₁ = Σ_pq o_pq a†_p a_q,   O₂ = ½ Σ_pqrs g_pqrs a†_p a†_q a_s a_r,

with these symbols (spin-orbitals ordered (α; β); L̂^a = C† $L̂^a_{AO}$ C,
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
Ĥ  :  o = h_eff  (the folded core field, §2.2)
      g = (pq|rs)   (the active-orbital Coulomb integrals)
```

All of it exact algebra — the L̂²/Ŝ² two-body symbols are products of
one-body matrix elements, nothing fitted or truncated. Ŝz carries no
certificate: Ĥ is spin-free, so [Ĥ, Ŝz] ≡ 0 identically; it only labels
the sector (in a fixed-M_S sector its matrix is exactly M_S·**1**), and
the symmetry certificates stay ‖[Ĥ, L̂z]‖ and ‖[Ĥ, L̂²]‖.

The ladder operators are built explicitly, for algebraic cross-checks.
Ŝ± are an **inter-sector map** — the one operator class that does not
preserve a sector (a bare ladder flips one β → α); L̂± are ordinary
one-body window operators (orbital ladders do not touch spin). They
have no analysis role, but they turn the operator algebra into pinned
tests: L̂² = L̂z² + ½(L̂₊L̂₋ + L̂₋L̂₊) and Ŝ² = Ŝz² + ½(S₋S₊ + S₊S₋) hold
elementwise against the delivered operators on every p² sector, and
[L̂z, L̂±] = ±L̂±.

### 1.3 The window: rotation closure, deformations, certificates

**We call a space rotation-closed if L̂W ⊆ W, and an operator a scalar if
it commutes with all rotations.** The members of an l-shell are its 2l+1
m-functions, indexed within the shell (0…2l) in AO order:

| p member | index | d member | index | f member | index |
|---|---|---|---|---|---|
| px | 0 | dxy | 0 | fy(3x²−y²) | 0 |
| py | 1 | dyz | 1 | fxyz | 1 |
| pz | 2 | dz² | 2 | fyz² | 2 |
| | | dxz | 3 | fz³ | 3 |
| | | dx²−y² | 4 | fxz² | 4 |
| | | | | fz(x²−y²) | 5 |
| | | | | fx(x²−3y²) | 6 |

(f abbreviations: yz² = y(3z²−r²), z³ = z(5z²−3r²); PySCF's own labels
for f members are the m values, `4f−3 … 4f+3`.) The ±m partners —
(px, py); (dxz, dyz), (dx²−y², dxy); (fxz², fyz²), (fz(x²−y²), fxyz),
(fx(x²−3y²), fy(3x²−y²)) — are the pairing behind the L̂z blind spots
below.

**Exact symmetry imposes two conditions on the treatment.**

1. **The window is rotation-closed.** Otherwise L̂ maps window functions
   outside the window, and the CAS Hamiltonian inherits the leak:
   [H, L̂] ≠ 0 even though the physical atom is symmetric. Violations
   arise in practice because canonical-MO windows and optimizer output
   carry member-dependent radial mixtures, because truncation-by-weight
   can drop members, and because even core-orthogonalization of a raw
   shell window can inject non-closed directions (the mechanisms, with
   numbers, are §4's subject).
2. **The core fold is a scalar.** The inactive core enters the active
   problem as a one-body field inside h₁; if that field is angularly
   anisotropic, it breaks L̂z directly. This happens because the SCF core
   polarizes in the anisotropic field of the open shell (§1.4).

**Rotation generators never touch radial coordinates.**

> **Theorem (closure).** For φ(r)·Y_lm one has L̂(φ·Y) = φ·(L̂Y) — the
> radial factor is untouched. Therefore a finite one-particle space W is
> rotation-closed **iff**, for each l present, W contains the complete
> set of 2l+1 members, all sharing the *same* radial subspace E:
> W = ⊕ₗ (E_l ⊗ {all members of l}). The radial subspace E itself is
> arbitrary.

**Exactly two failure modes exist: members missing, and radials not
shared across members.**

**Two projector norms check the window alone.** With no CI and no
optimizer involved, ‖[P, L̂z]‖ and ‖[P, L̂⊥]‖ — L̂⊥ = (L̂x, L̂y), the two
components perpendicular to z — P the orthogonal
projector onto the window span in the AO overlap metric — catch exactly
the two failure modes above.

**A member tilt is the elementary window deformation.** One member is
replaced by a radial admixture, cos θ·(n l m) + sin θ·(n′ l m) — same
l, same m, neighboring radial. In the real AO basis L̂z mixes only the
two functions of each ±m pair and annihilates m = 0 members, so
‖[P, L̂z]‖ has **blind spots** (an m = 0 tilt, or equal tilts of both
partners of a pair, are invisible to it); L̂⊥ rotates any member into
any other, so ‖[P, L̂⊥]‖ sees every member asymmetry — together, no
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

### 1.4 The diagnostics and their sensitivity

**ε measures the size of the symmetry breaking.** Throughout, it is the
weight of the angularly non-scalar part of the core field's density (a
scalar core field contributes nothing to any commutator), reported as
core anisotropy; ‖A‖_F = sqrt(Σ|A_ij|²) is the Frobenius norm.

| Diagnostic | What it measures | Sensitive to | Blind to |
|---|---|---|---|
| ‖[H, L̂z]‖_F, ‖[H, L̂²]‖_F | the assembled treatment (window + core fold) | everything, order ε | vacuous in single-irrep sectors (L̂²); ±pair structure (L̂z) |
| ‖[P, L̂z]‖, ‖[P, L̂⊥]‖ | the window span alone (members §1.3) | member asymmetry, missing members; no optimizer noise | — (as a pair) |
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
shell's quadrupole polarizes the core during the SCF (§1.3, condition
2); the cures — spherical averaging of the core density, the
fractional-occupation reference, core purification — are compared by
cost in §6.

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
### 1.5 Numerical conventions — what “exact” costs

| Class | Typical scale | Meaning | Always resolvable? |
|---|---|---|---|
| exact | ~1e-13 | round-off of the construction | — |
| baseline | 1e-8…1e-2 | a *named* mechanism: core-orthogonalization defect, core polarization, optimizer floor | each has a signature (checkable); all resolvable except the optimizer floor — certifiable and boundable, but with no restoring force to remove it (§4) |
| broken | 1e-1…10 | window identity change (swap/escape, §4) | yes — recognized and rejected by the diagnostics |

§4 maps every rung to a class.

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
  force (§2.4).


---

## 2. The method hierarchy

Each method gets the same triple: why it exists, what it computes, how we
run it. Two vocabulary items are defined here and used everywhere after:
the *substrate* (§2.3) and the *basin* (§2.3).

### 2.1 ROHF — the reference that defines "core"

**The reference calculation fixes what the core is: everything downstream
folds it.** Two of its aspects matter here: how an orbital's *type* is measured,
and how the reference is kept honest where energy-ordered filling
fails.
**An orbital's type is measured by projection.** For an AO span with
index set o, the S-metric projector is P_o = S[:,o] S_oo⁻¹ S[o,:], and an
orbital c carries weight w = c†P_oc ∈ [0,1] on that span (bounded — a
true projector weight; the same metric appears in §2.2's window
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
shell's quadrupole polarizes the core during the SCF (§1.4); the fold
variants of §2.2 are three ways of dealing with what this calculation
leaves behind.

### 2.2 CASCI — diagonalization at fixed core and window

**CASCI diagonalizes the CI problem inside one sector, at fixed core and
window.** The window is the raw shell AOs (rotation-closed by
construction, §1.3).

**Core, fold, and the exactness condition.** *Core* means the
closed-shell (doubly occupied) orbitals the treatment does not
correlate — a subset of the *electrons*; it is not an effective core
potential (where an ECP is in play — U — it has already replaced the
innermost shells before any partition is made). The core enters the
active Hamiltonian only through the **folded field** — the Ĥ one-body
symbol of §1.2:

  h_eff(p,q) = h(p,q) + Σᵢ[2(pq|ii) − (pi|iq)],  E_core kept separate,

sum over occupied core orbitals i. The partition is **exact iff the
window is orthogonal to the core orbitals**. A window vector leaning
into the core, w = w⊥ + w∥ with w∥ ∈ span(core), admits configurations
with extra density in *occupied* core orbitals that the antisymmetrized
problem forbids — the CI harvests deep-core stabilization the Pauli
principle would have priced, silently sub-variational. Hence every
window is built core-orthogonal (the §2.3 substrate step); core
averaging changes the field, never the orthogonality of the basis. The
three constructions below differ exactly in which core density feeds
this fold:

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
  multiplet splits by §1.4's mechanism.

With the averaged fold, the CASCI term tables of §6.1–6.2 come out exact
on all four atoms; the averaged–raw and averaged–frac differences are
the ΔE anchors of §6.3.

### 2.3 CASSCF — letting the orbitals move

**CASSCF minimizes the chosen root's energy over orbital rotations.** Its
motivation is visible in carbon's radial gain — the valence 2p functions
want a 3p admixture worth −250.8 mEh of window energy (§6.4 — in totals
this is the recovery of the raw-shell confinement cost, §6.3), a
relaxation CASCI cannot
express. How we run it is a grid of *starts × restrictions*, and both
axes matter.

**The substrate is the starting MO matrix the optimizer receives** —
core block, window block, complement. Ours is block-safe, and the window
step is explicit mathematics: a raw window column w is replaced by

  w̃ = (1 − P_core) w / ‖(1 − P_core) w‖_S,  P_core = C (C†SC)⁻¹ C†,

i.e. its component in the core span is subtracted in the AO overlap
metric and the remainder re-normalized. The projection is the leak guard
of §2.2: without it, a window that
leans into the core would double-count the core region. And the three
orthogonalizers of this project each preserve a different auxiliary
property (all in the S metric — §1.2's isometry turns them into the
Euclidean ones): symmetric **Löwdin** (shape-preserving — the core,
§2.2), two-pass S-**Gram–Schmidt** (identity-preserving — column j
is the (l, m) of AO j; the masks require it), and the residual
**completion** (span-filling — the complement). The construction is
forced by convergence — a non-orthonormal
substrate makes the optimizer diverge to sub-variational energies
(observed once; `elements/Si/NOTES.md`, the §1.4 guard).


**Restrictions are masks on rotation pairs.** *Window frozen*: only
core↔virtual pairs — correlation through the core alone, closure
untouched. *Radial restriction*: only pairs within a same-(l, m) radial
chain — the symmetry-preserving channel; where it exists (C, Si) it
recovers the full gain at exact closure. *Free*: no mask — the textbook
run. **A basin is the region of attraction of one local minimum of this
optimization**; the free landscape has several (§6.3's chaos paragraph).


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

### 2.4 SA-CASSCF — averaging over the multiplet

**State-averaging replaces the single root by the multiplet average:**

  E_SA(orbitals) = (1/(2L+1)) Σ_a E_a(orbitals),  a over the M_L
  components,

minimized over rotations (³P: three components; ²D: five; ³F: seven).
Where a single root represents one M_L component — and optimizes for it
alone — the average is the rotationally scalar object.

**The average anchors the term identity but not the window.** ⟨L̂²⟩
is restored, yet along flat directions the average is stationary with
no restoring force (§1.4: the deformation block is traceless), so
convergence tolerance, not physics, decides the endpoint; and other
basins remain reachable (the swap, §4). Exactness combines the tools:
state average + frozen window + isotropic substrate. Whether tilting
helps or hurts the average is atom-dependent — the frozen-scan
verdicts of §6.4 (per-atom values in `elements/`).


### 2.5 What textbook runs give you, and why

| Standard practice | What you get | Why | What we do instead |
|---|---|---|---|
| CASCI over canonical virtual MOs | closure defects at the 1e-1 level; term labels unreliable | canonical MOs diagonalize a Fock built from an anisotropic density; quasi-degenerate same-l members mix arbitrarily | raw shell windows, averaged fold (§2.2) |
| free CASSCF, default start | p²: full gain with a silently member-tilted window; d-shells: escape/swap basins, ⟨L̂²⟩ off by units | flat directions (endpoint = noise) and descending directions (§1.4); window identity free to change | restrictions + substrate control (§2.3) |
| free SA-CASSCF | ⟨L̂²⟩ ≈ L(L+1) restored, closure not (Sc: 6.002 with ‖[H,L̂z]‖ ≈ 0.8) | the average is scalar to first order only; no restoring force on flat directions; other basins still reachable | frozen-window SA on the isotropic substrate (§4) |

## 3. The pipeline tour

| Module | Owns |
|---|---|
| `system` | system validation — one atom at the origin, spherical basis (L̂ defined about r = 0) |
| `angular` | one-electron L̂x/L̂y/L̂z AO matrices — two independent constructions pinning each other (libcint; exact tesseral algebra) |
| `detbasis` / `secondquant` | many-body matrices over determinants — Slater–Condon and explicit a†/a, elementwise-identical twins (§1.2) |
| `operators` | operator assembly — the fold (§2.2), the H/L̂z/L̂²/Ŝ²/N̂ symbols, the window-closure certificates (§1.3) |
| `reference` / `casscf` | the engine — SCF reference with gates; the one deep optimizer (`constrained_casscf`) with its masks and adapters |
| `windows` / `core` | window construction, selection, the kind registry (`Window.columns`); the purified core substrate and fold energy |
| `stages` | the one uniform analysis pass per finished stage — the source of every §6 number |
| `scans` | frozen-landscape measurements — deterministic (§6.4) |
| `analysis` | term labels, configuration summaries, report rendering |
| `cases` / `treatments` / `study` | orchestration — the declared per-atom campaign (`cases`: CaseSpec → run_case); rungs as data (`treatments`, §6.3); the declared window × treatment study runner (`study`) |

**The drivers capability map.** Every treatment is a point in a
five-axis space:

| Capability | What it is |
|---|---|
| ROHF reference | character-forced filling + energy gate — the anchor (§2.1, §4) |
| CASCI | exact sector diagonalization — fold variants raw / averaged / fractional (§2.2) |
| CASSCF | orbital optimization — freedom by mask, start by substrate (§2.3) |
| SA-CASSCF | multiplet-averaged optimization — weights over the M_L components (§2.4) |
| core construction | the purified, exactly isotropic core span (§4) |

The 17-rung registry tabulated in `results/` enumerates the corners of
this space; §6.3 keeps the rows that differ in mechanism.

**The launch path.** A new atom is a declared case, executed top-down —
per-atom *data* (closed/open AO spans, the reference energy E_ref, the
combinatorial term tables), not per-atom code:

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
recorded — the gate then rejects the known-wrong basins on every later
run.

---

## 4. Where symmetry breaks

| Level | Symptom (order of magnitude) | Mechanism | Fix | Cost |
|---|---|---|---|---|
| 1. Basis bookkeeping | labels mislead; shells fuse | wrong inputs → wrong constructions | label + metric discipline | none |
| 2. Sector witnesses | ‖[H,L̂²]‖ = 0 with defects present | single-irrep sectors make L̂² = L(L+1)·**1**; ±pair blind spots (§1.3) | projector pair; multi-term sectors | sector size |
| 3. Core field | splittings and commutators 1e-5…1e-1 | polarization chain (§1.4) | fold cures / purification | ΔE per §6 |
| 4. Orbital optimization | tilted windows, floors, escape/swap | flat and descending directions (§1.4); identity change | restrictions + substrate | 0…~30 mEh above the deepest minimum |
| 5. SCF reference | ~1e2 mEh wrong configuration; ~10 mEh metastable basin | energy-ordered filling; degenerate-open-shell basins | toolbox of §2.1 | per-atom reference data |

Class mapping (§1.5): *exact* — raw-window CASCI, frozen-window rungs
on the isotropic substrate; *baseline* — the start-window defect and
the flat-direction floor; *broken* — the escape and the swap (L4).

**L1 — basis bookkeeping.** AO labels are truthful yet treacherous:
pattern matching is regex-flavored; shells fuse (one label, several
contracted functions — unpack through contraction data); same-l shells
of different n overlap (~0.3), so every projection and normalization
runs in the AO overlap metric, never the Euclidean one.

**L2 — sector witnesses.** A single-irrep sector makes ‖[H, L̂²]‖ vacuous
and L̂z has ±pair blind spots (§1.3); the projector pair is the
universal certificate, multi-term sectors the honest courts.

**L3 — the core field.** The reference leaves an anisotropic field
acting at first order on the multiplet (§1.4). Two cures act at the
density level — spherical averaging, the fractional-occupation
reference. The third, core purification, exists because of the
*window*, not the core: the start window is (1−P_core)-projected
(§2.3), so its closure inherits the core's; rebuilding the core span
as E_l ⊑ {all members of l} (§1.3's theorem applied to the core)
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

**L5 — the SCF reference.** The §2.1 toolbox; E_ref is extracted once
per atom by the diagnostics themselves (occupations by projection,
core purity, window closure) — self-calibrating, per-atom data, not
code.

---

## 5. The protocol card

**Run all five deliverables, report all five** — the certified treatment
and the direct values are complements, not alternatives. Each deliverable
is a named configuration in the §3 capability space; each optimizer
deliverable is reported as the pair *free variant / symmetric variant*.

| Deliverable | Configuration | Purpose | Acceptance |
|---|---|---|---|
| ROHF reference | forced filling + gate | the anchor every ΔE hangs on | occupations & core purity by projection; E = E_ref |
| CASCI anchor | raw shell window, averaged fold | the exact term tables | tables = combinatorial; commutators exact |
| Certified treatment | SA-CASSCF + window frozen + purified core | exact-class symmetry with correlation | commutators & window spans exact; tables intact |
| Direct CASSCF | free, shell substrate | the energy-seeking value + basin anatomy | classified (span overlap, window weights); totals quoted |
| Direct SA-CASSCF | free SA-CASSCF, shell substrate | best direct value; term-readable even in the swap basin | term label; basin classified; totals quoted |

For a new atom: (1) run the reference (spans and E_ref are the atom's
data, §2.1/§4); (2) the CASCI anchor — term tables must come out exact
(§6.1–6.2); (3) the certified treatment (expect the exact class, §6.3's
pur rows); (4) both direct runs, classified; (5) optionally, frozen
scans for the radial/tilt anatomy (§6.4).

---

## 6. Results tables

Five table kinds, each answering one question. **T1 term content**
(§6.1) — the combinatorial contract the exact stage must return;
**T2 spectrum** (§6.2) — the sharp CAS-level term energies; **T3
treatment × window** (§6.3) — what each orbital treatment buys and
costs, per window; **T5 frozen scans** (§6.4) — where the gain lives,
optimizer-free; **the window lattice** (§6.5) — window kinds at matched
determinant cost. Per-atom records live in `elements/<Atom>/TABLES.md`
(+ `NOTES.md` for the dated specifics); **Ti is the worked example**
below. All numbers are this project's measurements (C, Si, Sc, Ti on
cc-pVDZ; U on ECP + cut-Dyall); the artifacts live under `results/`
(regeneration in §7). On symmetry-broken basins the term table itself
scrambles — independent evidence beyond the commutators.

### 6.1 Term tables — what each sector must contain

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

### 6.2 Multiplet energies — the sharp, degenerate spectra

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

### 6.3 T3 — the treatment rungs, one row per mechanism

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

**Sc and Ti — the mechanism tables** live in
`elements/Sc/TABLES.md` and `elements/Ti/TABLES.md` (the T3b format:
E as run, ΔE, ⟨L²⟩, the commutator and span-certificate columns,
per-rung notes). Their headline on the unified ruler: the certified
rungs sit mEh-tight around the frozen value (CASCI +10.5 ≈
core-relax +10.4 ≈ SA-pur +11.7 on Ti) — orbital optimization buys
≈nothing on the as-run window; the negative free/SA numbers are
window-identity changes, classified by the span certificates. Ticket
25: the fold convention is machine-exact; the cross-check column is
the mEh-class averaging effect — the old ~90 mEh window-definition
offset was the retired raw-AO analysis convention's artifact.

Carbon and silicon: the summary table above is the whole story at this
precision — the p² confinement recovered in full by the symmetric
radial rotation, the pur family at the exact class; the dated censuses
(flat-direction draws, the free-CASSCF tilt anatomy) live in
`elements/C/NOTES.md` and `elements/Si/NOTES.md`.

### 6.4 T5 — the frozen-scan verdicts: where the gain lives

Method: freeze everything, deform the window along one coordinate, and
re-diagonalize the CASCI in it — no optimizer anywhere, so the numbers
are deterministic. Two coordinates, defined once: the **radial
rotation** (m-symmetric: the same rotation applied to every member —
closure-preserving) and the **single-member tilt** (one member only,
angle δ — closure-breaking; §1.3). What each can do to the energy is
derived in §1.4's deformations block.

| Atom | radial rotation (m-symmetric) | single-member tilt |
|---|---|---|
| C | recovers everything (θ* = 0.61, −250.8 mEh) | exactly flat |
| Si | gains everything (θ* = 0.85, −126.7 mEh) | exactly flat |
| Sc | rises | descends (~−0.3 mEh at δ = 0.1; deeper with more members) |
| Ti | rises | descends (~−160 µEh at δ = 0.1; pair patterns vary) |

Reading: *flat* → the optimizer may break m at zero cost (a gauge);
*descending* → it is driven to break m (the Jahn–Teller-like descent,
§1.4); *rises* → the direction buys nothing. On the p² atoms symmetry
is free (the symmetric radial rotation recovers the entire window
gain); on the d atoms the best symmetric treatments sit ~1 mEh from
ROHF while the escape basin reaches ~−30 mEh by redefining the
window — and the ticket-24 amendment holds: those ceilings are window
statements (the grown m-complete lattice crosses below them; §6.5).
Per-atom anatomy (multi-member patterns, blind-spot probes, substrate
channels): `elements/*/NOTES.md`.

### 6.5 The window lattice — exact symmetry vs raw MOs at matched cost

Three window kinds on the same mean field, compared at **equal
determinant count** (equal cost): **m-complete AO windows** —
contracted shells or primitive bands (the granularity axis is measured
weak; which granularity a row ran is a code detail); **raw-ROHF MO
windows** — the open reference MOs plus virtuals ranked by valence-band
S-weight; and **natural-orbital (NO) group windows** — occupation
groups of the term-averaged density (methodology below). Every window
row carries the span certificate ‖[P,Lz]‖ / ‖[P,L⊥]‖: the AO and NO
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
exact symmetry* — the radial-vs-escape gap (§6.4) is a WINDOW
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

## 7. Dictionary and file map

- **Per-atom records** (committed): `elements/<Atom>/` — `TABLES.md`
  (term content, spectrum, window × treatment, scan verdicts) and
  `NOTES.md` (dated specifics). New atom = new element file, not
  README growth.
- **Primary §6 sources** (regenerated, under `results/`):
  `*_report.md/json` — §6.1–6.2; `*_ladder.md/json` — §6.3 (the full
  17-rung registry); `*_scan.md` — §6.4's curves; `*_sa_decomposition.md`
  — the §4-L4 window-vs-core attributions.
- **Regeneration, per table kind**: T1/T2 — `scripts/run_c.py` and the
  per-atom `run_atom` calls (recorded in each report header); T3 — the
  per-atom `run_ladder` calls (recorded in each ladder header); T5 —
  `scripts/radial_scan.py` + the scan calls in `results/*_scan.md`;
  the §6.5 lattices — `scripts/run_tm_lattice.py`, `run_tm_methb.py`,
  `run_u.py`, `run_u_methb.py`, `run_u_sa.py` (ticket 23).
- **Executable documentation**: the test suite (128 tests) — term
  tables, exact-class closure pins, basin branches. Plain-language
  companions: `RUNGS.md`, the glossary appendix.

## 8. Status and directions

**Status.** C, Si, Sc, Ti delivered end-to-end; U campaign complete
through ticket 23 (the relaxation-vs-growth verdict: the deep
correlation is radial relaxation, held at exact ⁵L by the ensemble,
destroyed by free single-root). Per-atom records:
`elements/<Atom>/`; cross-atom story: §6.

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
  construction (§2.2).
- **core** — the orbitals we do not correlate; enters as the fold
  field.
- **fold** — the effective active-space Hamiltonian with the core's
  Coulomb–exchange field in h₁ (§2.2, explicit formula).
- **window** — the active space. Kinds: shell/primitive AO bands
  (m-complete), raw-SCF MOs (never closed), NO occupation groups
  (exact from a term-averaged density; §6.5).
- **sector** — determinants at fixed (n_α, n_β); label Ms =
  (nα−nβ)/2.

### The analysis

- **member** — the m-component of a shell (px, py, pz; §1.3).
- **m-completeness** — all members present, sharing the same radial
  mixture per l; l-purity alone is not enough (§1.3).
- **commutator norms** ‖[H, Lz]‖, ‖[H, L²]‖ — the symmetry meters
  (§1.4).
- **closure meter** ‖[P, Lz]‖, ‖[P, L⊥]‖ — the same question for the
  window span alone; as a pair, no blind spots (§1.3).
- **span overlap** — smallest principal cosine between two spans; the
  identity gauge (‖did the window become a different space?‖),
  independent of closure (§4-L4).
- **window purity** — l-content of each window orbital; necessary, not
  sufficient.
- **stage / stage facts** — a finished calculation and its one uniform
  analysis pass (§3).
- **ortho_err / DIVERGED** — max |C†SC − 1|; above 1e-6 every number
  of the rung is meaningless (§1.4).
- **PySCF cross-check** (max|dE|) — max over roots of |E_r(ours) −
  E_r(PySCF)|: zero = identical bookkeeping; jumps = a non-scalar raw
  core (§6.3).

### The experiment axis

- **ladder / rung** — the comparison series over orbital treatments /
  one step of it (`RUNGS.md`).
- **rotation mask** — which orbital pairs the optimizer may rotate;
  how physics questions are imposed (§2.3).
- **radial constraint** — the mask "same (l, member) only": pure
  radial breathing (§2.3).
- **gauge** — a description choice that does not change the physics
  (which member combination each window orbital is).
- **flat direction** — a deformation the energy does not see
  (second-order at a stationary point; §1.4).
- **escape / swap** — the two identity basins: window rebuilt from
  unimplied orbitals (no term) / identity traded 4s↔3p (a different
  partition; totals only) — §4-L4.
- **input-window baseline** — the ~1e-3 closure defect inherited from
  core-orthogonalizing the start window against a non-scalar core
  (§4-L4); fixed by the closed core span.
