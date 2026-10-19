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
summary card (ROHF reference −848.406 Eh = −23 086 eV; ΔE vs ROHF, eV):

| deliverable | free variant | symmetric variant |
|---|---|---|
| CASCI, raw shell | — | +0.29 (exact) |
| CASSCF | **−0.85** (escape — broken) | +0.28 frozen (baseline) / −0.06 radial (single root, gauge-chaotic) |
| SA-CASSCF | −0.27 (swap — broken) | +0.03 radial purified (floor) / +0.32 frozen purified (exact) |

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
(§3.3, §3.6). The active orbital set — the *window* (§2.3) — can then
deform in three ways with different consequences: some deformations
leave the energy unchanged, and the optimizer stops at an arbitrary
point of them; some *lower* the energy by tilting orbitals (on d-shells;
the Jahn–Teller mechanism, §3.6); and some change *which* orbitals are
active at all — trading a valence 4s for a core 3p (§5-L4c).

State-averaged CASSCF exists because a single optimized state represents
just one component of the multiplet; the average over all components is
the rotationally symmetric object (§3.4). It restores sharp L values,
but it cannot alone guarantee the full symmetry, and on heavier atoms it
can still trade orbitals.

Heavier atoms add failure modes (§5): the energy-lowering tilts, the
orbital trades, and the reference calculation itself can converge to a
wrong or metastable solution (§5.5).

Is the minimum correct? Where the energy is insensitive to a
deformation, the energy is still well defined — the orbitals are not,
and that is fine (§3.6). Across genuinely different solutions, every
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

Why L̂² has a two-body part at all: total angular momentum is a *sum over
electrons*, L̂ = Σᵢ L̂ᵢ, and the square of a sum contains cross terms —
L̂² = Σᵢ L̂ᵢ² + 2Σᵢ<ⱼ L̂ᵢ·L̂ⱼ. N̂ and L̂z are plain sums (no square, no
cross terms). Ŝ² is also a square of a sum, but spin operators never
touch spatial coordinates, so their cross terms reduce to pure
pair-combinatorics over the occupation tuples.

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
shared across members.** A window where dxy is built from one radial
mixture and dx²−y² from another is not closed, however pure its
l-content.

**Window operators never leave the window Fock space.** Constructing
"L̂²" from a†_p a_q (p, q ∈ W) therefore gives the *projected* L̂²; it
equals the true L̂² on W's determinants iff W is rotation-closed —
otherwise the truncated symbol silently reports wrong quantum numbers
while everything looks perfectly well-defined.

**Two projector norms check the window alone.** With no CI and no
optimizer involved, ‖[P, L̂z]‖ and ‖[P, L̂⊥]‖ — P the orthogonal
projector onto the window span in the AO overlap metric — catch exactly
the two failure modes above.

**A member tilt is the elementary window deformation.** One member is
replaced by a radial admixture, cos θ·(n l member) + sin θ·(n′ l member)
— same l, same m, next radial shell; for instance, tilted
dxy = cos θ·(3dxy) + sin θ·(4dxy). The tilt is the unit in which member
asymmetry is analyzed: a general window deformation decomposes into
l-content changes, member structure, and radial content, and the tilt is
the cleanest coordinate of the member-radial part. Real optimizer output
need not be a single tilt — on Ti it demonstrably is not (§5) — but its
member structure is read in these terms. Visibility of tilts is
asymmetric in a way that can be stated exactly. PySCF's AO functions are
real combinations of the spherical harmonics; in this basis L̂z mixes only
the two functions of each ±m pair (px↔py, dxy↔dx²−y²) and annihilates the
m = 0 members (pz, dz²). So a tilt of pz alone, or *equal* tilts of both
partners of a ±pair, leave the window L̂z-closed: ‖[P, L̂z]‖ has **blind
spots** — defect classes it cannot detect. L̂⊥ = (L̂x, L̂y) rotates any
member into any other, so ‖[P, L̂⊥]‖ sees every member asymmetry;
together the two norms have no blind spots. Both directions of this
asymmetry have been observed. On Ti, a deliberately constructed
four-member tilt (all but dz² — equal tilts of both ±pairs, applied to a
frozen window as a deterministic probe) kept ‖[H, L̂z]‖ at 1e-14 while
‖[H, L̂²]‖ fired at ≈ 1–2. On carbon, by contrast, the deformed window
was the optimizer's own product: an optimized window tilted in an
L̂z-invisible way kept even the many-body commutators quiet
(‖[H, L̂z]‖ ≈ 5e-7) while ‖[P, L̂⊥]‖ ≈ 1.1 exposed it.

**A flat direction is a deformation along which the energy is constant to
numerical precision.** Along such a direction the optimizer's endpoint is
decided by numerical noise alone; the endpoint is also called a *gauge
choice* — energetically indistinguishable alternatives, distinguishable
only through the closure diagnostics. The member tilts of the p² atoms
are the working example (§7.3: five fresh optimizer draws share one
energy while ⟨L̂²⟩ ranges over 1.69–2.00).

**Sector content decides what its commutators can testify.** Any sector
decomposes into irreducible L-representations, and L̂² acts as L(L+1)·**1**
on each copy. If a sector carries a single irrep — carbon's ³P Ms=1
sector, three determinants — then L̂² is a multiple of the identity on
the whole sector and ‖[H, L̂²]‖ vanishes *whatever is wrong*: the check
is vacuous there, and multi-term sectors are the honest courts. The
four-member tilt above, invisible to L̂z, is caught here precisely
because this sector's many irreps give L̂² room to act. These examples
come from
the four atoms analyzed so far (C, Si, Sc, Ti); testing the picture on
further atoms is part of the roadmap (§9). The projector pair, being
sector-independent, is the universal certificate.

**The two conditions fail independently, and raw-shell CASCI probes both.**
CASCI never breaks window closure — the raw shell window is exactly
closed, so the L̂² we build is the true L̂² (on carbon's raw window its
eigenvalues are exactly {2, 2, 2}); a drifting ⟨L̂²⟩ like 1.72 requires a
*deformed* window (optimizer output) and cannot occur here. CASCI over
the raw ROHF fold, however, does carry the ε defect (§2.4) through the
core field, not the window.

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

**All five diagnostics respond to one object: the symmetry-breaking part
of H.** Write H = H₀ + εV, with H₀ exactly symmetric, ε the size of the
breaking and V its shape. ε is measured, not assumed: core anisotropy is
2.4e-3 (C), 2.1e-3 (Si), 7.8e-3 (Ti), 1.3e-2 (Sc) — percent-level at
worst — and the bookkeeping is verified: commutators and splittings
track ε to first order across the atoms, ⟨L̂²⟩ shifts stay orders below.

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

**Core anisotropy has a physical origin and three cures.** The origin is
condition 2 of §2.3: the open shell's anisotropic charge polarizes the
closed core during the SCF, and the folded field acts unevenly on the
valence M_L components. The cures — spherically averaging the core
density, the fractional-occupation reference, or the purified core span
— are compared in §7 by what each costs.

**Orthonormality is a validity guard, not a symmetry diagnostic.** Every
CAS stage consumes an MO matrix, and if it is not orthonormal in the AO
overlap metric the construction is invalid: on Si we watched a
non-orthonormal input make the optimizer diverge and return energies
*below* the variational bound while the symmetry diagnostics still looked
plausible. Each stage therefore reports its orthonormality error.

**The chain from quadrupole to splitting closes quantitatively.**
Inversion symmetry forbids any dipole component in the atom's charge
distribution, so the leading anisotropic term is the quadrupole; this
fixes the *shape* of the breaking, and carbon's ³P components indeed
split in the pattern {M_L = 0} vs {±1} — the quadrupole signature, by
which its form was identified in this project. The *size* comes from the
induced response: the valence quadrupole polarizes the core with
amplitude set by the coupling over the core excitation gap, and the
measured core anisotropy ε is precisely that induced quadrupole weight.
The splitting then follows at order ε, multiplied by a coupling prefactor
reflecting how strongly the anisotropic field couples to the particular
valence manifold:

| Atom | ε (core anisotropy, relative) | splitting (Eh) | splitting/ε (Eh) | ⟨L̂²⟩ of the ground term (raw fold) |
|---|---|---|---|---|
| C | 2.405e-3 | 4.373e-5 (³P) | 0.02 | 2.0000000000 |
| Si | 2.109e-3 | 2.431e-4 (³P) | 0.12 | 2.0000000000 |
| Ti | 7.754e-3 | 4.144e-3 (³F septet) | 0.53 | 12.0000000000 |
| Sc | 1.259e-2 | 1.857e-3 (²D quintet) | 0.15 | 6.0001 |

The ratio splitting/ε is the effective coupling strength in Eh per unit
of relative anisotropy — an atomic property, not a universal constant,
since the valence coupling geometry differs from atom to atom. The order
relationships, however, hold throughout: commutators and splittings are
the detectors; a numerically pretty ⟨L̂²⟩ is not.

### 2.5 What "exact" means numerically

| Class | Typical scale | Meaning | Always resolvable? |
|---|---|---|---|
| exact | ~1e-13 | round-off of the construction | — |
| baseline | 1e-8…1e-2 | a *named* mechanism: core-orthogonalization defect, core polarization, optimizer floor | each has a signature (checkable); all resolvable except the optimizer floor — certifiable and boundable, but with no restoring force to remove it (§5) |
| broken | 1e-1…10 | window identity change (swap/escape, §5) | yes — recognized and rejected by the diagnostics |

§5 maps every rung to a class.

---

## 3. The method hierarchy

Each method gets the same triple: why it exists, what it computes, how we
run it. Two vocabulary items are defined here and used everywhere after:
the *substrate* (§3.3) and the *basin* (§3.3).

### 3.1 ROHF — the reference that defines "core"

**The reference calculation fixes what the core is: everything downstream
folds it.** Two of its aspects matter here: how an orbital's *type* is
measured, and how the reference is kept honest on transition metals.

**An orbital's type is measured by projection.** For an AO span with
index set o, the S-metric projector is P_o = S[:,o] S_oo⁻¹ S[o,:], and an
orbital c carries weight w = c†P_oc ∈ [0,1] on that span (bounded — a
true projector weight; the same metric appears in §3.2's window
construction). Occupations can then be assigned by *character* rather
than energy order.

**On transition metals, energy-ordered filling fails twice, in both
cases measurably.** First, it runs to the wrong configuration: on Ti,
default ROHF converges to a doubly-occupied 3d inside the core with a
{4s, 4p} open block, +123.07 mEh above the physical 3d²4s² ³F state, and
every energy-ordered start tested falls back into that trap — the repair
is character-forced filling (occupy, each cycle, the orbitals with the
largest weights on the closed and open spans). Second, even with forced
occupations the SCF is multi-basin: a metastable +11.2 mEh solution with
the same occupations exists, whose polarized core smears later
constructions. The reference is therefore gated on its energy, with
*cycled initial guesses* — on missing the target energy, the SCF is
restarted from a different built-in starting density (PySCF guesses:
minao, atom, huckel, …) until the reference basin is reached. Mechanisms
in §5-L5.

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
  ensemble. Invariance holds self-consistently (the anisotropy never
  develops), but the reference itself differs — open-shell exchange
  against smeared rather than localized valence — by +8.7 mEh (C) to
  +38 mEh (Ti).
- **raw** — the core field as the SCF left it: ε unscreened, the
  multiplet splits by §2.4's mechanism.

With the averaged fold, the CASCI term tables of §7.1–7.2 come out exact
on all four atoms; the averaged–raw and averaged–frac differences are
the ΔE anchors of §7.3.

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
metric and the remainder re-normalized. The construction is forced by
convergence — a non-orthonormal substrate makes the optimizer diverge to
sub-variational energies (observed on Si; the §2.4 guard exists for
this) — and it has two measured prices: an energy offset (the projected
window is a slightly different, higher window definition; the ΔE columns
of §7.3 carry it at the ~88 mEh class on Sc/Ti) and a closure defect
(the projection against a non-scalar core span injects non-closed
directions — the input-window baseline, §5-L4a).

**Restrictions are masks on rotation pairs.** *Window frozen*: only
core↔virtual pairs — correlation through the core alone, closure
untouched. *Radial restriction*: only pairs within a same-(l, m) radial
chain — the symmetry-preserving channel; where it exists (C, Si) it
recovers the full gain at exact closure. *Free*: no mask — the textbook
run. **A basin is the region of attraction of one local minimum of this
optimization**; the free landscape has several (§7.3's chaos paragraph).

### 3.4 SA-CASSCF — averaging over the multiplet

**State-averaging replaces the single root by the multiplet average:**

  E_SA(orbitals) = (1/(2L+1)) Σ_a E_a(orbitals),  a over the M_L
  components,

minimized over rotations (³P: three components; ²D: five; ³F: seven).
Where a single root represents one M_L component — and optimizes for it
alone — the average is the rotationally scalar object.

**The average anchors the term identity but not the window.** Measured:
⟨L̂²⟩ is restored (Sc 6.002 where the single root reads 3.5), yet closure
is not guaranteed. Two reasons, both mechanisms of §3.6: along a flat
direction the average is stationary at the m-symmetric window (all
members of each l sharing one radial) but with **no restoring force** —
meaning zero curvature: the energy does not rise away from the symmetric
point in every direction, so nothing pushes the optimizer back to it,
and convergence tolerance, not physics, decides where it stops. The
remedy is to remove the freedom (freeze the window / restrict rotations),
not to tighten convergence. And on Sc/Ti the free state average can
still reach the swap basin. The symmetric endpoint combines the tools:
state average, window frozen, isotropic core substrate (§5.3).

**Is the symmetric window stationary for the average? Measured, not
assumed.** A member tilt (one window column rotated toward its next
radial, §3.6) is *not* a rotation of the atom — a rotation maps a closed
window onto itself; a tilt deforms it into a non-closed one, and the
term components split at first order. What the tilt does to the average
therefore depends on which physics the tilted radial carries. Measured
(frozen-window scans, §7.4: mean over the ground term's components, one
member tilted by δ = 0.1 from the raw shell): on C and Si the mean
*descends* (−28.1 and −11.0 mEh) — tilting toward 3p/4p is a fraction of
the beneficial radial rotation θ* for every component containing that
member; on Sc and Ti it *rises* (+1.6 and +4.9 mEh) — there the same
rotation only loses nuclear attraction. In both regimes the m-symmetric
radial direction is orders of magnitude steeper (the full −250.8 mEh on
C), which is why converged state-averaged runs land m-symmetric
(measured: ‖[P, L̂⊥]‖ ~ 1e-4, §7.3): the optimizer follows the symmetric
direction first, and on the d atoms the average turns the single-root
tilt descent — the Jahn–Teller-like instability of §7.4, for which the
symmetric window would be a saddle for one state — into an uphill
direction.

### 3.5 What textbook runs give you, and why

| Standard practice | What you get | Why | What we do instead |
|---|---|---|---|
| CASCI over canonical virtual MOs | closure defects at the 1e-1 level; term labels unreliable | canonical MOs diagonalize a Fock built from an anisotropic density; quasi-degenerate same-l members mix arbitrarily | raw shell windows, averaged fold (§3.2) |
| free CASSCF, default start | p²: full gain with a silently member-tilted window; d-shells: escape/swap basins, ⟨L̂²⟩ off by units | flat directions (endpoint = noise) and descending directions (§3.6); window identity free to change | restrictions + substrate control (§3.3) |
| free SA-CASSCF | ⟨L̂²⟩ ≈ L(L+1) restored, closure not (Sc: 6.002 with ‖[H,L̂z]‖ ≈ 0.8) | the average is scalar to first order only; no restoring force on flat directions; other basins still reachable | frozen-window SA on the isotropic substrate (§5.3) |

### 3.6 Why optimization breaks m — the variational surface

**The optimizer moves on an energy surface over window deformations;
§2.3's decomposition (l-content, member structure, radial content) reads
it.** Both symmetry-breaking mechanisms live in the member part. A tilt
is applied as

  |m̃(θ)⟩ = cos θ · |n l m⟩ + sin θ · |n′ l m⟩  (S-normalized),

which preserves l and m — the determinant sector and its M_L content are
unchanged; only one member's radial content differs. All other members,
the core and the diagonalization stay fixed; we read E₀(θ).

**The key structural fact: over m-labeled determinants, H is
block-diagonal in M_L.** The one-body fold is m-diagonal (a scalar
conserves m); the two-body Coulomb term is a scalar and conserves total
M_L. Each diagonal CI element depends only on the members its own
determinant occupies; couplings exist only between determinants of equal
M_L.

**The dodge, derived: one-dimensional M_L blocks make single-member
tilts exactly flat.** Carbon's ³P Ms=1 sector has three determinants
with three different M_L values — three 1×1 blocks:

  |M_L=+1⟩-block: |+1, 0⟩;  |M_L=0⟩-block: |+1,−1⟩;  |M_L=−1⟩-block:
  |−1, 0⟩  (m-labeled members),

so H is diagonal there, each determinant's energy depending only on its
own two members. Tilt any member: the determinant avoiding it is
untouched, and since all three are degenerate at θ = 0,

  E₀(θ_m) = min over blocks = const  for every θ_m

— exactly flat, measured to 1e-12. The optimizer may break m at zero
cost: nothing drives it to the m-symmetric window, nothing penalizes
leaving it, and its endpoint is numerical noise (§7.3's five-draw
census).

**The descent, derived: multi-dimensional M_L blocks open a coupling
channel.** In d-shells the sectors have M_L blocks of dimension ≥ 2 — a
d² M_L = +1 block, for instance, couples |+2,−1⟩ ↔ |+1, 0⟩
through the exchange integral ⟨+2,−1| r₁₂⁻¹ |+1, 0⟩ (Δm₁ = −1, Δm₂ = +1:
allowed, because the Coulomb term only fixes Δm₁+Δm₂ = 0). A tilt now
moves both diagonal elements and couplings within the block; on the
degenerate multiplet the deformation's block is traceless (§2.4), its
eigenvalues split symmetrically about zero, and the lowest descends
linearly:

  E₀(δ) = E₀ − |V(δ)| + O(δ²).

The CI ground state follows the descending level. Measured on Sc:
single-member tilts descend (−0.3 mEh at δ = 0.1), deeper with more
tilted members, with the ground ⟨L̂²⟩ drifting as neighboring terms join
in (2.00 → 1.87); Ti shows the same at half depth (§7.4). This
mechanism — a degenerate electronic state, a symmetry-lowering
coordinate, a guaranteed linear descent — is that of the Jahn–Teller
effect, where the coordinate is a nuclear displacement; here it is the
window's member asymmetry, hence the shorthand "JT-like" in the data
sections.

**Why the state average helps — and where it stops.** The first-order
shifts of the components sum to the trace of the deformation's block —
zero (traceless). So Σ_a ∂E_a/∂δ = 0: the average is stationary at the
m-symmetric window. But stationarity with zero curvature in some
directions (the dodge directions) is not a minimum with restoring force;
and a *different* basin, if reachable, is still found (the swap).
Symmetry is bought by the average only up to the flat directions;
exactness requires the constraints (§3.3) and the substrate (§5.3).

---

## 4. The pipeline tour

| Module | Owns |
|---|---|
| `system` | validates the system: one atom at the origin, spherical basis (L̂ is defined about r = 0) |
| `angular` | one-electron L̂x/L̂y/L̂z AO matrices, two independent constructions (libcint integrals; exact algebra over the real spherical harmonics) that pin each other |
| `detbasis` / `secondquant` | the many-body matrices over determinants: Slater–Condon rules and explicit a†/a — elementwise-identical twins (§2.2) |
| `operators` | operator assembly: the fold (§3.2), H/L̂z/L̂²/Ŝ²/N̂ symbols, and the window-closure certificates ‖[P, L̂z]‖, ‖[P, L̂⊥]‖ (§2.3) |
| `drivers` | the deep module — see the capability map below |
| `stages` | one uniform analysis pass per finished stage — the source of every number in §7 |
| `scans` | the frozen landscape measurements of §7.4 (deterministic) |
| `analysis` | term labels from ⟨L̂²⟩/⟨Ŝ²⟩, configuration summaries, report rendering |
| `pipeline` / `ladder` | orchestration: the per-atom CASCI sectors (`pipeline`) and the rung registry (`ladder` — treatments as data; §7.3's tables are its output) |

**The drivers capability map.** Every treatment in this project is a
point in a five-axis space:

| Capability | What it is |
|---|---|
| ROHF reference | character-forced filling + energy gate → the anchor energy and orbitals (§3.1, §5.5) |
| CASCI | exact sector diagonalization; fold variants raw / averaged / fractional (§3.2) |
| CASSCF | orbital optimization; freedom by mask, start by substrate (§3.3) |
| SA-CASSCF | multiplet-averaged optimization; weights over the M_L components (§3.4) |
| core construction | the purified, exactly isotropic core span (§5.3) |

The 17-rung registry tabulated in `results/` enumerates the corners of
this space; §7.3 keeps the rows that differ in mechanism.

---

## 5. Where symmetry breaks

| Level | Symptom (measured) | Mechanism | Fix | Cost |
|---|---|---|---|---|
| 1. Basis bookkeeping | labels mislead (2p\* hits 2s); shells fuse | wrong inputs → wrong constructions | label + metric discipline | none |
| 2. Sector witnesses | ‖[H,L̂²]‖ = 0 with defects present | L̂² = L(L+1)·1 on a single-irrep sector; ±pair blind spots | projector pair; multi-term sectors | sector size (3 → 120 dets) |
| 3. Core field | splitting 4.4e-5…4.1e-3 Eh; commutators 1e-5…1e-1 | polarization chain (§2.4) | fold cures / purification | ΔE per §7 |
| 4. Orbital optimization | tilted windows, floors, escape/swap | flat and descending directions (§3.6); identity change | restrictions + substrate | ΔE above deepest minimum: ~0 (C/Si), +24.3 mEh (Sc), +30.9 mEh (Ti) |
| 5. SCF reference | +123.07 mEh wrong configuration; +11.2 mEh metastable | energy-ordered filling; degenerate-open-shell basins | projector-weight filling + energy gate (§3.1) | per-atom reference metadata |

Class mapping (§2.5): *exact* — raw-window CASCI, frozen-window rungs on the
isotropic substrate; *baseline* — the L4a defect, the L4b floor; *broken* —
the L4c escape and swap.

### 5.1 Level 1 — basis bookkeeping

**Reading the basis correctly is a precondition, not physics.** AO labels
are truthful yet treacherous: pattern matching is regex-flavored; shells
fuse (one label, several contracted functions — unpack through
contraction data); same-l shells of different n overlap (⟨3p|2p-core⟩ ≈
0.3), so every projection and normalization runs in the AO overlap
metric (§3.1's P_o), never in the Euclidean one.

### 5.2 Level 2 — the sector and its witnesses

**A wrong window can pass a check the sector cannot testify against.**
On a single-irrep sector L̂² = L(L+1)·**1**, so ‖[H, L̂²]‖ ≡ 0 whatever
is wrong (carbon's tilted window: ⟨L̂²⟩ = 1.72 with ‖[H, L̂²]‖ = 0); and
L̂z cannot see m = 0-member or equal-±pair deformations (§2.3). The
universal certificate is the projector pair ‖[P, L̂z]‖, ‖[P, L̂⊥]‖; for
the sector commutators, prefer multi-term sectors (Ti's 120-determinant
sector catches at ‖[H, L̂²]‖ ≈ 1–2 what its ‖[H, L̂z]‖ = 1e-14 cannot
see).

### 5.3 Level 3 — the core field

**The reference leaves an anisotropic field: ε acts at first order on
the multiplet (§2.4).** Symptoms: raw-fold splittings 4.4e-5 Eh (C) to
4.1e-3 Eh (Ti). Two cures operate at the density level — one-shot
spherical averaging (ρ → ρ̄) and the fractional-occupation reference (a
different SCF, self-consistently invariant). **The third cure, core
purification, exists because of the window, not the core.** The window
an optimizer starts from is (1−P_core)-projected (§3.3), so its starting
closure is bounded by the *core's* closure — the L4a inequality below
has the core's nonclosure as a factor. Purification rebuilds the core
span to factorize as E_l ⊗ {all members of l} for each closed l —
§2.3's closure theorem applied to the core itself — which zeroes that
factor: the projection preserves closure, and the fold field becomes
scalar by construction of the orbitals. We clean the core because that
is the only way to clean the window's start. Energies and bookkeeping
in §7.

### 5.4 Level 4 — orbital optimization

**Three independent mechanisms, measured by two tools**: window-closure
certificates for (a), frozen scans for (b)–(c) — deform the window along
one coordinate and re-diagonalize the CASCI in it; no optimizer anywhere
(deterministic; §7.4).

**(a) The start-window defect.** The required orthogonalization (§3.3,
w̃ = (1−P_core)w/‖·‖_S) subtracts non-closed directions when the core
span is not rotation-closed:

  ‖[P_w̃, L̂]‖ ≲ ⟨window|core⟩ · ‖[P_core, L̂]‖,

with the window–core overlap large on transition metals (4s vs the 3p
core ≈ 0.3): Sc's input windows carry 5.7e-3/2.1e-2 defects that persist
into convergence; Ti's are reference-basin-dependent (3e-8 gated, 1.6e-2
metastable). The purified substrate zeroes the second factor — the same
rungs then close at the exact class.

**(b) Flat directions and the floor.** Where M_L blocks are
one-dimensional, E₀(θ_m) is constant in θ_m (§3.6) — *constant*, hence
all its derivatives vanish; nothing prefers the m-symmetric window.
Symptoms: member-tilted results with quiet many-body commutators (C:
‖[P, L̂⊥]‖ ≈ 1.1 at ‖[H, L̂z]‖ ≈ 5e-7); the state-averaged floor (Sc
4.5e-07, Ti 1.1e-04), which tightening the convergence tolerance does
not move (measured: conv_tol 1e-7 → 1e-12 changes nothing) — the
optimizer rests on a plateau, not in a well. The fix removes the
freedom: the radial restriction recovers the full gain at exact symmetry
on C/Si; the frozen window guarantees closure by construction. The
optimizer can also exceed any single coordinate (Ti: all raw
one-parameter rotations rise, the radial-masked rung still gains ~11 mEh
via the substrate's s-slot) — the scans map the landscape, they do not
enumerate it.

**(c) Window identity basins.** The optimizer may change what the window
*is*, measured by the span overlap σ = σ_min(C_a†SC_b) (σ = 1: same
span; σ < 1: identity change). Two basins on Sc/Ti: the escape (L ≈ 0,
window borrows half a p shell, 38–46 mEh below the frozen window) and
the swap — the window's 4s traded for a core 3p. The swap solution
deserves respect: it is variational, it still identifies the term
(⟨L̂²⟩ = 6.002), it is the lowest solution we find, and it exposes the
mechanism (an m-incomplete window is energetically reachable). Its
energy is not comparable to ours for one reason only: it has *redefined
the problem* — the core/active partition changed (a 3p core orbital
became active, the [Ar] core no longer exists as such), so a different
CAS is being solved, with core correlation legitimately included.
Wanting that correlation honestly means enlarging the window
deliberately, not letting the optimizer silently redraw the partition.
On a multi-basin landscape every free run is classified — not merely
converged — by σ, window weights and the projector pair.

**A rebuilt window is not a broken window — an m-incomplete one is.**
Closure (§2.3) requires complete member sets with shared radials per l,
whatever the orbitals' origin: a mean-field-relaxed radial span, and
even a complete multi-l window ({3d, 4s, 4p} on Sc — ~25 mEh of 4p
correlation at exact closure), are legitimately symmetric. What no
rebuilt window obtains for free is the optimizer's respect: closure is
imposed (masks, substrate) or certified (the projector pair), never
converged to — a complete start window still escapes (§7.3).

### 5.5 Level 5 — the SCF reference

**The bottom of the stack fails first on transition metals.** Energy-
ordered filling runs away (+123.07 mEh on Ti: a doubly-occupied 3d in
the core; every energy-ordered start lands there), and with occupations
fixed, the degenerate open shell still leaves a multi-basin SCF (+11.2
mEh metastable whose polarized core degrades everything downstream —
window defects 1.6e-2 against 3e-8 in the good basin). The repair:
filling by projector weight w = c†P_oc (§3.1) plus an energy gate
|E − E_ref| ≤ 1e-4 with cycled starts. **E_ref is the energy of the
correct reference basin, extracted once per atom by the diagnostics
themselves**: a first run is classified by character (occupations by
projector weight, core purity, window closure); the basin passing all
checks is recorded as the reference, and its energy becomes E_ref —
self-calibrating, no external data. The benefit propagates: the gated
reference is what puts Ti's window-frozen rung at 2.3e-08. The cost is
per-atom reference metadata (AO spans, E_ref) — data, not code.

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
data, §3.1/§5.5); (2) the CASCI anchor — term tables must come out exact
(§7.1–7.2); (3) the certified treatment (expect the exact class, §7.3's
pur rows); (4) both direct runs, classified; (5) optionally, frozen
scans for the radial/tilt anatomy (§7.4).

---

## 7. Results tables

All numbers below are measurements of this project (C, Si, Sc, Ti;
cc-pVDZ). The artifacts themselves live under `results/` (regeneration
commands in §8); every rung of every stage receives the same full-sector
analysis — term table of the whole sector, not only the ground root — so
the tables here are summaries of a uniform diagnostic pass. On the
symmetry-broken basins the term table itself scrambles (Ti's swapped
state-averaged rung counts ³F: 38 instead of 35, ³G: 25 instead of 18):
independent evidence beyond the commutators.

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

### 7.3 The ladders — one row per mechanism

The full 17-rung registry per atom (every orbital treatment, including
start-variants and single-root radial runs) is tabulated in
`results/*_ladder.md`. Here each atom is presented twice: a summary card
(totals against the ROHF reference, in eV, broken and symmetric variants
side by side), then the mechanism tables (ground energy in our fold; ΔE
in total energies against ROHF, mEh; ⟨L̂²⟩; many-body and window-span
commutators).
"Substrate" means the starting MO matrix the optimizer receives — core
block, window block, complement (§3.3).

**How to compare.** The reference is the ROHF total energy — one
unambiguous number per atom, the starting point of everything. Two rules
remain. *Fold signs*: the E_ground column is our fold; µEh…mEh signed
differences inside the window-frozen family are fold bookkeeping (the
max|dE| column of the full tables is the detector). *Partitions*: rows
that moved the core/active boundary (the swap) are compared in totals
only — in totals the swap SA-CASSCF reads −13.3 mEh (Sc) and −10.1 (Ti)
below ROHF, not the −730/−786 mEh our fold seems to show.

**The summary cards.**

**The reference table** (all four atoms; totals in Eh; the anchor every
ΔE hangs on. The averaged fold leaves the ROHF total unchanged — measured
raw-vs-avg ground shifts: −29 µEh / −0.31 mEh / −1.5 mEh / −2.8 mEh, the
multiplet-splitting scale. The frac-occ SCF is a different functional,
not a variational neighbor of ROHF):

| reference | C | Si | Sc | Ti |
|---|---|---|---|---|
| ROHF (raw core = averaged fold) | −37.682418 | −288.846437 | −759.736255 | −848.406489 |
| frac-occ SCF (ΔE vs ROHF, eV) | −37.678694 (+0.10) | −288.867793 (−0.58) | −758.014967 (+46.8) | −847.296157 (+30.2) |

**C / Si rungs** (ΔE vs ROHF, eV; CAS(2,3), Ms=0):

| rung | C | Si |
|---|---|---|
| CASCI raw shell, averaged fold | +6.82 (confinement; exact) | +6.03 (confinement; exact) |
| CASCI raw core fold | +6.82 (comm 9.8e-6 — the detector) | +6.03 (comm 5.6e-4) |
| CASCI frac reference | +6.59 | +5.93 |
| CASSCF free | 0.00 (member-tilted window, baseline sector) | 0.00 (lands symmetric) |
| CASSCF radial restriction | 0.00 (full recovery, exact) | +0.001 (exact) |
| SA-CASSCF | +0.001 (lands symmetric, exact) | +0.001 (exact) |

**Sc** — CAS(3,6)/{3d,4s}

| deliverable | free variant | symmetric variant |
|---|---|---|
| CASCI, raw shell, averaged fold | — | +0.18 (exact) |
| CASCI, raw core fold | — | +0.18 (raw split: commutator 2.8e-2 — the polarized-core detector) |
| CASCI, frac reference | — | +1.00 (invariant; a different reference) |
| CASSCF | **−0.86** (escape: L ≈ 0, identity changed — broken) | +0.18 frozen (baseline) / −0.18 radial (single root, gauge-chaotic) |
| SA-CASSCF | −0.36 (swap; ⟨L̂²⟩ = 6.002, term-readable — broken) | +0.20 frozen purified (exact) / 0.00 radial purified (floor) |

**Ti** — CAS(4,6)/{3d,4s}

| deliverable | free variant | symmetric variant |
|---|---|---|
| CASCI, raw shell, averaged fold | — | +0.29 (exact) |
| CASCI, raw core fold | — | +0.29 (raw split: commutator 4.9e-2) |
| CASCI, frac reference | — | +1.31 (invariant; a different reference) |
| CASSCF | **−0.85** (escape — broken) | +0.28 frozen (baseline, gated basin) / −0.06 radial (single root, gauge-chaotic) |
| SA-CASSCF | −0.27 (swap — broken) | +0.03 radial purified (floor) / +0.32 frozen purified (exact) |

One sentence per atom: where the raw-shell CASCI total sits against
ROHF (+6.82 / +6.03 eV on the p² atoms — there it is the confinement
cost; +0.18 / +0.29 eV as run on the d atoms, where the window already
correlates — our raw-span analysis fold reads −2.2 / −2.0 eV there, the
window-definition offset), what breaking buys (−0.86 / −0.85, or the
swap's −0.36 / −0.27), and what symmetry costs against ROHF (≈ 0 on p²;
−0.02 to +0.03 on the d-atoms for the best symmetric treatment).

None of these numbers is a fine paid for symmetry. The p² CASCI row is
confinement, recovered in full by an m-symmetric rotation; the escapes
are window-identity changes — a different question, not a cheaper
answer to this one. What symmetry itself costs is ≈ 0 eV on the p²
atoms and −0.02 to +0.03 eV on the d atoms — and it buys term tables
exact to machine precision with correlation included.

The escape and the swap fail differently. The escape rebuilds the
active space from orbitals the configuration never implied; its states
are not L̂² eigenstates at all (⟨L̂²⟩ wanders across draws) — it
describes no term. The swap keeps the window's size but changes its
identity (a valence 4s for a core 3p); its ground root keeps a readable
label only approximately — the commutator is O(1) and the term table
scrambles (Ti: ³F 38 instead of 35, ³G 25 instead of 18). The terms are
described exactly by the certified treatment; the direct runs are
reported for their energies and their anatomy, not their labels.

**The confinement finding.** The +250.8 on carbon is *purely* the
raw-shell confinement, verified directly: CASCI over the ROHF's own
2p-like orbitals (canonical window) reproduces the ROHF total exactly
(ΔE = 0.00 mEh) — the ³P pair carries no net correlation the CAS can
express, and the famous "−250.8 mEh CASSCF gain" is the recovery of the
confinement cost, not correlation won. The optimized window (~2/3 2p +
1/3 3p at θ\* = 0.61, §7.4) simply rebuilds the relaxed ROHF-like
orbital inside the window. The same holds on Si (+221.5 recovered); on
the d-atoms the as-run confinement reading is small (+6.7 / +10.6 mEh)
and the state-averaged radial channel nets −0.8 (Sc) / −0.4 (Ti) mEh —
the scan verdict of §7.4 in totals.

**Chaos.** Rungs marked † have an optimizer landscape with several
distinct local minima, each with its own region of attraction (a
*basin*); which one a run reaches is decided by last-digit SCF noise.
Two cases are measured. On a *flat direction* the noise picks the
stopping point: five fresh draws of carbon's single-root radial rung
give one common energy (spread < 1e-6 mEh) while ⟨L̂²⟩ ranges over
1.69–2.00 — one draw landed the symmetric window itself. On a
*multi-basin* landscape the noise picks the basin: the same census on Sc
spreads 2.2 mEh with ⟨L̂²⟩ 3.5–5.7.

Carbon's free CASSCF is the instructive case: the full −250.8 mEh
window gain (the confinement recovery, above) with ⟨L̂²⟩ = 1.69–1.72
while ‖[H, L̂²]‖ = 0. Both are
window effects on this smallest sector. The optimized window is
member-tilted (‖[P, L̂⊥]‖ ≈ 1.1); each tilted member is still an l=1,
m-eigenfunction, so the one-electron parts of L̂² are unchanged, but the
ladder cross-terms are weighted by radial overlaps between members with
different radials — the window-projected L̂² is no longer 2·**1**
(eigenvalues 1.723, 1.862, 1.862, against exactly {2, 2, 2} on the raw
shell), and the ground expectation drifts to 1.723. The vanishing
commutator is the single-irrep caveat of §2.3: with one irrep on three
determinants, H and the window-L̂² are simultaneously diagonal (verified
— the expectations in the H eigenstates reproduce the L̂² eigenvalues
exactly), so ‖[H, L̂²]‖ carries no closure information here. (C's full
mechanism table, like Si's, is in `results/`; its summary card above
carries the essentials.)

**Si:** analogous to carbon throughout — raw 5.6e-04; the +221.5 mEh
confinement recovered by every optimizer rung; the pur family at the
exact class. Full table in `results/Si_ccpvdz_s2_ladder.md`.

**Sc — CAS(3,6), mechanism table:**

| rung | E_ground | ΔE vs ROHF (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | annotation |
|---|---|---|---|---|---|---|---|---|
| rohf | −759.736255 | 0 | — | — | — | — | — | the anchor: character-forced filling (§5.5); raw core density anisotropy 1.3e-2 |
| frac-scf | −758.014967 | +1721 | — | — | — | — | — | the smeared-shell ensemble: a different functional (§3.2); only its isotropic core field is used |
| casci_avg | −759.817560 | +6.7 | 6.0 | 5.6e-15 | 6.4e-13 | 7e-15 | 9e-15 | spherically averaged core: reference CASCI (§3.2) |
| casci_raw | −759.819083 | +6.7 | 6.0004 | 2.8e-02 | 1.4e-01 | 7e-15 | 9e-15 | raw field: percent-level quadrupole (§2.4) |
| casci_frac | −759.787295 | +36.8 | 6.0 | 3.4e-14 | 3.5e-13 | 7e-15 | 9e-15 | average-of-configuration: a different reference (§3.2) |
| CASSCF (free) † | −759.765900 | **−31.5** | 0.05 | 3.2e-01 | 7.2e-01 | 2.3 | 2.7 | unconstrained single root: escapes to the L≈0 basin — the deepest minimum, identity-destroying (§5-L4c) |
| casscf_relax | −759.728374 | +6.6 | 5.9998 | 2.8e-03 | 4.0e-02 | 5.7e-03 | 2.1e-02 | window frozen, core relaxes; carries the start-window defect (§5-L4a) |
| sa † | −760.548082 | −13.3 | 6.002 | 8.3e-01 | 3.8 | 1.4 | 1.4 | SA-CASSCF: the 4s↔3p swap basin; the fold ΔE is partition bookkeeping (§5-L4c) |
| sa_radial | −759.735593 | −0.8 | 5.9997 | 3.3e-03 | 5.4e-02 | 6.0e-03 | 2.2e-02 | SA-CASSCF + radial restriction: at the input-window baseline (§5-L4a) |
| pur family | −759.7284 … −759.7362 | +6.7 frozen / −5.8 radial / +0.0 SA-radial | 6.0 | 3.5e-14 (frozen) / 4.5e-07 (SA-radial floor) | ≤4.8e-13 | 7e-15 | 9e-15 | isotropic core substrate (§5.3): the input defect gone; the swap basin survives (one draw: `sa_pur`, 9.2e-02) |

**Ti — CAS(4,6), mechanism table:**

| rung | E_ground | ΔE vs ROHF (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | annotation |
|---|---|---|---|---|---|---|---|---|
| rohf | −848.406489 | 0 | — | — | — | — | — | the anchor: character-forced filling + E_ref basin gate (§5.5) |
| frac-scf | −847.296157 | +1110 | — | — | — | — | — | the smeared-shell ensemble: a different functional (§3.2); only its isotropic core field is used |
| casci_avg | −848.481312 | +10.6 | 12.0 | 2.2e-14 | 4.1e-13 | 6e-15 | 8e-15 | spherically averaged core: reference CASCI (§3.2) |
| casci_raw | −848.484093 | +10.6 | 12.0000 | 4.9e-02 | 2.6e-01 | 6e-15 | 8e-15 | raw field: the quadrupole at first order (§2.4) |
| casci_frac | −848.443405 | +48.0 | 12.0 | 1.2e-13 | 7.5e-13 | 6e-15 | 8e-15 | average-of-configuration: a different reference (§3.2) |
| CASSCF (free) † | −848.434973 | **−31.3** | 0.05 | 4.3e-01 | 1.1e+00 | 2.8 | 4.0 | escape to the L≈0 basin — the deepest minimum, identity-destroying (§5-L4c) |
| casscf_relax | −848.393394 | +10.4 | 12.0 | 2.3e-08 | 1.7e-07 | 3.2e-08 | 5.9e-08 | window frozen, core relaxes; the gated reference basin keeps the start defect small (§5-L4a, §5.5) |
| sa † | −849.267643 | −10.1 | 12.003 | 1.1e+00 | 6.0e+00 | 1.4 | 1.5 | SA-CASSCF: the exact 4s↔3p swap (window p = 1/6); fold ΔE is partition bookkeeping (§5-L4c) |
| sa_radial | −848.404105 | −0.4 | 12.0 | 8.7e-04 | 3.5e-03 | 1.2e-03 | 7.4e-04 | SA-CASSCF + radial restriction: the substrate s-slot gain, at the input baseline |
| pur family | −848.3934 … −848.4099 | +10.4 frozen / −3.5 radial / +1.0 SA-radial floor | 12.0 | 1.0e-13 (frozen) / 1.1e-04 (floor) | ≤1.1e-12 | 6e-15 | 8e-15 | isotropic core substrate (§5.3); sa_pur = a partial swap (s-slot ⅔ s + ⅓ p, −3.8) |

### 7.4 The scan anatomy — where each atom's gain really lives

Method and motivation in §5-L4; in brief: freeze everything, deform the
window along one coordinate — the m-symmetric radial rotation
(closure-preserving) or a member tilt (closure-breaking) — and
re-diagonalize; no optimizer anywhere, so the numbers are deterministic.
Member sets below use the §2.3 ordering (dxy, dyz, dz², dxz, dx²−y² =
0–4).

| Atom | m-symmetric radial rotation | single-member tilt | where the real CASSCF gain lives |
|---|---|---|---|
| C | recovers everything: θ\*=0.61, −250.8 mEh of window energy — the confinement cost of §7.3 | exactly flat (five draws, one energy) | the θ\* rotation itself |
| Si | gains everything: θ\*=0.85, −126.7 mEh | exactly flat | same |
| Sc | rises | descends: −0.3 mEh @ δ=0.1; deeper with more members (−0.9 with four) | nowhere radial; the ceiling sits +25 mEh above the escape basin |
| Ti | rises (+9.8 mEh @ θ=0.05) | descends: −163 µEh @ δ=0.1; the pair {0,1} doubles it; {0,1,2} and {0,1,3,4} rise | the orthogonalized substrate's s-slot (4s→6s admixture), ~11 mEh |

Reading: *flat* → the optimizer may break m at zero cost (§2.3, flat
directions); *descending* → it is driven to break m (mechanism §5-L4);
*rises* → the direction buys nothing. Sc and Ti share the qualitative
structure but differ in three measured ways: tilt depth (Sc ≈ 2× Ti per
member at δ = 0.1); multi-member patterns (Sc descends further with more
tilted members; Ti's {0,1} pair doubles the single-member descent while
{0,1,2} and {0,1,3,4} rise steeply); and the substrate channel (Sc has
no radial gain anywhere; Ti's radial-masked rung gains ~11 mEh through
the orthogonalized s-slot — visible as the sa_radial-minus-relax
difference in §7.3).

**What the prices say.** On p² atoms symmetry is free: the radial
restriction recovers the entire window gain (−250.8 mEh on C — the
confinement cost) at exact closure, and state-averaging costs 0.00 mEh.
On d-shells it has a price: against ROHF, the best symmetric treatments
sit at −0.8 (Sc) to −0.4 (Ti) mEh while the escape basin reaches −31.5 /
−31.3 — a real, identity-destroying gain that no symmetric window can
express on these atoms (the table above), and the swap follows at
−13.3 / −10.1, reachable and term-readable rather than deeper.
Purification is the cheapest real repair of all: the start-window
baseline falls to the exact class at −0.05 mEh.

---

## 8. Dictionary and file map

**Primary sources of §7** (regenerated, not committed): per-atom
`results/*_report.md/json` — the pipeline CASCI multiplet tables
(§7.1–7.2); `results/*_ladder.md/json` — the full 17-rung registry
(§7.3's source, including the per-runk totals used in the cards).
**Provenance artifacts** — kept and regenerable, but investigation
records: `results/*_scan.md` (the source of §7.4's curves),
`results/*_sa_decomposition.md` (the source of the §5-L4 window-vs-core
attributions).

**Regeneration, per table.** §7.1–7.2: `python scripts/run_c.py` (and
the analogous `run_atom` calls per atom, recorded in each report
header). §7.3 cards and tables: the four `run_ladder` calls (recorded in
each ladder header; the canonical-window confinement check is
`run_casci(mf, ncas, sector, mf.mo_coeff, …)`). §7.4:
`python scripts/radial_scan.py` and the scan calls recorded in
`results/*_scan.md`. The executable
documentation is the test suite (91 tests): the term tables, the
exact-class closure pins, the basin branches. Plain-language companions:
`RUNGS.md` (the rung ladder) and the glossary appendix (vocabulary).

---

## 9. Roadmap

*Scalar relativity* (DKH/X2C): h₁ changes, nothing else — a scalar
Hamiltonian is rotationally invariant, so every theorem and diagnostic
survives; the reference basins (E_ref, the good-basin census) must be
re-measured once. *ECP-based full relativity*: the pseudopotential's
angular projectors are l-resolved, so windows stay closed, but the
AO-label and member bookkeeping (§5.1) must be re-verified per
pseudopotential. *Spin–orbit coupling*: the natural next symmetry step —
it splits the multiplets this note certifies; the machinery extends in
place (Ĵ² = (L̂+Ŝ)², Ĵz as determinant-basis matrices; term symbols
become levels). *Further atoms* (V, Cr, …): data rows, not code — the
L4/L5 phenomena are expected to recur, and the protocol card's checks
are the guard.

---

## Appendix. Plain-language glossary

Vocabulary for AtomTerm. Every term the reports and code use is defined
here in plain language; the precise statement follows the plain one and
lives in the sections of this note.

### The objects

- **atom config** — the input card: which atom, which basis, which spin
  (2S), which charge. Everything else is derived.
- **reference** — a solved mean-field (an SCF). The ROHF reference is
  *lopsided*: it occupied one specific open-shell orbital (say px), so
  its density remembers that choice. The fractional-occupation reference
  puts equal fractions of an electron in every open-shell orbital and is
  round (spherical) *by construction*.
- **core** — the orbitals we do NOT correlate (everything below the
  window). The core enters the active problem as an electric field: the
  "fold".
- **fold** — building the effective Hamiltonian for the active electrons
  with the core's field included (one-body part) instead of treating
  core electrons explicitly.
- **window** — the active space: the few orbitals whose electrons we
  correlate. The natural choice for an atom: the raw functions of one
  valence shell (e.g. the three 2p). A window made of *complete shells*
  is rotation-closed, which is what makes the angular-momentum analysis
  exact.
- **sector** — the determinant space for fixed numbers of alpha and beta
  electrons; the label is Ms = (na − nb)/2. Atoms need several sectors
  to hold all terms of a configuration.

### The analysis

- **member** — the m-component of a shell: px, py, pz are the three
  members of p. Under rotation they mix into each other.
- **m-completeness** — the real closure condition: a window is safe under
  rotations when, for each l it contains, it holds ALL members and those
  members share the SAME radial mixture (all built from the same
  combination of 2p/3p/…). l-purity alone is not enough.
- **commutator norms** ‖[H, Lz]‖, ‖[H, L²]‖ — the symmetry meters. Exactly
  zero when core+window is rotation-closed; their size is a measured
  number for every calculation ("rung") we do.
- **one-particle closure meter** ‖[P, Lz]‖, ‖[P, L⊥]‖ — the same question
  asked of the window *span* alone (P = its projector): no core, no CI.
  [P, Lz] is exactly the part the many-body [H, Lz] can see (blind to
  m=0-member and equal ±-pair tilts); [P, Lx/Ly] sees *every* member
  asymmetry. Both zero = the complete closure certificate, no blind
  spots.
- **span overlap** — smallest principal cosine between two spans (1.0 =
  identical). The *identity* gauge of the ladder ("did the window become
  a different space?"). Independent of closure: an overlap of 0.9855 can
  hide a member tilt of ~0.17 rad, worth closure defects O(0.5).
- **window purity** — how much of each window orbital lies in its
  nominal l (1.000 = no leakage). Necessary, not sufficient: a pure-p
  window can still be m-incomplete (see above).
- **stage / stage facts** — a finished calculation (CASCI, CASSCF, …) is
  a *stage*; `stage_facts` performs the single analysis pass every stage
  gets: fold, classify terms, compare energies with PySCF, check purity
  and orthonormality.
- **ortho_err / DIVERGED** — max |C†SC − I| of a rung's orbitals. Above
  1e-6 the optimizer lost orthonormality and every number of that rung
  is meaningless; the report says DIVERGED.
- **max|dE|** — our energies vs PySCF's own for the same calculation.
  Zero means identical bookkeeping; known physical offsets (e.g. we
  average a polarized core, PySCF uses it raw) show up here on purpose.

### The experiment axis

- **ladder** — the comparison series: the same atom and window analyzed
  under every orbital treatment, from frozen to fully optimized.
- **rung** — one step of the ladder (see RUNGS.md for what each rung
  holds fixed and releases).
- **rotation mask** — a rule declaring which pairs of orbitals the
  optimizer may rotate. Masks are how we impose physics questions on
  PySCF's optimizer ("only radial mixing", "no core exchange", …).
- **radial constraint** — the mask "same (l, member) only": 2px may mix
  with 3px, never with py or any s/d. Pure radial breathing.
- **gauge** — a choice of description that does not change the physics
  (like a phase). For us: *which* linear combination of the members each
  window orbital is — the span is the same, the energy is the same, only
  the description differs.
- **flat direction / flat member gauge** — a deformation of the window
  along which the energy does not change (measured < 1e-6 Eh on C/Si).
  Giving each member a *different* radial mixture is such a direction:
  it breaks the commutators (first order) but not the energy (second
  order at a stationary point). Optimizers wander along it for free —
  which is why free CASSCF "breaks symmetry" without gaining energy.
- **escape / escape basin** — on Sc the *free* optimizer rebuilds the
  window from ~half p-character (borrowed from virtual 4p and core 3p):
  the active space is no longer the 3d/4s shell at all. The "basin" is
  that solution's region of attraction; it is a real local minimum
  (variational, orthonormal), just not the window we intended. Not
  reachable by radial mixing (the radial ceiling sits ~25 mEh above it).
- **input-window baseline** — the ~1e-3 commutator level that every Si/Sc
  optimized rung inherits from its STARTING window: to make it
  orthonormal we must orthogonalize the shell AOs against the core, and
  the canonical core carries small asymmetric tails, so the projected
  window leaks under rotations before any optimization runs (§5.3).
  Not core polarization; not the optimizer. Fix = the closed core span.
- **closure repair** — (designed, then dropped) post-processing that
  would re-share the radial mixture equally among members, restoring
  exact closure at ~µEh energy cost, because that direction is flat.
  Dropped: on the isotropic core substrate the free state-averaged rungs
  close exactly by themselves (§7.3) — nothing left to repair.
