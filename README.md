# AtomTerm

Atomic term-symbol analysis on top of [PySCF](https://pyscf.org): the
pipeline runs ROHF → CASCI/CASSCF/SA-CASSCF, then classifies every
state by **exact** angular-momentum operators (L̂z, L̂², Ŝz, Ŝ² built
as plain matrices in the Slater-determinant basis) and certifies
rotational-symmetry preservation by ‖[Ĥ, L̂z]‖_F and ‖[Ĥ, L̂²]‖_F.
Atoms covered: **C** (2p²), **Si** (3p²), **Sc** (3d¹4s²), **Ti**
(3d²4s²), the uncontracted twins **Sc-unc**/**Ti-unc**, and **U14e**
(f³d¹; ECP) — the uranium ground term's parent ⁵L certified at
machine-precision closure.

## Quick start

```bash
conda activate pyscf    # PySCF 2.13; bare shell: source <conda-root>/etc/profile.d/conda.sh first
pytest -m fast          # 282 tests, seconds — the per-commit gate's base
python scripts/run.py C # term tables -> elements/C/
python scripts/regen_cells.py Ti --method sa-casscf   # one slice of the methods grid
```

Run from the repo root — `pyproject.toml` preconfigures
`pythonpath = src`; no install step. Everything human-readable lands
committed in `elements/<Atom>/`.

## What it is

The exact Hamiltonian of an atom is SO(3) × SU(2) invariant, so every
term (2S+1)L is exactly (2L+1)(2S+1)-fold degenerate. A CAS run breaks
this in exactly two places, and the whole design answers them:

- **The folded core field.** CAS takes its core from a prior SCF and
  represents the omitted electrons by f_tu = h_tu +
  Σ_core[2(tu|ii) − (ti|iu)] — the folded part carries whatever
  symmetry the core density has. A polarized core makes Ĥ_CAS merely
  axial; the cure is a symmetric core field (`core.py`).
- **The window span.** The projected algebra closes iff the active
  space is rotation-closed ([P, L̂] = 0) — the condition is
  **m-completeness**, imposed by construction, never converged to
  (`windows.py`).

The discipline: every window is built core-orthogonal and
rotation-closed, and the commutators + ⟨L̂²⟩ *measure* the result —
construction and verification are separate instruments.

## Methods

| method | what it does | what it answers |
|---|---|---|
| **ROHF / UHF** (the reference) | character-forced occupations (S-projector weight on declared spans, never energy order), basin gate, cached trusted rungs | the anchor — ΔE's zero; the filling's honesty |
| **CASCI** | exact diagonalization of one M_S sector at frozen core + window | the pricing class — what PT2 sees as built |
| **CASSCF** | one root, orbitals free | the escape disclosure — a single root breaks symmetry by m-incompleteness |
| **SA-CASSCF** | the ensemble over the declared term's M_L set; the two-pass readout in the frozen window | the term — the rotationally scalar object; exact ⟨L̂²⟩ |

Two routes run the same methods: **free** (symmetry-blind; ⟨L̂²⟩
*measured* afterwards) and **pinned** (the declared term's D2h
branching *constructed*: one `wfnsym` solver per component — the
frozen pinned CI reads term identity at 1e-13 and is the
term-identity certification instrument). **PT2** (strongly contracted
NEVPT2, via `pyscf.mrpt`) adds the dynamic-correlation axis: per root,
per Koopmans class.

## Core and window

**The core** is N chemical shell pairs (declared, never guessed from
labels). Three frozen fields are compared: `raw` (the SCF's density —
the anisotropy detector), `avg` (its l = 0 part), `frac` (a separate
average-of-configuration SCF); two substrates: `pol` (the SCF's own
core span) and `pur` (the rebuilt rotation-closed span, machine
precision). Every window is built orthogonal to the core span.
Definitions and construction: `core.py`, `cells.py`.

**The window** is m-complete by construction. The radial sources:
the *incumbent* (the reference's occupied-span Gram), `@r<k>` (the head
configuration's density blocks, capture-maximizing), `@s<k>` (the
blocks averaged over a configuration series), and the *mixed source*
(a closed pair takes the reference's own MO, sector-projected; open
shells take density-block radials). The device per shell: the AO
chain → the density block → the generalized eigenproblem (block,
S_r⁻¹) → m-replication. Definitions and construction: `windows.py`
(`window_distilled`), `basis.py` (`series_window`, `metric_radials`).
The `MO`/`MOC` starts take the SCF columns directly — m-incomplete by
construction, kept as certified controls.

The prescription the record states: **per-channel sources by the
channel's character** — a closed pair takes the reference's own MO
(sector-projected); open shells take density-block radials rebuilt
m-replicated (a raw-ROHF MO span never rotation-closes).

## What we measure

- **Closure certificates** — ‖[P, L̂z]‖/‖[P, L̂⊥]‖ certify the window
  span; ‖[Ĥ, L̂z]‖/‖[Ĥ, L̂²]‖ certify the assembled treatment (the
  ψ-residual form past the dense-meter cap). Zero certifies; a value
  localizes the breaking.
- **Term identity** — ⟨L̂²⟩, ⟨Ŝ²⟩ per root; the multiplet splittings;
  the component spread (ptp).
- **Spin honesty** (UHF) — ⟨S²⟩ by the overlap formula + the
  corresponding-orbital pairing spectrum.
- **PT2** — per root, per Koopmans class (norm, energy), on
  certified solutions only; a class norm far from its peers is the
  meter of a mispriced object.
- **The invariance structure** — E and ⟨L̂²⟩ are orbital-gauge
  invariant (span/Casimir); ⟨L̂z⟩ is not, and E(2) is state-specific
  rather than an eigenvalue: compare it only between certified
  solutions of the same term (redraws of a closed term agree to
  ≤ 0.003 mEh; a quasi-degenerate manifold's root partition moves it
  far more). Compare energies at matched partitions.

The worked example (Ti, the 3d4s window; the full table:
`elements/Ti/TABLE.md`): the free CASSCF escapes the declared window
(⟨L²⟩ collapses to 0.05) while the SA restores the term (12.0029) —
the pipeline **certifies** the atomic ground term, not just labels it.

## Peculiarities (measured, generic)

- A raw core field splits multiplets axially; the averaged core cures
  at sub-mEh cost.
- A single-root CASSCF breaks symmetry by m-incompleteness — its
  converged windows can lose rotation closure entirely; the SA
  restores the term.
- The MO/MOC windows read ⟨L²⟩ ≈ 0 by construction — the controls
  that prove the disease.
- The optimizer can escape or swap window identity; the partition
  bookkeeping catches it.
- On quasi-degenerate manifolds the class partition is window- and
  core-sensitive — the window decides the frozen class; the ensemble
  decides the optimized one.
- An E(2) that dwarfs its neighbours is first an instrument suspect:
  price only a certified object — the orbitals and the CI vectors of
  one solve.  The class norms are the meter that catches it.
- The static/dynamic correlation trade converts *within* a window
  (the frozen → the optimized): the static gain rises while the E(2)
  remainder falls; across windows the E(2) is nearly constant while the
  static gain carries the window economics — the trade is not an
  across-window law. See the V&V record (`scripts/u14e/VV.md`).
- A band-named shell is its dominant radial, not the whole shell; the
  core column names the core's shells when the noble guess would lie.

The uranium record carrying all of this measured and classified:
`scripts/u14e/VV.md` (the V&V family) and
`elements/U14e/TABLE.md` (the campaign methods table).

Reading `scripts/u14e/VV.md` (its row store is `scripts/u14e/vv.json`;
a key outside the declared table is refused):

- **Method column** — `free@pur`: PySCF-direct, started on the
  purified core; `pinned@pur`: the wfnsym-pinned state average, started
  on the same core (PySCF's optimizer lets the core relax; every
  solution's certificate reads how far it moved). One optimizer class
  per block; inside a block the plain windows come first, then the
  `@r1` windows.
- **Rows** — every dE is against the ROHF row. A free row carries
  root 0; a pinned SA row carries the mean over the term's component
  blocks.
- **conv** — `<achieved |dE|>/<macro-cycle count>`; a trailing `*`
  marks a cycle-cap draw; `…/—` a CI solve (no iteration structure).
- **E(2)** — strongly contracted NEVPT2 on the row's own solution:
  the orbitals and the CI vectors of one solve.
- **Dashes** in ⟨L²⟩, E(2) or E_tot — a cell withdrawn pending
  re-measurement; the row's other cells stand.

## Analysis (PT2)

**Why E(2) is state-specific.** The variational energy is an
eigenvalue — identical for every gauge. E(2) is a functional of one
state's density matrices, and inside a degenerate term block which
state a root IS depends on the partition. Its invariance is a
measured property, never an assumption: redraws of a closed term
agree to ≤ 0.003 mEh.

**The invariance ledger.** Invariant: the block projector, E_CAS,
⟨L²⟩, ⟨S²⟩. Not invariant: ⟨Lz⟩, the density matrices, and E(2) per
root. Spin is safe by construction (the Dyall partition is
spin-free).

**Reading the ladder.** *R* = E(2) priced on a window against the
dynamic correction the optimized solution of the same span carries: a
window whose R dwarfs its neighbors is *fueling*; a ladder whose
corrected totals are flat has SATURATED. The class norms are the
honesty meter of the OBJECT, not of a gauge: a one-hole (Sr) norm far
above its peers means the CI vector was not solved in the orbitals it
is priced on — a certified solution reads ~0.004 there, a mispaired
one 14–15. Absolute eV values are basis-bound; ladder *differences*
are the usable signal.

## Files and regeneration

| Module | Owns |
|---|---|
| `system` | validation — one atom at the origin, spherical basis |
| `angular` | one-electron L̂ AO matrices — two independent constructions |
| `detbasis` / `secondquant` | many-body matrices over determinants — elementwise twins |
| `operators` | operator assembly, the closure certificates, both PySCF convention boundaries |
| `reference` / `casscf` | the engines — gated SCF; the one optimizer driver with masks |
| `windows` / `core` | window construction + the kind vocabulary; the purified core substrate |
| `natural` | NO mathematics |
| `stages` | the one analysis pass per finished stage — the source of every number |
| `scans` | frozen-landscape measurements — deterministic |
| `analysis` | term labels, configuration summaries, report rendering |
| `cases` / `cells` | the declaration (`CaseSpec` → `run_case`); the cell addressing layer |

Everything human-readable lives in `elements/<Atom>/`: `TABLE.md` +
`TABLE.csv`, `report.md/json` (term tables), `cells.json/md` (the
treatment grid), plus frozen one-shot artifacts. New atom = new element
file, not README growth. Regeneration: term tables — `scripts/run.py`;
cells — `scripts/regen_cells.py`; the methods table —
`scripts/build_methods_tables.py`; scans — `scripts/radial_scan.py`.
The U14e V&V family — the library `scripts/u14e/vv` (the reference
substrate, the active spaces, the routes, the pricing, the table; each
module's docstring owns its contract), its legs
`scripts/u14e/u14e_vv_ladder.py` (the table's one writer: every case
of the CAS ladder through the free trio and the pinned SA) and
`u14e_vv_uhf.py` (the reference rows), the driver
`u14e_vv_orchestration.py`, the HF-virtual channel screen
`u14e_vv_screen.py` (print-only); its table `scripts/u14e/VV.md`. The
test suite is the executable documentation.

## Status and directions

C, Si, Sc, Ti, the uncontracted twins, and U14e delivered end-to-end;
the uranium ground term's parent ⁵L certified at machine-precision
closure, with the term-identity readout (the pinned CI) and the V&V
family re-deriving every number standalone. Open physics: a NIST
comparison round; inner core-pair valence. Forward: spin–orbit
(extends in place — Ĵ², Ĵz as determinant-basis matrices), transition
moments, g-factors, scalar relativity. Further atoms = element files +
data rows.
