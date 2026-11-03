# AtomTerm

Atomic term-symbol analysis on top of [PySCF](https://pyscf.org): the
pipeline runs ROHF → CASCI/CASSCF/SA-CASSCF, then classifies every
state by **exact** angular-momentum operators (L̂z, L̂², Ŝz, Ŝ² built
as plain matrices in the Slater-determinant basis) and certifies
rotational-symmetry preservation by ‖[Ĥ, L̂z]‖_F and ‖[Ĥ, L̂²]‖_F.
Atoms covered: **C** (2p²), **Si** (3p²), **Sc** (3d¹4s²), **Ti**
(3d²4s²), **U** (f³d¹; ECP).

## Quick start

```bash
conda activate pyscf    # PySCF 2.13; bare shell: source <conda-root>/etc/profile.d/conda.sh first
pytest                  # 185 tests, ~10 min; -m fast = the seconds tier
python scripts/run.py C # term tables -> elements/C/
python scripts/regen_cells.py Ti --method sa-casscf   # one slice of the methods grid
```

Run from the repo root — `pyproject.toml` preconfigures
`pythonpath = src`; no install step. Everything human-readable lands
committed in `elements/<Atom>/`.

## What you get

**Term tables** — every state classified, per Ms sector. Titanium,
CAS(2,5) over {3d}, the Ms = +1 sector (the d² triplet manifold):

```
 #           E (Eh)      <L^2>  L      <S^2>    S  term
 0      -1.46866251  12.000000  3   2.000000    1    3F   (x7)
 7      -1.41396932   2.000000  1   2.000000    1    3P   (x3)

Lz eigenvalues (M_L content): [-3.0 … 3.0]
||[H, Lz]||_F = 1.386e-14    ||[H, L^2]||_F = 5.819e-14
```

The textbook d² multiplet comes out exactly (Hund: ³F ground) with
commutators at round-off — the {3d} window is exactly rotation-closed.
Every atom delivers this per sector; any atom's tables:
`elements/<Atom>/TABLES.md`.

**The methods grid** — every orbital treatment priced against the
atom's ROHF anchor and classified by its diagnostics (Results).

---

## Motivation

The exact Hamiltonian of an atom is SO(3) × SU(2) invariant, so every
term (2S+1)L is exactly (2L+1)(2S+1)-fold degenerate. A CAS run breaks
this in exactly two places.

**The restriction.** CAS takes its core from a prior SCF reference,
confines configurations to a small active set — the **window** — and
represents the omitted core electrons by a folded mean field, f_tu =
h_tu + Σ_core [2(tu|ii) − (ti|iu)]. The partition is exact iff the
window is orthogonal to the core (a window vector leaning into
occupied core orbitals silently harvests Pauli-forbidden
stabilization) — every window is therefore built core-orthogonal.

**Injection point 1 — the fold.** The active–active Coulomb part is a
rotational scalar; **the folded part f carries whatever symmetry the
core density ρ_c has**. An open-shell reference polarizes its core
(the valence quadrupole field mixes the s core orbital with the m = 0
d member), so a raw fold makes Ĥ_CAS merely axial: the term branches
by |M_L| — the σ/π branching of a linear molecule — instead of
staying degenerate. The labels survive (the fold is spin-free);
the degeneracy does not. The cure is a symmetric ρ_c — three core
treatments are compared: **raw** (the SCF's density as-is — the
detector), **averaged** (its l = 0 part — anisotropy removed at
sub-mEh cost), **fractional** (a separate average-of-configuration
SCF — a different functional, not a variational neighbor).

**Injection point 2 — the window.** With P the projector onto the
active space, the projected generators close the algebra iff the space
and its complement are rotation-closed ([P, L̂] = 0):

    [PĤP, PL̂P] = PL̂(1−P)ĤP − PĤ(1−P)L̂P.

This is what the code checks: ‖[P, L̂z]‖ and ‖[P, L̂⊥]‖
(`window_closure`). The condition is **m-completeness** — every (n, l)
shell present with its whole 2l+1 m-set on one shared radial; closure
is imposed by construction, never converged to. A d-shell atom is the
stark case: its CASSCF windows stay d-pure yet break — radial member
mixtures lacking the full m-set; l-purity is necessary, not
sufficient.

**The failure inventory** (both injection points; the SCF reference
sits upstream of both; the instrument can lie on its own):

| Where | Failure | Symptom | Mechanism | Cure |
|---|---|---|---|---|
| the instrument | basis bookkeeping | labels mislead; shells fuse | wrong inputs | label + metric discipline |
| the instrument | sector witnesses | ‖[Ĥ,L̂²]‖ = 0 with defects present | single-irrep sectors; ±pair blind spots | the projector pair |
| upstream | SCF reference | wrong configuration; metastable basins | energy-ordered filling; degenerate open shells | character-forced occupations + energy gate |
| **IP 1** | core field | splittings and commutators 1e-5…1e-1 | the polarization chain | averaged fold / SA / purification |
| **IP 2** | orbital optimization | tilted windows, escape/swap | flat and descending directions; identity change | shell-complete starts + masks + substrate |

Every cure removes the anisotropy at its source: state averaging over
rotation-complete sets (IP 1), promoting the core into the active
space (dissolves IP 1), m-complete windows (IP 2), and verification
itself — ⟨L̂²⟩, the ±M_L pairing, the L̂₊ ladder tests. The detectors
are the commutators and the splittings; a numerically pretty ⟨L̂²⟩ is
not one.

The operator formalism (sectors, the (o, g) symbol table, the ladder
algebra): the Symbols section of `atomterm.operators`' docstring.

## Methods

**ROHF — the reference.** Everything downstream folds it. Where
energy-ordered filling fails (runaway configurations; metastable
polarized basins) it is kept honest by character-forced occupations
(orbital type measured by S-projector weight, never energy order), an
energy gate over cycled starting guesses, and cached trusted
references.

**CASCI — at fixed core and window.** Exact diagonalization of one
sector; the window is a stated start's span; closure is a property of
the span.

**CASSCF — letting the orbitals move.** Minimizes one root over
orbital rotations, run as a grid over starts × masks × substrate.
Masks restrict rotation pairs: window **frozen** (correlation through
the core alone), **radial** (same-(l, member) chains — the
symmetry-preserving channel), **free**. A free run can change window
identity: the **escape** (a rebuilt window matching no term) and the
**swap** (a valence orbital exchanged for a core one); start and
converged windows are compared by span overlap.

**SA-CASSCF — the multiplet average.** One window optimized for the
average over the ground term's M_L components — the rotationally
scalar object that restores the term identity. The kernel returns the
ensemble energy; per-root quantities come from an exact readout pass
in the frozen converged window.

**Partition bookkeeping.** A folded total is an energy *of a stated
partition*; a free mask may exchange core and window orbitals, and the
analysis split then differs from the run's partition — such a number
is an energy of no partition and reads below the kernel's own total
(the fold shift catches it; compare kernel totals only).

**Fold notation.** A fold is named by a **noble-gas element** when it
is exactly that element's closed-shell prefix — He {1s}, Ne
{1s,2s,2p}, Ar {1s..3p}. `semi` marks the **semi-core fold**: for an
open-(n)d atom, the noble-gas prefix with the (n−1)p shell promoted to
the active window and the valence ns pair folded in — on the 3d
metals, Ar − 3p + 4s = {1s,2s,2p,3s,4s}. Every other fold is written
as an explicit shell list, e.g. {1s,2s,2p,4s}. The machine form
(addresses, declarations) is the concatenated list — `1s2s2p4s` —
resolved by `windows.fold_shells`; a lane window carries its fold as
a tag on the label (`3d3p@semi`).

**The treatment grid.** One declared grid, priced and classified at
fixed window: fold (raw/avg/frac) × start (AO/MO/NO) × mask
(free/radial/frozen) × objective (single root / M_L average) ×
substrate (polarized/purified — the core span). Every cell is priced
against the ROHF anchor; the grid localizes each failure by which axis
flips it. The spine (default subset): the three folds; frozen ×
purified under both objectives — the certified class; one free row
per substrate — the basin detectors. The ungated variants (stock MO
picks, free optimization, the raw fold) run through the same analysis
— classified, not discarded.

**Gates.** Character-forced reference + energy gate; checkpoint trust;
the orthonormality guard (C†SC = 1 before any analysis); block-safe
construction (near-linear dependence raises, never NaNs); basin
verdicts read window weights, never noise.

### The pipeline

A new atom is a declared case, executed top-down — per-atom *data*,
not per-atom code:

```
CaseSpec — the declaration (cases.CAMPAIGN)
│  atom, basis/ECP, spin, nelec_act, ncore, shell_patterns, ref, scan
└─ run_case(spec) — the staged pipeline (the term tables)
   ├─ 1. ROHF reference (gated; cached checkpoints trusted as-is)
   ├─ 2. CASCI per sector × fold — term tables vs the combinatorial contract
   ├─ 3. single-root CASSCF (flagged off on multi-solution atoms)
   └─ 4. SA-CASSCF over the ground term's M_L set

The methods grid is the OTHER deliverable path (cells.py): every axis
combination is one addressed CellKey derived from the same
declaration; scans (scans.py) are the optimizer-free anatomy layer.
```

### The variety axes

Every point of the campaign's variety is **data** — an axis value in a
`CellKey` address, never a new script. The declaration says which
slices are interesting; the runners execute the declared set:

| Axis | Values | Run by |
|---|---|---|
| atom × basis × ECP | any `AtomConfig` (all-electron or `ecp=`) | both runners |
| window | any shell set the basis carries; `planned_windows` = declared, not yet run | both runners |
| sector (M_S) | `(na, nb)`, from spin + electron parity | `run.py` |
| method | casci / casscf / sa-casscf | `regen_cells.py` |
| mask | free / frozen / radial / virtmix | `regen_cells.py` |
| start | ao / mo / no | `regen_cells.py` |
| substrate | pol (canonical core) / pur (closed span) | `regen_cells.py` |
| fold | raw / avg / frac | `regen_cells.py` |
| reference | default / fractional-occupation / gated | both runners |

Two runners, two deliverables — `run.py <atoms>` runs the fixed staged
pipeline (the term tables); `regen_cells.py <atoms>` runs the methods
grid and selects within it by key identity:

```bash
python scripts/regen_cells.py Ti                        # full grid
python scripts/regen_cells.py Ti --method sa-casscf     # one methodology
python scripts/regen_cells.py Ti --window 3d4s --mask radial
python scripts/regen_cells.py --group diagnostics       # a declared group
```

Filtered runs **merge** into the existing store (untouched cells keep
their payloads; a merge never prunes); only the unfiltered run
replaces it. Undeclared exploration stays open — `run_atom(...)` for a
one-off staged look, `run_cell(...)` for a single cell: the
declaration records what is *interesting*; it does not gate the
machinery.

---

## Results

All tables live in `elements/<Atom>/TABLES.md` (committed,
regenerable). This section is the reading key and one excerpt.

**Tables grammar** — every element file, same section families in
the same order (headings carry no type labels; the family is the
content):

- **Term table** — the combinatorial term census per CAS/window: no
  energies, route-independent (hand, checked against the solver).
- **Multiplet energies, window lattices, method ladders** — the
  declared campaign as run: `report.md` / `cells.json` — regenerated,
  never hand-copied.
- **CAS ladder** — window families × routes; generated between
  markers, the ladder pins certify.
- **Basin map and studies** — readings and frozen history,
  session-sourced, cited.

Every energy table states its **route** (free — our symmetry-blind
native driver, symmetry measured and restored by averaging; pinned —
symmetry by D2h construction, one `wfnsym` solver per term
component), its **partition** (the fold — Fold notation; a column
when it varies, stated once when uniform), and its **source**, as a
one-line tag under the heading. Per-section legends live nowhere.

Columns of the CI window survey (one row per declared CI cell):

- **window / CAS / dets** — the span and its active size; **fold** —
  the partition: the pinned legs' own tags (`semi`, explicit lists),
  the prefix name for free rows, `·raw`/`·frac`/`·pur` marking the
  free variants' deviation; **start** — the construction (AO shell
  span / canonical MO slice); **route** — free (the native driver) or
  pinned (D2h construction).
- **E** — the free rows' own total; the pinned rows' TERM energy (the
  branch components averaged + the fold-independent bookkeeping
  constant — legs compare across folds). **ΔE** — anchored on the
  primary window's pinned AO row (the tag names it).
- **⟨L²⟩** — the free rows' ground root; the pinned rows' components
  (one value when they agree, a range when they scramble). The
  commutator columns are the assembled meters, em-dash past the
  dense-meter cap.

Columns of a methods table:

- **E** — the run's variational total (Eh; eV in parentheses).
- **ΔE** — gap to the ROHF anchor (the anchor: the caption).
- **Span (Lowdin)** — the converged window's composition per l
  (S-metric fractions of the window's columns — the span, not the
  occupation: the swap reads s → 0).
- **Elec (Lowdin)** — the ground root's electron count per l: what
  the electrons *do* inside the span.
- **⟨L̂z²⟩** — gauge diagnostic, not an observable: an exactly closed
  fold gives a degenerate M_L ground and the solver returns an
  arbitrary member mix; it quantizes to M² only when something lifts
  the degeneracy.
- **Terms** — from ⟨L̂²⟩/⟨S²⟩; the M_L-resolved reading (L̂z
  diagonalized inside the ground eigenspace) is the multiplet
  demonstration.
- **cross-check** — not a table column: the store's
  `max_de_vs_pyscf`, the analysis-ruler agreement with the run's own
  totals. ~0 by construction; an mEh-class residue is core-field
  anisotropy + window definition; a huge value (the swap row's ~800
  mEh) is the partition-bookkeeping shift (Methods) — compare kernel
  totals there, never the analysis split. Values live in
  `cells.json/md`.

Columns of a ladder table (window families × routes):

- **route / method** — the treatment identity: pinned (D2h
  construction, one `wfnsym` solver per term component) or free (the
  symmetry-blind native driver); CI today, CASSCF/SA with the route
  axis.
- **CAS / dets** — the leg's active size and sector determinant
  count; **fold** — the leg's partition (the column that varies).
- **E** — the TERM energy: the mean over the 2L+1 branch components
  plus the fold-independent bookkeeping constant — legs compare
  across folds (not sector grounds).
- **ΔE** — gap to the leg spanning the case's primary window (the
  tag line's anchor), in mEh and eV.
- **⟨L²⟩** — across the branch components: one value when they agree
  (the term reading), a min–max range when they do not (the scramble;
  **ptp** carries its energy side — the component spread, the
  degeneracy certificate).
- **‖[H,Lz]‖, ‖[H,L²]‖** — the assembled-treatment meters, exact up
  to the dense-meter cap (~6k dets); beyond it the columns carry the
  stored STATE RESIDUAL ‖[H,Lz]ψ‖ (the per-root certificate — exact
  zero iff the state is a symmetry eigenstate) where a payload has
  one, else the em-dash — structural zeros for the certified fold
  class.

**Matrix-free meters.** Meters evaluate through *application*, not
matrices: ⟨ψ|Ô|ψ⟩ and ‖[Ĥ,L̂z]ψ‖ need only Ôψ. The algebra fixes the
nonzero pattern — a (one-, two-body) operator connects determinants
differing by at most two spin-orbitals; beyond that the element is
exactly zero (Pauli vanishing, then orthogonality), never computed:

| ⟨D_i|Ô|D_j⟩ | same det | 1 excitation | 2 excitations | ≥ 3 |
|---|---|---|---|---|
| | diagonal rules | single rules | double rules | **0** |

`detbasis.operator_apply` accumulates the same signed couplings the
matrix builder scatters, straight into the output vector (the N²
container never forms; one-body operators ride a singles-only
enumeration); `operators.pyscf_apply` is the same Hamiltonian through
PySCF's own C contractions at any sector size. Below the cap the
matrix forms remain as the cross-certified twins.

The meters, and what each is blind to (ε = ‖ρ_c − ρ̄_c‖/‖ρ̄_c‖, the
core anisotropy — `core_anisotropy()`, no CI involved):

| Diagnostic | Measures | Blind to |
|---|---|---|
| ‖[Ĥ, L̂z]‖, ‖[Ĥ, L̂²]‖ | the assembled treatment | single-irrep sectors (L̂²); ±pair structure (L̂z) |
| ‖[P, L̂z]‖, ‖[P, L̂⊥]‖ | the window span alone | — (as a pair) |
| ⟨L̂²⟩, ⟨Ŝ²⟩ per root | term identity | small defects (order ε²) |
| multiplet splitting | fold anisotropy (order ε) | — |

Class key (row-level discriminators, never assigned by hand):

- **exact** — commutators and span certificates at the 1e-13 class.
- **baseline** — a named non-closure mechanism at its signature scale
  (the raw fold's percent-class commutators at unchanged energy make
  it the polarized-core detector).
- **escape** — ⟨L̂²⟩ ≈ 0 with the Span column showing massive
  borrowing.
- **swap** — SA rows whose Span reads s → 0 with p ≈ 1/6 and
  ‖[Ĥ,L̂z]‖ at O(1): the window's s member exchanged for a core p.
- **gauge** — single-root radial rows: exact shell purity,
  basin-dependent ⟨L̂²⟩, large commutators (the flat member
  direction).
- **construction-class** — MO-start rows whose span certificates never
  close (not an optimizer outcome).
- **NOT-CONVERGED** — the optimizer failed; the row waits for a re-run.

The worked example (Ti, the 3d4s window; excerpt — the full table and
every other atom: `elements/<Atom>/TABLES.md`; frozen, this draw):

| treatment | E, Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ | verdict |
|---|---|---|---|---|---|---|
| ROHF (anchor) | −848.4065 | 0 | — | — | — | — |
| CASCI, averaged | −848.3959 | +10.5 | 12.0000 | 4.2e-09 | 2.3e-08 | exact |
| CASCI, raw | −848.3959 | +10.5 | 12.0000 | 4.9e-02 | 2.6e-01 | baseline (detector) |
| CASSCF, free | −848.4377 | −31.3 | 0.0503 | 4.3e-01 | 1.1e+00 | escape |
| SA-CASSCF, free | −848.4166 | −10.1 | 12.0029 | 1.1e+00 | 6.0e+00 | swap |
| SA-CASSCF, frozen, pur | −848.3948 | +11.7 | 12.0000 | 1.0e-13 | 4.0e-13 | exact |

The **window lattice** is the grid's complement — method fixed
(CASCI, symmetric folds), window grown: the symmetric solution
reaches the broken ones' energies by growing the span (both
transition metals cross below their escaped CASSCF basins at exact
⟨L̂²⟩); the raw-MO spans never close (the certificate disqualifies
them without running any CI); NO spans match the AO value at equal
determinant count and exceed it at the top. Per-series tables: the
element files. Two construction facts, measured: CASCI over the ROHF
MO window equals ROHF exactly (singles Brillouin-dead, doubles dead
unless same-l — the only live channel is the same-l double shell);
and NOs of the idempotent ROHF density are the MOs themselves — the
meaningful NO source is the term-averaged CI density.

## Files and regeneration

| Module | Owns |
|---|---|
| `system` | validation — one atom at the origin, spherical basis |
| `angular` | one-electron L̂ AO matrices — two independent constructions |
| `detbasis` / `secondquant` | many-body matrices over determinants — elementwise twins |
| `operators` | operator assembly, the closure certificates, both PySCF convention boundaries |
| `reference` / `casscf` | the engines — gated SCF; the one optimizer driver with masks |
| `windows` / `core` | window construction + the kind vocabulary; the purified core substrate |
| `natural` | NO mathematics (the NO columns of the lattice) |
| `stages` | the one analysis pass per finished stage — the source of every number |
| `scans` | frozen-landscape measurements — deterministic |
| `analysis` | term labels, configuration summaries, report rendering |
| `cases` / `cells` | the declaration (`CaseSpec` → `run_case`); the cell addressing layer |

Everything human-readable lives in `elements/<Atom>/`: `TABLES.md`,
`report.md/json` (term tables), `cells.json/md` (the treatment grid),
plus frozen one-shot artifacts. New atom = new element file, not
README growth. Regeneration: term tables — `scripts/run.py`
(the per-atom calls recorded in each report header); cells —
`scripts/regen_cells.py` (chkfiles land under `rungs/`, gitignored);
methods tables — `scripts/build_methods_tables.py`; lattice sections —
`scripts/build_lattice_tables.py`; scans — `scripts/radial_scan.py`;
the U lattices — the `run_u*` scripts (the frozen record until the ECP
wave). The test suite is the executable documentation.

## Status and directions

C, Si, Sc, Ti delivered end-to-end; U (ECP) complete as the frozen
`run_u*` record, pending declaration in the campaign (the ECP wave).
A second route to L²-symmetric states is opening: D2h-pinned CI
(one `wfnsym` solver per term component, state averaging over the
complete component set) — symmetry by construction rather than by
measurement. Head-to-head with the classic SA on the p² atoms: the
ensembles are identical to numerical noise, and the pinned route's
optimized windows close orders of magnitude tighter (the
irrep-blocked rotations cannot leak the m-set; the measured record
lives in the tests and ticket 32). On the open d² shell the routes
diverge: the free SA escapes the declared window into the swap
basin (lower ensemble, split multiplet), the pinned route stays in
the declared window with degenerate components — symmetry kept at
a few mEh of ensemble energy.
Open physics: the L = 0 reference mechanism; inner core-pair valence;
the f³d¹/f²d² ordering inside one NO window; a NIST comparison.
Verification program planned: toy axial-p² Λ-branching, the L̂₊ ladder
amplitudes, the stabilizer test, polarization-shell tracking. Forward:
spin–orbit (extends in place — Ĵ², Ĵz as determinant-basis matrices),
transition moments, g-factors, scalar relativity. Further atoms =
element files + data rows.
