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
`elements/<Atom>/TABLE.md`.

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

**Injection point 1 — the core.** The active–active Coulomb part is a
rotational scalar; **the folded part f carries whatever symmetry the
core density ρ_c has**. An open-shell reference polarizes its core
(the valence quadrupole field mixes the s core orbital with the m = 0
d member), so a raw core makes Ĥ_CAS merely axial: the term branches
by |M_L| — the σ/π branching of a linear molecule — instead of
staying degenerate. The labels survive (the core is spin-free);
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
| **IP 1** | core field | splittings and commutators 1e-5…1e-1 | the polarization chain | averaged core / SA / purification |
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

**ROHF — the reference.** Everything downstream builds on it. Where
energy-ordered filling fails (runaway configurations; metastable
polarized basins) it is kept honest by character-forced occupations
(orbital type measured by S-projector weight, never energy order), an
energy gate over cycled starting guesses, and cached trusted
references. **UHF** is the spin-broken sibling of the same discipline
(`run_uhf`: the occupations forced per spin, at the maximal-M_S
sectors the same counts as ROHF — any breaking is closed pairs
polarizing around the open shell); `uhf_contamination` reads the
reference's own spin honesty: ⟨S²⟩ by the overlap formula plus the
corresponding-orbital pairing spectrum (1 = paired; measured: the
[Ar] core stays paired while the valence s pair polarizes most).
⟨Sz⟩ is exact by construction (each spin block conserves its
electron count) — the meter is ⟨S²⟩+σ, never Sz. The basin gate
anchors on `occ_e_ref_uhf`, measured once per atom by the guess
census (`uhf_basin_census`: the production start, the named PySCF
starts, one strong spin break — one call, any atom; the census
measures, the declaration picks). The UHF object is a
reference-class product — the methods grid does not run on it (the
cores assume one closed-shell core density; the broken-core partition is
the declared next step); the anchor row renders from the store's
`uhf` record beside the ROHF anchor.

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
(the core shift catches it; compare kernel totals only).

**The treatment grid.** One declared grid, priced and classified at
fixed window: core field (raw/avg/frac) × start (D/MO/NO) × mask
(free/radial/frozen) × objective (single root / M_L average) ×
substrate (pol/pur). Every cell is priced
against the ROHF anchor; the grid localizes each failure by which axis
flips it. The spine (default subset): the three core fields; frozen ×
purified under both objectives — the certified class; one free row
per substrate — the basin detectors. The ungated variants (stock MO
picks, free optimization, the raw field) run through the same analysis
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
   ├─ 2. CASCI per sector × core — term tables vs the combinatorial contract
   ├─ 3. single-root CASSCF (flagged off on multi-solution atoms)
   └─ 4. SA-CASSCF over the ground term's M_L set

The methods grid is the OTHER deliverable path (cells.py): every axis
combination is one addressed CellKey derived from the same
declaration; scans (scans.py) are the optimizer-free anatomy layer.
```

### The variety axes

Every point of the campaign's variety is **data** — an axis value in a
`CellKey` address, never a new script. Each axis is defined once,
here, in the same shape (values as `name` — plain phrases); the
declaration says which slices are interesting, the runners execute
the declared set, and Results *reads* the axes — it never re-teaches
them:

| Axis | Values | Run by |
|---|---|---|
| atom × basis × ECP | any `AtomConfig` — all-electron or ECP-pseudopotentialized | both runners |
| window | any shell set the basis carries — `planned_windows` = declared, not yet run | both runners |
| sector (M_S) | `(na, nb)` — derived from spin and electron parity | `run.py` |
| method | casci — exact diagonalization of the sector; casscf — one root, orbitals free to move; sa-casscf — the term average over the ground term's M_L set | `regen_cells.py` |
| mask | free — all rotations; radial — same-(l, member) pairs only; frozen — the window pinned; virtmix — the core-mixing probe | `regen_cells.py` |
| start (the window kind) | `ao` → **D** — the distilled window: the reference's radial content expressed in AO form within the declared bands (a band-named shell is its dominant radial, not the whole shell), m-complete (the production kind); `mo` → **MO** — the SCF's canonical slice (a control: m-incomplete by construction, the hybrid-open gauge); `moc` → **MOC** — the SCF's character-selected columns (a control, the same disease); `no` → **NO** — the source window's term-averaged density's occupation groups (store-only: its span *is* the source's, re-gauged — the CI degeneracy and the optimizer equality are its certificate) | `regen_cells.py` |
| source (the radial prior) | untagged — the incumbent: radials distilled from the reference's own occupied orbitals by their span Gram — *which radials the occupied span uses* (on contracted bases ≈ the basis authors' 1:1 contraction); `@r<k>` — the head's *density block*, occupation-weighted — *where the valence density sits* — k radials per shell; `@s<k>` — the same density block averaged over the **configuration series** (one occupation-forced SCF per valence configuration; its first rung, the **head**, is the campaign's anchor reference) — every source is MO-informed (the signal always comes from solved SCF orbitals) and every product is AO-form (explicit basis columns, transferable without the SCF they came from); the two distillation families (span-signal vs density-block) never share a comparison row | the case declaration (`series` + `prior_bands`) |
| radials (k) | k = 1 — the usual case, one radial per shell (the contraction analog); k = 2 — only the series source (the second radial *is* the configuration spread): double budget, the CAS column states it; the tail certificate refuses barren k | the case declaration (`shell_radials`) |
| core (the frozen shells) | a noble-gas prefix — exactly that closed shell (Ar = {1s..3p}); `semi` — the prefix with (n−1)p promoted to the window and the valence ns pair folded in; an explicit shell list — everything else ({1s,2s,2p,4s}) | both runners |
| core field (the frozen density) | `raw` — the SCF's core density as converged (the anisotropy control); `avg` — its spherical average (the default); `frac` — a separately converged fractional-occupation reference (a *different* SCF, not an averaging limit); `lane` — the pinned route's core (pinned-only): the reference's own occupied orbitals projected onto the core span, so free and pinned compare at one construction | `regen_cells.py` |
| substrate (the core span) | `pol` — the SCF's canonical core span (polarized by the open shell); `pur` — the rebuilt rotation-closed span (a raw SCF span never rotation-closes; this one does, at machine precision) | `regen_cells.py` |
| field | neutral (the default) / folded — the reference the cell's substrate lives on; a ninth address slot, written only when folded | `run_cell` / `run_lane_cell` |
| reference | default / gated / UHF / fractional-occupation | both runners |

**Window labels** — the address language, one grammar: shells, then at
most one `@`-suffix, then at most one kind suffix. `3d4s` — the
shells; `@ar` / `@semi` / `@1s2s2p4s` — the core declaration (lane
windows); `@r1` / `@s2` — the source axis: prior family + k; `no2` —
the natural-orbital kind at k = 2. So `3d4s@s2no2` reads: shells 3d4s,
series-average prior at k = 2, natural-orbital kind. Cores count
**electrons chemically**: a core of N
pairs freezes N pairs whatever its AO realization — on multi-band
realizations (multi-zeta, uncontracted) each shell names its declared
band list, and the frozen core is the span's own closed-shell
solution over those bands.

**Reading the source axis — two controlled pairs.** Every source asks
solved SCF orbitals a question and hands the answer back as basis
columns ("the radial directions are data, never a basis rewrite").
The incumbent asks the *span* — the occupied MOs' unweighted Gram;
`@r` asks the *density* of the head configuration — and the head IS
the incumbent's own source SCF, so **D vs @r1** isolates the question
(span vs density) at fixed source. **@r1 vs @s1** asks the same
question of different sources — one configuration vs the series
average — and `@s2`'s second radial *is* the configuration spread the
averaging carries. The `MO`/`MOC` controls take the canonical MOs
themselves as the window — no distillation, m-incomplete by
construction; the `NO` kind (store-only) asks the *correlated*
density — the term-averaged CI's natural orbitals.

**The prior comparison's baseline.** On contracted bases the untagged
incumbent is the authors' 1:1 contraction — the equal-budget baseline
the priors must beat; on uncontracted bases no authors' contraction
exists — the head class is the baseline, and `@s1`/`@s2` measure what
the series averaging adds; the ROHF anchor row is ΔE's absolute zero
in both.

**The prior windows' mean-field quality.** The subspace-ROHF
certificate prices each prior class at the mean-field level: the
valence solved *confined* to the window against the free reference —
same configuration, both sides seeded and gated — so a window that is
expensive here is compromised before any correlation method runs. The
verdict on the uncontracted twins, every valence configuration:
**@s2** is mean-field-sound everywhere (1.3–2.2 mEh — the confinement
costs nothing a correlation method would notice); **@s1** is
filling-dependent (3.3–38 mEh — one series-average radial misses the
d-richer configurations' radials); **@r1** is head-filling-only —
excellent for its own head (+1.4/+1.9 mEh), impoverished for the
d-richer members (+148.6/+199.0 mEh: a d¹-derived radial cannot host a
d³ density). The full 12-row record per twin:
`elements/Sc-unc/logs/subspace.log`, `elements/Ti-unc/logs/subspace.log`.

**Address guarantees** — the address IS the store identity (chk name,
payload key, table row) and never renames: a new axis serializes as
an extra slot only when non-default (an eight-slot address is a
neutral cell and re-serializes unchanged), lane windows carry their
`@core` tag and NO windows their `<source>no<k>` prefix as part of
the window identity, and illegal combinations are rejected at the key
(`frac` on `pur`, `folded` without a core tag, a mask on `casci`).

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

All tables live in `elements/<Atom>/TABLE.md` (committed,
regenerable). This section is the reading key and one excerpt; every
axis it names is defined once, in Methods.

**Tables grammar** — the element file is `TABLE.md`: the ONE methods
table, generated between `methods` markers, its CSV mirror `TABLE.csv`
beside it (same rows, ASCII headers). The table renders the
COMPARISON classes -- the distilled block (the incumbent and the
prior windows), the controls, the treatments;
repeats stay store-only (the prior windows' control rows repeat
the incumbent's at equal width, the frac reference class belongs
to the reference study). Families (the lattice spans, the primary
window's core groups, the lane families, the prior comparison
blocks) are separated by rule rows; the ROHF anchor row leads — ΔE's zero — and the UHF anchor
follows where the store carries it. The tables render what ran;
per-section legends live nowhere.

**The window-composition approaches, reviewed (the U14e record).** The
window is the campaign's central object: which columns enter the CAS
decides both the economics (how much of the reference's relaxation the
window captures) and the symmetry (rotation closure ⇒ definite L). The
approaches, compared on the uranium record:

| approach | columns' source | m-complete | capture | U14e record |
|---|---|---|---|---|
| raw band rows | declared AO bands, identity realization | by declaration | none | the pre-distillation baseline |
| `mo` slice | canonical SCF columns | **no** (hybrid-open gauge) | full | control: ⟨L²⟩ = 0.0019 — rejected |
| `moc` pick | character-selected SCF valence | no | partial | control: the same disease |
| incumbent distillation | occupied-block span Gram, k × m-set | yes (m-replicated) | 10f 0.82 / 13d 0.92 / 7s 0.27–0.71 | (6,13): frozen ⁷F at +1065 mEh — the weak-window record |
| prior axis (`@r`/`@s`) | configuration-density blocks, capture-maximizing | yes | 5f 0.93 / 6d 0.98 / 7s 0.71 | (4,12): ⁵L frozen at +149 mEh |
| lane core route | D2h SCF + the reference's occupied core S-projected onto the declared span | lane-native | full | (4,12) certified 17×72.000000; (4,15) refuses (no window) |
| pin-on-pur (unified) | no SCF — the C1 rung's integrals on a D2h-tagged wrapper + the purified core | yes | as the window's | (4,12): 17×72.000000 at dev 1.1e-13 |
| mixed source | `@r` radials + the reference's own pair MO, s-sector-projected | exact (‖[P,L⊥]‖ 3.3e-13) | pair MO at cos² 0.987 | (6,13): frozen ⁵L at +140.5 mEh — the declared term wins frozen |

The prescription the record states: **per-channel sources by the
channel's character**. A closed pair (l = 0, one m) takes the
reference's own MO, sector-projected — the C1 gauge's m-mixed tail
measured 1.3 % and removed (an l = 0 direction is L-invariant, so
rotation closure survives the projection); open shells (f/d, many m)
take density-block radials rebuilt m-replicated — their SCF orbitals
are irreducibly m-delocalized (a raw-ROHF MO span never
rotation-closes). The single-band fallback for the 7s pair measured
structurally insufficient (the pair MO is s-diffuse across ≥6
primitives; the best band alone carries ~0.27 of it, the declared
three-band span's distilled radial cos² 0.71).

**What buys ⟨L²⟩ = L(L+1) exactly, and what buys 72.0.** Rotation
closure is the exactness certificate: a closed window's CI commutes
with the many-body L², so *every* root reads some L(L+1) at solver
noise — frozen CASCI included (measured 3.6e-06; the pinned block
construction pins to 1e-13). The optimizer adds nothing to exactness.
The *value* 72.0 is the window's economic answer: the ⁵L must beat its
rivals, which takes (i) enough relaxation captured — AO bands froze
the ⁵I at +923 mEh, the prior axis re-ordered the classes (+149), the
pair-MO channel brought the (6,13) ⁵L to +140.5, the frozen winner —
and (ii) no cheaper foreign-configuration class in the sector: the
(4,15) p-promotions are cheaper frozen, so the block grounds leave the
⁵L and the SA's optimization recovers it (71.79–71.96). On
quasi-degenerate manifolds the class partition is window- and
core-sensitive — the (6,13) block reads 12..42 (as-built window,
pur core), 16×90 + one×72 (lane), **17×72.000000** (`@r1`+MOs, pur
core) — the honest window makes the *declared* term the exact ground
of every irrep block.

**Route (free vs pinned).** The *free* route is the native
symmetry-blind driver: one Davidson CI over the whole sector;
symmetry is *measured* afterwards (⟨L²⟩, the commutators, ptp). The
*pinned* route builds symmetry by construction: the declared term's
2L+1 components are distributed over D2h irreps (`branch_d2h`), one
PySCF `wfnsym` solver per component, each component's identity
guaranteed by its irrep — it runs on `lane` core fields, and its
*free* twin (same window, same core, free driver) demonstrates route
agreement at the same partition. The E column's semantics follow the
route: **free rows carry the run's own total; pinned rows carry the
TERM energy** — the branch components' mean at CI plus the
core-independent bookkeeping constant (their ΔE residual vs the
anchor is that constant's projected-density form, mEh-class,
disclosed in the note).

Columns of the methods table (the axis columns — window, CAS, core,
method, route — read their values off Methods' axes; these are the
reading rules):

- **window / CAS / dets** — the span and its active size; the CAS
  shape is the **budget** (a k = 2 prior window doubles the per-shell
  columns — its dominance is variational headroom, never same-budget
  superiority).
- **E** — **rows compare within a family** (same span, same sector):
  free rows the run's own total, pinned rows the TERM energy.
- **ΔE** — against the anchor row that leads the table (the TRUE
  reference — ROHF, or UHF where the campaign anchors there). CASCI
  rows read the partition's own mean-field price: the window cannot
  relax the field it was cut from, so a small or band-identity window
  sits above the anchor by its poverty — the number is the window's
  quality meter, not an error (U14e's band-identity 5f6d reads
  +923 mEh with exact closure and exact ⟨L²⟩).
- **conv** — the optimizer rows' final macro-cycle |ΔE| (order of
  magnitude) per cycle count, `1e-9/7`; CI rows carry the em-dash.
- **⟨L²⟩** — the free rows' ground root; the pinned rows' components
  (one value when they agree, a range when they scramble), with
  **ptp** carrying the energy side of the scramble — the component or
  root spread, the degeneracy certificate (the dash = no set).
- **‖[P,Lz]‖, ‖[P,L⊥]‖** — the one-particle closure certificate of
  the window span (gauge-invariant: rotations within the span cannot
  move it) — zero certifies the span itself is m-closed; cheap at any
  size.
- **‖[H,Lz]‖, ‖[H,L²]‖** — **zero certifies symmetry; the value
  localizes the breaking.** Plain values are the exact operator norms
  (whole matrices — sectors under the dense-meter cap, ~6k dets);
  values marked **ψ** are the STATE RESIDUAL ‖[H,Lz]ψ‖ / ‖[H,L²]ψ‖ —
  the commutator applied to that row's converged solution (exact zero
  iff the state is a symmetry eigenstate) — the only form past the
  cap; the em-dash where neither exists.
- **note** — terse, data-derived flags only (`cycle cap` at the
  optimizer's cycle limit; `m-incomplete` on the MO/MOC rows — an
  explained expectation, never a failure; `bookkeeping form` on the
  pinned rows — their anchor residual's identity; `swap` where the
  optimizer rotated the core block away from the reference core and
  the active orbitals re-pointed at a semi-core-dominated subspace —
  the row's energy is that solution's own (the raw core field, exact
  against PySCF), its ΔE is not comparable to the declared window's
  rows, and the meters are the swap's fingerprint); never fresh
  verdicts.

**Further variants** — the tables render what ran; other combinations
of the treatment axes are one regen away, the first levers when a run
misbehaves (bad convergence, a suspected basin escape). One reading
rule on uncontracted bases, where the axis values shift meaning:
**window names are chemical shells**, but the basis counts radial
primitives — one shell spans several bands (U14e: 5f³ reads on
9f–11f, 6d¹ on 12d–14d, 7s² on 18s–19s, 6p⁶ on 13p–15p; the
certificate lives in the element's `logs/probe_bands.log`), and the
core column names the core's shells when the noble guess would lie
(`core_chemistry`: U14e cores 6s6p7s / 6s6p, not Ne). A band-named
shell in a window is its DOMINANT radial, not the whole shell — the
AO 5f6d7s window's s radial is the compact band-7s function, not the
physical 7s orbital (the MO rows carry the true 7s column). When a
window's table name is chemical, its basis bands demote to the note
column (`bands 10f+13d`) — the address grammar keeps the counters,
the table reads chemically.

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
| multiplet splitting | core-field anisotropy (order ε) | — |

Class key (row-level discriminators, never assigned by hand):

- **exact** — commutators and span certificates at the 1e-13 class.
- **baseline** — a named non-closure mechanism at its signature scale
  (the raw core field's percent-class commutators at unchanged energy
  make it the polarized-core detector).
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
every other atom: `elements/<Atom>/TABLE.md`; frozen, this draw):

| method | route | window | CAS | core | dets | E, Eh | ΔE, mEh | conv | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −848.4065 | 0 | — | — | — | — | forced filling |
| CASCI | free | 3d4s / D | (4,6) | Ar | 120 | −848.3959 | +10.5 | — | 12.0000 | 4.2e-09 | 2.3e-08 |  |
| CASSCF | free | 3d4s / D | (4,6) | Ar | 120 | −848.4377 | −31.3 | 1e-8/38 | 0.0503 | 4.3e-01 | 1.1e+00 |  |
| SA-CASSCF | free | 3d4s / D | (4,6) | Ar | 120 | −848.4166 | −10.1 | 1e-9/33 | 12.0029 | 1.1e+00 | 6.0e+00 |  |

(ΔE against the ROHF anchor; the ellided columns — eV twins, ptp, the
span certificates — are in the element file. The reading: the free
CASSCF escapes the declared window — ⟨L²⟩ collapses to 0.05 — while
the SA restores the term; the classes and their row-level
discriminators are the reading key above.)

The **window lattice** is the grid's complement — method fixed
(CASCI, the averaged core field), window grown: the symmetric solution
reaches the broken ones' energies by growing the span (both
transition metals cross below their escaped CASSCF basins at exact
⟨L̂²⟩); the raw-MO spans never close (the certificate disqualifies
them without running any CI); NO spans match the distilled value at
equal determinant count and exceed it at the top. Per-series tables: the
element files. Two construction facts, measured: CASCI over the ROHF
MO window equals ROHF exactly (singles Brillouin-dead, doubles dead
unless same-l — the only live channel is the same-l double shell);
and NOs of the idempotent ROHF density are the MOs themselves — the
meaningful NO source is the term-averaged CI density.

## Analysis

**PT2 corrections on certified solutions** — the dynamic-correlation
axis on top of any classified solution (CASCI, CASSCF, SA-CASSCF, any
window), via PySCF's own `pyscf.mrpt` (strongly contracted NEVPT2),
consumed as an adapter: one seam takes the substrate
(mol, mf, mo, ncore, ncas, nelecas) and the roots; SA solutions enter
through the multi-root readout pass (the kernel refuses SA solver
objects — the readout IS the required CASCI). Output: per-root,
per-Koopmans-class rows (norm, energy), `E_tot = E_CAS + E(2)`.

**Why E(2) is gauge-visible.** The variational energy is an eigenvalue —
identical for every gauge. E(2) is state-specific: a functional of the
root's own density matrices, and inside an exactly degenerate term
block the root IS a gauge (the solver's arbitrary basis choice within
the pinned subspace). Degenerate perturbation theory diagonalizes the
effective interaction first; state-specific PT evaluates at the vector
handed to it. The correction's invariance is a measured property,
never an assumption.

**Gauge ledger.** Invariant: the block projector, E_CAS, ⟨L²⟩, ⟨S²⟩
(Casimir-type, ∝ identity on a pure-L block). Gauge-visible: ⟨Lz⟩, the
density matrices, E(2). Spin is safe by construction (the Dyall
partition is spin-free); the question lives on the orbital side only.

**The two probes** (they ship with the runs):

| Probe | Measures | Reading |
|---|---|---|
| gauge redraw — a random unitary inside the block | E(2) spread; E, ⟨L²⟩ as controls | the error bar every per-root comparison carries; at the margins of an interpretation → the canonical gauge becomes primary |
| component identity — the pure M_L components | E(2)(M_L) spread | equal by covariance; deviation = a window-closure defect, never gauge |

**Per-L expectations** (exact arithmetic, m-complete basis):

| Block | E(2) across gauges |
|---|---|
| L = 0 | no gauge exists |
| L = 1 | rotations act transitively — constant over all gauges; measured spread = a covariance defect (machinery certificate) |
| L ≥ 2 | rotation orbits do not cover the block — genuine spread is possible (a formulation property, not a bug); the canonical gauge is the primary row, the spread the disclosure |

**Reading the ladder (R and the axes).** The ladder read is a ratio
rule: *R* = E(2) priced on a window against the dynamic correction
the optimized solution of the same span carries — a window whose R
dwarfs its neighbors is *fueling* (the unoptimized span misses
static radial content the correction cannot replace); a ladder whose
PT2-corrected totals are flat across windows has SATURATED — the
CAS-side gain was gauge, not physics.  **The reference-class rule**:
E(2) is orbital-sensitive and rotation-insensitive — broken
references (⟨L²⟩ far from 2S(S+1)) price the same correction as
their kept twins; a refusal class means the store and the readout
disagree about the solution, and the recovery is the regen loop,
never a tolerance change for the occasion.  **The basis-thinness
caveat**: cc-pVDZ-class bases inflate E(2) absolutely (the
per-class compression meter `e_per_norm` reads how hard each class
works); absolute eV values are basis-bound, ladder *differences*
are the usable signal.  **The lane rule**: lane cells are a different payload class (no
ci_stack/core_energy) — their rows price through the lane store's own
readout twin (`cells.readout_lane`: the runner's construction
re-solved and asserted against the record, per branch component).
The lane rows are the D2h-pinned components (the blocks hold the
irrep's own roots; the meters carry their state class) — the core-
partition axis: the Koopmans classes move with the core, and the
row's block spread is the term's across-component E(2) disclosure.
The route-invariance anchor: the lane's ground block vs the free CI
twin's ground root — agreement certifies the irrep pinning is
pricing-inert (band 1e-4 Eh; the observed values live in the gate
log).  **The CAS-size axis** (U14e's promoted trio
(4,12)/(6,13)/(4,15)): the flagship R-test — the window ITSELF is
the axis; a flat corrected total across the trio reads saturation,
a moving one reads the missing static part.  (Ti's same-l CAS-size
windows are VOID-declared — the distilled same-l construction waits
on window_distilled's future seam.)

**The variant axis.** The engine is strongly contracted NEVPT2:
parameter-free, intruder-free by construction (Dyall's H₀ makes
denominators sign-definite orbital-gap sums; the strong contraction
averages each class over its own well-separated Dyall spectrum), and
strictly size-consistent. The alternatives differ on three axes at
once — the H₀ construction, the contraction of the perturber space,
and the presence of empirical dials. CASPT2 (internally contracted,
non-diagonal H₀) is structurally intruder-exposed: its level shift
and IPEA push denominators off zero by hand — mEh-scale dials that
mask exactly the compression information the class meters read; it
joins only as a disclosed variant on an env upgrade (and needs only
the 3-RDM where SC-NEVPT2's contractions scale as the 4-RDM class —
the decisive axis at large CAS). The uncontracted variant (UC:
every Dyall-sector configuration its own perturber) prices the
contraction approximation and is the natural carrier of an explicit
first-order wavefunction — not implemented here (not shipped
in-env); it is the resurrection path should the wavefunction-level
analyzation ever be wanted.

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

Everything human-readable lives in `elements/<Atom>/`: `TABLE.md` +
`TABLE.csv` (the one methods table and its mirror),
`report.md/json` (term tables),
`cells.json/md` (the treatment grid),
plus frozen one-shot artifacts. `elements/ecp/` is the ECP seam's
certification workdir (`tests/integration/test_ecp.py`'s rungs —
gitignored chks, regenerated by the tests). New atom = new element
file, not README growth. Regeneration: term tables — `scripts/run.py`
(the per-atom calls recorded in each report header); cells —
`scripts/regen_cells.py` (chkfiles land under `rungs/`, gitignored;
filtered runs merge into the store, the unfiltered run replaces it);
the one methods table — `scripts/build_methods_tables.py`; scans —
`scripts/radial_scan.py`. The test suite is the executable
documentation.

## Status and directions

C, Si, Sc, Ti, the uncontracted twins, and U14e delivered end-to-end.
The ECP atom is the certification showcase: the free-mask SA legs read
the uranium ground term's parent ⁵L with ⟨L²⟩ = 72.00 exactly and
window closure + commutators at machine zero — the pipeline
**certifies** the atomic ground term, not just labels it (the record:
`elements/U14e/TABLE.md`; single-root optimizations break the same
symmetry by m-incompleteness, the failure the instruments exist to
catch). The D2h-pinned route ran its campaign (C/Si/Sc/Ti on the one
table): the p² ensembles identical to SA to numerical noise, the
pinned windows closing orders of magnitude tighter; on the open d²
shell the routes diverge — the free SA escapes into the swap basin
(lower ensemble, split multiplet), the pinned route stays in-window
with degenerate components. The uncontracted-twin campaign rides the
seeded-reference machinery (the projected 6-31G source + MOM-follow
occupations): fresh SCF on those bases diverges without it, lands the
declared anchors exactly with it.
Open physics: the m-complete 6p promotion (needs the char-selected
fold seam); a NIST comparison round; inner core-pair valence.
Forward: spin–orbit (extends in place — Ĵ², Ĵz as determinant-basis
matrices), transition moments, g-factors, scalar relativity. Further
atoms = element files + data rows.
