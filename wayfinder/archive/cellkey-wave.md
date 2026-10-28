# The CellKey wave — plan (2026-10-27, plan-mode agreed)

Self-contained. Design agreed with the user in conversation (this file is
the record). Execute on the code switch.

PROGRESS 2026-10-27 (b): P0 DONE (5f17df2), P1 DONE (09bdc28), P2 DONE
(2a6132e) — 145 tests green, full suite. The bridge found a real bug:
PySCF CI ordering != detbasis ordering for norb >= 4 (fact 24);
configuration_summary / study diagnostics / run_u_sa were consuming the
wrong order on Ti/Sc/U windows — fixed through the bridge. run_cell
takes mf as REQUIRED (reference construction stays with cases/reference);
virtmix lives on the MASK axis (it is no_core_window_exchange_mask);
start="no" raises NotImplemented pending the parent-prefix policy (P3).
NEXT: P3 (gated on the user's run enable), then P4.
P3 DONE 2026-10-27 (b) (6f22dfe + 63adf28): the no-start policy
(window grammar <source>no<k>, prefix = window identity); the queue
ran -- regen_cells.py regenerated C/Si/Sc/Ti (17 cells each; two run
bugs found+fixed: the budgets->weights leak, the ci serialization ->
ci_stack); build_methods_tables.py renders the Ti manifest FROM
cells.json (unified method-column tables per the user's spec; CASCI
blocks in planned windows; no C/S pairing; no starred gauges).
Suite 148 green. REMAINING: P4 + the campaign continuation through
the seam (lattice cells, optimizer manifest cells, pass curves),
gated per run.
QUEUE UPDATE 2026-10-27 (c) (aa21071..1d51de7): steps (a)+(b) DONE.
(a) the d2 lattice -- 9 CASCI cells (CAS(2,5/8/13) x AOs/MO/NO,
sector (2,0), 4s2 closed in core), the TABLES section rendered from
cells.json by scripts/build_lattice_tables.py. Physics: parity keeps
{3d,4p}=={3d}; {3d,4p,4d} crosses below the anchor; the MO slice is
LOWER in E on O(1)-defective spans (the inversion the lattice
exposes). The seeded anchor row was wrong on both counts (10 dets
not 50; -848.3973 not the campaign number). (b) the FULL Ti manifest
grid -- regen_cells.cell_set is the declared product (radial
ao-only); 33/33 cells + AUDIT verdicts rendered. Headlines: the
escape attractor is start-invariant; free SA from closed starts
swaps, from the MO start escapes; frozen-MOs rows = 'defective
window' (below the closed rows in E -- the inversion inside one mask
group). THREE seam bugs found by the queue, fixed+pinned: the
occupation-grouping m-set split (scale-aware tol, fact 23 amended),
the SA store completion (store_name=), the mask wiring (mo/no-masked
cells had silently run FREE). Suite 150 green.
NEXT cheapest-first: C/Si/Sc grid completion (one regen each; then
their manifests can migrate to the grid tables), the (4,10)/(4,14)
planned windows, pass curves, P4.

## The design

```python
@dataclass(frozen=True)
class CellKey:
    atom: str                    # "Ti" -- the element directory
    window: str                  # "3d4s" -- the NOMINAL span label (a table row)
    sector: tuple[int, int]      # (na, nb) -- fixes Ms and nelec
    method: str = "casci"        # "casci" | "casscf" | "sa"
    start: str = "ao"            # "ao" | "mo" | "no" | "virtmix"  (SINGULAR)
    mask: str | None = None      # "free" | "radial" | "frozen"    (optimizer only)
    substrate: str = "pol"       # "pol" | "pur"
    fold: str = "avg"            # "raw" | "avg" | "frac"
```

- **substrate = the core SPAN** (which core columns exist: canonical vs
  purified); **fold = the core FIELD convention** (raw detector / avg
  ruler / frac = a different reference SCF). Composition matrix:
  (raw,pol) detector; (avg,pol) ruler; (avg,pur) certified; (frac,pol)
  own reference; (raw,pur) DEGENERATE with (avg,pur) -- on pur the
  density is exactly m-scalar (fact 20), so the key CANONICALIZES
  fold="avg" when substrate="pur" (a validator, not a crash).
- Extensibility: every new axis = a field WITH a default; old keys keep
  constructing and keep their address. address(key) <-> key.from_address(s)
  round-trip is the A <-> B mapping (candidate A subsumes B).
- Address grammar (the derived serialization, also the chk name and the
  table-row id):
  `<atom>_<window>_<na>-<nb>_<method>_<start>_<mask|->_<substrate>_<fold>`
  e.g. `ti_3d4s_3-1_sa_no_free_pur_avg`, `c_2p_1-1_casci_ao_-_pol_raw`.
  Order = the manifest's reading order: atom -> window -> sector ->
  method -> start -> mask -> substrate -> fold (where / the physics box /
  the objective / the optimizer freedom / the core). SINGULAR axis
  values throughout (ao, mo, no -- never aos/mos/nos).

## The interface (one deep module, src/atomterm/cells.py)

```python
def run_cell(key, *, mf=None, budgets=None, fresh=False) -> StageFacts
    # resolve each axis through its OWNER (never re-implement):
    #   window+start -> windows (shell spans / mo slice / natural.no_group_windows)
    #   substrate    -> core.purified_core_orbitals | canonical core
    #   mask+method  -> casscf masks + adapters
    #   fold         -> stage_facts(spherical_core=...)
    # cached by address under elements/<Atom>/rungs/<address>.chk;
    # payload embeds payload["cell"] = address (the A<->B bridge)

def load_cell(key, chkdir=...) -> CellData
    # PURE DATA, no live mol needed: mo_coeff (pyscf mcscf group),
    # ci + energies (atomterm groups), and the ordering bridge
    # (pyscf FCI strings -> detbasis sorted tuples).

def analyze_cell(key, mol, mf, chkdir=...) -> StageFacts
    # load_cell + stage_facts: re-analysis of a stored run without
    # re-running the optimizer (the two-function lean, agreed).
```

## Implementation phases (strangler, not big-bang)

- **P0 (independent)**: delete RUNGS.md (stale: ladder naming, dead
  run_ladder API, dissolved-glossary pointer, superseded {M_L=0} vs
  {+/-1} split claim); migrate the per-rung "Question answered" lines
  into treatments.py's module docstring; drop the README pointer.
- **P1**: cells.py -- CellKey + validator (the pur-fold degeneracy rule)
  + address round-trip + CellData + load_cell + the FCI-ordering bridge.
  Tests (fast, pure): address <-> key round-trip incl. the degenerate
  canonicalization; the bridge == detbasis elementwise (regenerate
  pyscf strings for a small norb/nelec, map, compare against a direct
  determinant_basis enumeration).
- **P2**: run_cell -- the resolution table over the owning modules;
  checkpointed by address. Tests (slow, one cheap atom): run_cell on a
  carbon casci cell == the direct run_casci + stage_facts path
  (energies, commutators, payload); load_cell/analyze_cell on the same
  cell == the in-session result (this test IS the bridge certificate);
  one post-hoc property computed from a stored run (e.g. <Lz^2> from
  the stored ci == the payload value).
- **P3**: the queue runs through it -- the spine as list[CellKey] (the
  three folds; frozen x pur both objectives; one free row per
  substrate), regen_rungs.py writes keyed payloads. GATED on the user's
  enable per run (cheapest-first, unchanged order).
- **P4 (after the queue proves it)**: fold the older orchestration in:
  StudySpec cell iteration -> product over key-sets (run_study stays as
  the renderer); treatments._registry composite rung tags -> axis
  values + the resolution table; retire the orphaned campaign scripts
  (diag_s5, diag_mdecomp, run_tm_methb, run_u_* one-offs).

## Vocabulary unification (user direction: 'not at once, but in general')

- **cell** = THE unit of computation/analysis (one grid point: run,
  load, address, table row).
- **study** = a DECLARED COLLECTION of cells + per-cell budgets + the
  render contract (study markers in TABLES). StudySpec redefined over
  key-sets; its name is honest (a study, not a unit).
- **rung** = RETIRES AS A CODE CONCEPT: a rung was a cell with
  window/start pinned (composite tags like sa_relax_pur = method+mask+
  substrate in a string). P4 deletes the Rung dataclass + registry tags
  into axis values. 'Rung' survives as PROSE only (the README
  ladder-of-freedom narrative).
- End state: cell (unit) -> study (declared collection + render) ->
  manifest/table (a view over a key-set). Artifact naming SETTLED
  (user, 2026-10-27): cells.json / cells.md — the rename lands at P3
  together with the keyed payloads (one rename, at the moment the
  content becomes key-addressed); README · Files and the regen script
  move with it. "Rung" then lives only in prose.
- Naming: 'slug' REJECTED (user) -- the derived serialization is the
  cell ADDRESS: key.address / CellKey.from_address(s) /
  elements/<Atom>/rungs/<address>.chk.

## Decisions recorded

- window = NOMINAL label; start = the construction (agreed).
- sector pair (na, nb) -- fixes Ms AND nelec; address form "3-1".
- Singular axis values (user call).
- Two-function lean: load_cell (data) vs analyze_cell (needs mol/mf).
- The code's existing fold variants (name-keyed run_casci stores) stay;
  the key's fold axis names them (user: "ok to keep them").
- Old chks invalidated by the address rename -- acceptable (round-2
  flattening precedent; the queue regenerates).
- RUNGS.md: DELETE (user-approved pending this plan's P0).

## Out of scope this wave

- The run queue itself (gated, executes through P3 once landed).
- stage_io beyond load_cell/analyze_cell (RDM caching etc. -- only when
  a property needs it).
- U(14) (next wave; consumes this seam).
