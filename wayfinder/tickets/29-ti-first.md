# Ti-first: completion, declaration unification, deep generalization

Agreed in conversation 2026-10-28 (plan mode). This file is the record;
it supersedes the queue/P4 ordering of `cellkey-wave.md` (archived at
`wayfinder/archive/cellkey-wave.md`; its QUEUE UPDATE is the state up
to df6ec4f/dde6071).

## The protocol (user-settled)

1. Work out Ti COMPLETELY first — including deep refactor and
   generalization of the code wherever the Ti work exposes shallowness
   or duplication. Ti is the hard atom on purpose (d², chaos, basins,
   escape/swap); everything the machinery needs to survive, it must
   survive here.
2. Then migrate C/Si/Sc (their grids through the final runner; their
   tables to the grid style). Watch for problems.
3. If problems surfaced: clarify/reformulate the development protocol
   in AGENTS.md. The treatments/study retirement (the last P4 piece)
   is its OWN ticket after the migration — its deliverable is the
   physics test-pin migration (test_treatments' rung fingerprints
   re-pointed at cells), not a silent delete.

## T1 — code only, no gated runs

### T1a: the (mo, pur) builder — the gated branch

- NO SIBLING (user-settled): do not duplicate `build_full_mo`'s
  projection treatment. GENERALIZE it to take window COLUMNS — the
  existing AO-index callers become thin adapters passing the AO
  columns; project-renormalize-assemble stays ONE treatment; the
  mo-pur path is its second consumer. This executes T1c candidate #1
  immediately instead of deferring it. Strangler check: the existing
  paths' behavior must stay byte-equal (the carbon pins already cover
  them).
- The mo-pur consumer: canonical MO columns (window + virtuals),
  re-orthogonalized against the purified core span (S-metric,
  renormalized), assembled `mo_full = [pur core | projected columns]`.
- Docstrings carry the PRINCIPLES: what the pur substrate is (the
  exactly rotation-closed core span, E_radial ⊗ m-set; the column-Gram
  construction and its S⁻¹-wardping trap — fact 20), why MO columns
  need re-projection (orthonormal against the canonical core, not the
  pur one — the fact-17-class input defect), and the honest
  expectation (the window's own m-delocalization survives, fact 22d;
  the likely verdict stays "MO start was not closed" — but as a
  measured statement on an honest substrate).
- Validator: mo×pur flips from void to VALID (pur-fold canonicalization
  applies). The four cells declared in `grid_keys` → render `pending`
  until the gated batch. "We can always reduce their number."
- Fast pins: assembled matrix S-orthonormal (C†SC = I); core block ==
  c_pur; window S-weights shift at the input-defect scale; the
  test_cells validator pin flips with the rule.

### T1b: the declaration unification — one declaration, many renderers

- `cases.CaseSpec` stays the INTENT owner and absorbs the gated
  reference specs (regen_cells.REFS dies).
- `cells` exports the derived key-set FROM CaseSpec, as named groups —
  plain key lists, NO role tags:
  `spine_keys` (established term; keep) / `grid_keys` / `diagnostics`
  (virtmix) / `planned` ((4,10)/(4,14) CASCI blocks).
- THE BRANCH RULE (the switch lives once): renderers may branch on key
  IDENTITY (fields: method, sector, start…), never on presentation
  facts (group, order, emphasis) — those live once in the declaration
  builder's grouping; renderers receive ready-made lists and format.
- AUDIT (the verdict column) stays renderer-side human data — it is
  interpretation (which basin/defect class), not derivable from
  measurements; the renderer asserts every declared cell has one else
  renders pending. (Moving it to a data file next to cells.json is
  optional cosmetics, not structure.)
  [SUPERSEDED 2026-10-30, user-directed: the verdict column is
  ELIMINATED -- the renderer reads payloads only; the classes and
  their row-level discriminators are the reading key in README ·
  Results; the reader classifies. Kills the split source, the
  new-atom audit cost, and the per-draw staleness liability (a regen
  with a flipped basin would leave committed verdicts wrong).]
- regen_cells + both renderers consume the groups; validation =
  byte-diff of the rendered tables vs. current: ONLY the four new
  pending rows may differ. AGENTS module rows updated (cells: the
  groups; cases: the references).

### T1c: the deepening mandate (user note: "deeply refactor and
generalize while working with Ti")

Generalize where Ti's completion exposes it — with guardrails:

- Concrete candidate #1: the projection/orthogonalization family —
  EXECUTED IN T1a (the generalization IS the mo-pur builder; the
  second consumer arrives with it, making `windows.py`'s internal
  substrate seam real — the 2026-10-28 survey's one-consumer
  threshold crossed).
- Candidate #2: whatever the (4,10)/(4,14) windows expose (window
  declarations, bigger sectors, renderer row machinery shared between
  build_methods_tables and build_lattice_tables).
- GUARDRAILS: a generalization needs its real second consumer (one
  adapter = hypothetical seam — no speculative generality); strangler
  discipline (byte-diff renders before deleting the old path); full
  suite before physics-bearing commits; targeted/-m fast otherwise.

## T2 — the gated Ti batch (cheapest-first, on user enable)

### PROGRESS (T1, 2026-10-28) — DONE

- T1a landed, NO sibling: `build_full_mo` takes window columns (AO
  rows stay the shorthand — the one adapter inside; AO-rows ≡
  AO-columns pinned byte-equal on a synthetic substrate). The
  SECOND consumer arrived with it (run_cell's mo-pur branch), and the
  collapse went one step further than ticketed: `Window.columns`
  routes BOTH representations through the one treatment — the c_win
  sibling assembly (hstack + complement, with its latent
  `mf.mo_coeff[:, :ncore]` U-trap default) is gone; for NO columns
  the projection is a machine-precision no-op (pinned). mo×pur valid
  at the key; the four cells declared, rendering pending.
- T1b landed: `cases.CAMPAIGN` is the one declaration (REFS/CASES/
  cell_set in regen_cells retired); `cells.window_grid` owns the one
  render order, `grid_keys`/`diagnostics`/`planned` the named groups;
  the renderer formats by key identity only (address construction
  deleted; WINDOWS-vs-CAMPAIGN coverage asserted). Byte-diff of the
  Ti render: 12 insertions, 0 deletions — exactly the 4 mo-pur
  pending rows per table.
- Found en route (the unification paying for itself): the retired
  cell_set hard-coded `<window>no2` — wrong for single-shell windows
  (C/Si would have raised "no2 is not among them" at their pending
  grid completion); the derived label (`k = shell count`) fixes it —
  carbon's existing `2pno1` pin confirms.
- Suite green (158); AGENTS rows + tests section updated.

1. The four (mo, pur) cells (optimizer runs; expectations above).
2. The (4,10)/(4,14) CASCI anchors, then their optimizer cells.
3. Pass curves (the basin landscape program).

## M — migration (gated)

C/Si/Sc grid completion through the FINAL runner (post-T1 shape);
their Methods tables migrate to the grid style; payloads gain
`electrons` at their regens. Record every friction point — that list
feeds A.

## A — protocol + retirement

If M surfaced problems: reformulate the development protocol in
AGENTS.md (the Ti-first loop made explicit, plus whatever M taught).
Then the treatments/study retirement ticket (pin migration as the
deliverable).

## Open / parked

- The (4,19)/(4,21) lattice rows: ndet² analysis-cost ceiling —
  separate problem, deferred.
- /tmp/opencode pre-session artifacts: user decision pending.
- Verdict vocabulary lives in the AUDIT dict; "MO start was not
  closed" confirmed for the frozen-MO rows. [The dict died with the
  verdict column, 2026-10-30 — the vocabulary lives in README ·
  Results' reading key.]

## Pointers

- map.md · Post-wave survey 2026-10-28 (the review + design round
  that produced this plan's findings).
- message.txt handoff (h) = state at dde6071 (review + survey landed).
- cellkey-wave.md (archived: `wayfinder/archive/`) = the wave's
  history (P0–P3 done, queue steps a/b
  done, its NEXT superseded by this file).
