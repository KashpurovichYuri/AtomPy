# Ticket 06: Generalize selection & pipeline for arbitrary atoms

Status: open

## Question

How to de-duplicate so new atoms (Si, Sc, Ti, ...) need no code copies and
new tests are pure data?

## Resolution (to implement)

**`drivers.py`**

- New `select_valence_ao(mol, patterns: list[str]) -> list[int]`: AO indices
  concatenated per pattern, prefix-matched labels (never `search_ao_label`,
  AGENTS.md fact 1). Keep `select_shell_ao(mol, pattern)` as a one-pattern
  convenience wrapper (update call sites or keep back-compat).
- Everything else (`build_full_mo`, `purified_s_core`, drivers) already
  generic over atom.

**`pipeline.py`**

- New `run_atom(cfg, shell_patterns, nelec_act, workdir="results",
  sa_weights=None, run_casscf=True, nroots_primary=9)`:
  - ROHF (chk-cached) unchanged.
  - `act_idx = select_valence_ao(mol, shell_patterns)`; active = raw shell AO
    columns (identity) placed at `[ncore:ncore+ncas]` — unchanged semantics.
  - **Ms sectors** (derive from cfg.spin = 2S):
    primary = `(na, nb) = ((nelec_act + spin)//2, (nelec_act - spin)//2)`;
    Ms=0 sector = `(nelec_act//2, nelec_act//2)` only when `nelec_act` even
    and `spin > 0` (skip duplicate when spin == 0; skip entirely when
    nelec_act is odd, e.g. Sc).
    Check: C/Si/Ti nelec_act=2, spin=2 -> (2,0) + (1,1); Sc nelec_act=3,
    spin=1 -> (2,1) only.
  - Ms=0 CASCI: nroots = full sector size (complete multiplet table: 9 dets
    C/Si, 25 dets Ti). Primary CASCI: nroots = nroots_primary (default 9).
  - CASSCF single root: run when `run_casscf` (default True; TMs pass False).
  - SA-CASSCF only when `sa_weights` given.
- `run_carbon(cfg=None, workdir="results")` becomes a thin wrapper:
  `run_atom(cfg or AtomConfig(), ["2p"], 2, workdir, sa_weights=[1/3]*3)` —
  test_pipeline.py must stay green unchanged.

**`tests/test_atoms.py`** (new): table-driven, parameterized by case dict:

```python
CASES = {
    "C":  dict(cfg=AtomConfig("C", "ccpvdz", spin=2), patterns=["2p"],
               nelec_act=2, ground=("3P", 1, 1),
               primary={"3P": 3},
               ms0={"1S": 1, "1D": 5, "3P": 3}),
    ...
}
```

Shared helpers: per-atom session fixture (ROHF once), assert term counts
(Counter over states), commutator norms < 1e-8, N diagonal = nelec_act,
ground term & (L, S). Existing carbon-specific tests in test_operators.py
stay (they validate against PySCF FCI); only duplicative ones migrate.

## Acceptance

- All 32 existing tests green, run_carbon output byte-identical.
- New test_atoms.py runs C through the generic path.
