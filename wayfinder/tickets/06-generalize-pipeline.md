# Ticket 06: Generalize selection & pipeline for arbitrary atoms

Status: resolved

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

## Resolution (closed)

Implemented as specified, with three deviations/notes:

1. `nroots_primary` is **clamped to the sector size** (`min(9, ndet)`):
   carbon's primary sector has 3 determinants, so the default 9 would be
   illegal; clamping reproduces the historic `nroots=3` exactly while Sc/Ti
   get the full 9. The Ms=0 sector always diagonalizes the whole space.
2. Sector chkfile tags derive from M_S = (na−nb)/2 (`_ms_label`): "1", "0"
   for the historic carbon names, "0p5" for a doublet primary sector (Sc) —
   plain na−nb would collide ("casci_ms0" meaning both Ms=0 singlet and
   Ms=1/2 doublet sectors).
3. The `run_casscf` flag parameter shadows the imported driver function of
   the same name; the import is aliased `casscf_stage` inside pipeline.py
   (public parameter name kept as specced, callers unchanged).

`run_atom` also records `n_dev = max|diag(N) − (na+nb)|` per stage in its
result dict (not in the reports), so tests can pin the particle number.

**Equivalence evidence**: with the SAME ROHF chkfile as input, old and new
code produce byte-identical reports except signed zeros (`0.0` vs `-0.0`)
on two noise-level Lz eigenvalues of the orbital-optimized stages. Two
from-scratch runs of the *unchanged* old code differ far more (CASSCF root
splittings ~1e-2 Eh): threaded BLAS makes SCF nondeterministic in the last
digits, and orbital optimization amplifies that chaotically. Recorded as an
AGENTS.md fact — equivalence checks must pin the chkfile input.

Tests: 33 green (32 + `test_atoms.py::test_atom_term_tables[C]`, the
table-driven generic-path test; Si/Sc/Ti land as new CASES rows in
tickets 07–09).
