# Ticket 27: launch & naming polish — conv_accelerator, sa-casscf, study runner

Status: OPEN

Origin: the 2026-10-24 architecture discussion (approved item by item).

## 1. ReferenceSpec.conv_accelerator

Rename `scf_method` → `conv_accelerator`, values
{"diis", "newton", "diis+newton"}:
- "diis" — the plain SCF (PySCF's default DIIS);
- "newton" — the second-order solver from the start (current U
  behavior);
- "diis+newton" — DIIS linearly first, then the Newton switch: run the
  plain-DIIS SCF, continue with `newton()` from that state (verify the
  exact continuation idiom at implementation: mf.run() then
  mf2 = mf.newton(); mf2.run() vs newton-attribute setting).
Update: reference.py, the U CaseSpec data (run_u.py REF), tests.

## 2. sa-casscf naming everywhere

The treatment name is "sa-casscf" (hyphenated) in: CaseSpec.want keys
({"terms", "casci-variants"?, "casscf", "sa-casscf", "treatments",
"scans"}), treatment value sets, T3b row names, table column headers,
README. The engine stays ONE module (casscf.py): SA-CASSCF = the same
constrained_casscf driver + weights (a six-line adapter) — the split
was considered and declined (would duplicate the shared driver).

## 3. The runner: sweep.py dissolves into treatments.py

- treatments.py = the registry (was ladder.py: Rung/_registry) + the
  runner over windows + the render. ladder.py stays a compat shim
  (like pipeline.py/drivers.py).
- Names per approval: StudySpec, plan_study, run_study (NOT
  "sweep"/"grid" — both rejected).
- Treatment value set (the axis table, approved):
  | axis | values |
  |---|---|
  | window | contracted-shell set · primitive band · raw-SCF-MO span · NO group |
  | treatment | CASCI · CASSCF core-relax · CASSCF radial · CASSCF free · SA-CASSCF free · SA-CASSCF core-relax · SA-CASSCF radial |
  | fold | averaged · raw · fractional (applies to EVERY treatment — CASSCF's H carries the core field) |
  | substrate | shell start · purified-core start |
  | budgets | per-cell optimizer knobs (declared) |
  Per cell: E, ⟨L²⟩ (own column), ⟨S²⟩, ‖[H,Lz]‖, ‖[H,L²]‖, span pair,
  ortho_err, PySCF cross-check; light mode = E + span pair.
- "frozen" is NOT a treatment attribute (CASCI is frozen-orbital by
  definition; "frozen-window" is a state any method can have).
- l0..l3 span vocabulary stays as is (metric-exact via lm_sectors; NO
  s/p/d/f aliases — decided).
- The NO-window dependency: derived from the parent CASCI (stage 2);
  the runner sequences it — an intrinsic, documented in code, not in
  the structure diagram.

## 4. README additions

- The launch diagram (v4, approved) into §4:

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
   │      (treatments.py): the windows are walked for frozen CASCI as
   │      well as for the optimizer treatments.
   └─ 6. if "scans" → radial / tilt scans (scans.py)
```

- The T3b mechanism-table format (ticket 26) replaces §7.3's current
  mechanism tables.
- occ_e_ref provenance note: the self-calibration protocol (several
  guesses → classify by character → the passing basin's energy
  recorded as atom data).

## Acceptance

Suite green (123 tests; update the sweep→study test file names and the
AGENTS module table rows accordingly); the diagram + format in README;
no "sweep"/"grid"/"sa"/"scf_method" tokens left in user-facing names.
