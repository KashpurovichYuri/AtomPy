# Ticket 28: completion runs — the measured window × treatment cells

Status: OPEN (after 25/26/27; needs the study runner of 27)

Origin: ticket 23's measured gaps + the 2026-10-24 systematization
decision.

## The declared studies (cells of the product table)

U (sector (4,0), M_S = 2; NO substrate; exact-eigh labels):
- SA-CASSCF 24 col (NO-24 = the full window): saturation check of the
  SA curve (SA-12 = −50.423810 exact ⁵L; SA-17 degraded) — budget:
  long (17 CI solves × 10,626 dets per macro); the certificate stop
  rule applies (span open ⇒ report, don't push).
- CASSCF core-relax, NO-24: the relax curve's 24-col point
  (relax-12 = −49.899729, relax-17 = −49.900366, both machine-exact).
- SA-CASSCF on CAS(6,13) (7s² window, (5,1) sector): the ⁵L block
  measured by exact eigh (17 roots); does the 7s²-active window change
  the SA verdict (CASCI said inert, −0.58 mEh)?

Sc/Ti (NO-start ladders — the user's comparison interest):
- NO-group windows (from the term-averaged density of each atom's
  largest-window CASCI ground term) × {CASCI, core-relax, SA-CASSCF}
  on the pur substrate; read against the shell-lattice anchors
  (Sc −759.729…−759.769; Ti −848.395…−848.446).

## Stop rules (approved 2026-10-24)

(a) certificate stop — when a bigger window's optimizer gauge breaks
m-completeness (span certificate opens), bigger is worse: report;
(b) saturation stop — when the quantity of interest (term identity,
⟨L²⟩, energy) is window-independent at mEh precision;
(c) cost stop — per-cell budgets; declare long budgets explicitly.

## Execution record (2026-10-24/25)

The studies run through the GENERAL machinery (ticket 27's runner): the
NO arm `atomterm.study.no_group_windows` + `run_case(want={"treatments"})`
with the `spec.study` declaration. No bespoke per-atom scripts (the
first-cut `run_completion.py` was deleted on review -- its knowledge
lives here as declarations).

Sc/Ti NO study (the declaration; U analogous with the band source):

    CaseSpec(symbol="Sc", basis="ccpvdz", spin=1, nelec_act=3, ncore=9,
             shell_patterns=["3d", "4s"], want={"treatments"},
             study={"label": "sc_no",
                    "source": {"patterns": ["3d", "4d", "4s", "4p"]},
                    "windows": ["no-groups"],
                    "treatments": ("casci", "relax", "sa-casscf"),
                    "budgets": {"sa-casscf": {"max_cycle_macro": 60}},
                    "diagnostics": "full"})

U studies (run_case declaration + the one-off six13 probe on the raw
StudySpec -- it needs the custom core, pur minus the 7s column):

    CaseSpec(..., tag="U_stuttgart78_dyall", nelec_act=4, ncore=5,
             ref=REF, want={"treatments"},
             study={"label": "u_no",
                    "source": {"bands": {3: (0.5, 2.0), 2: (0.1, 0.5)}},
                    "windows": ["no-groups"], "no_prefixes": [2, 4],
                    "treatments": ("sa-casscf", "relax"),
                    "budgets": {"sa-casscf": {"max_cycle_macro": 60}},
                    "diagnostics": "full"})

    # six13: ncore=4 (7s out of the pur core), sector (5,1),
    # window = NO {d5, f7} + the 7s pur column; inject tag="study-six13"
    # -> its own section of elements/U/TABLES.md (user decision).

## Cost-stop record (2026-10-25): the U SA-24 cell

Landed: Sc/Ti NO studies (all cells); U relax-24 = -49.918895 (the
relax curve deepens -18.5 mEh from 12->24 cols -- growth, not
saturation; no2 reproduced the campaign's relax-12 to 1 uEh; span
certificates 1e-13 everywhere).

Cost-stopped: U SA-24. Three attempts, zero completed macros:
(a) 4h20m blind (verbose=0, default Davidson) -- runner improvements
came out of this (verbose field, the converged flag captured per
cell); (b) the exact-CI branch (pspace_size=ndet, ticket 23's
straggle fix) -- no macro in 40 min; (c) the exact branch + a
loosened micro solver (conv_tol 1e-6, max_cycle_micro 30) -- no
macro in 31 min. Diagnosis: the 17-clustered-root CI solve in
10626 dets needs >= 1 h per macro in pyscf's python-level direct-CI
-- the macro cadence, not the budget, is the wall. A saturation
verdict needs a specialized SA solver (root-tracking Davidson or a
selected-CI path) -- POSTPONED with NIST/semicore. six13 (17 roots
x 16731 dets) is the same cost class: postponed with it.
Budgets that DID work: relax (single-root, minutes); the small NO
windows' SA cells (Sc/Ti, minutes).

## Acceptance

- Every cell: E (as run), ΔE, ⟨L²⟩ (exact eigh), commutators, span
  pair, cross-check — the T3b column set.
- Cells land in elements/*/TABLES.md between the study markers; the
  T3c product section grows from measured data only (never "—" walls).
- Labels never from PySCF ci vectors (the straggle rule); SA weights
  from the measured ground block (never hardcoded term values).
