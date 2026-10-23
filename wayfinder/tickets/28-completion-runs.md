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

## Acceptance

- Every cell: E (as run), ΔE, ⟨L²⟩ (exact eigh), commutators, span
  pair, cross-check — the T3b column set.
- Cells land in elements/*/TABLES.md between the study markers; the
  T3c product section grows from measured data only (never "—" walls).
- Labels never from PySCF ci vectors (the straggle rule); SA weights
  from the measured ground block (never hardcoded term values).
