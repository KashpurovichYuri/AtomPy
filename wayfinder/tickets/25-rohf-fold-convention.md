# Ticket 25: ROHF fold-convention diagnostic — our h_eff vs PySCF's h1eff

Status: OPEN (task zero; blocks ticket 26's cross-check column meaning)

Origin: the 2026-10-24 README/tables review. The user flagged the
two-energy presentation; verification against
`results/Ti_ccpvdz_s2_ladder.json` exposed a REAL open bookkeeping gap.

## The observation (verified)

`stage_facts` computes `e_mine = eigvalsh(our H) + cas.core_energy`
with `cas.core_energy` = PySCF's OWN core constant. If our folded
one-body field reproduced PySCF's `h1eff`, `max|dE|` would be ~1e-10.

- Closed-shell p^2 atoms (C, Si): cross-check = 1e-13 — the
  bookkeeping contract HOLDS.
- ROHF open-shell d-atoms (Sc, Ti): ~90 mEh **on the ground root**
  (Ti casci: our ruler −848.484093/−848.481312 vs PySCF's own
  −848.395940) — the contract does NOT hold; our h_eff and PySCF's
  h1eff differ as OPERATORS by something with a ~90 mEh ground shift.
- The two analysis rulers (raw density vs spherically averaged) agree
  with each other to 2.8 mEh on the ground (the anisotropy effect, the
  splitting scale of §2.4) — the 90 mEh is NOT the averaging.

Suspect: an ROHF-specific fold-convention difference (open-shell Fock
½-factors, density-matrix vs orbital-pair J−K construction), possibly
mixed with the fact-14 window-definition offset.

## Method

One Ti CASCI rung (the gated reference, CAS(4,6) shell window):
1. build our `h_eff` (operators.build_operator_parts fold, raw and
   averaged, with the run's own core columns);
2. build PySCF's `mcscf.CASCI.get_h1eff()` for the same mo (fresh mc
   object — fact 6: get_h1eff re-processes mo_coeff);
3. subtract element-wise; identify the operator difference's structure
   (is it a pure density-matrix difference? an exchange scaling? a
   constant?);
4. separate the three contributions: (i) h_eff operator term,
   (ii) E_core constant term, (iii) the window-definition offset
   (fact 14: raw-vs-orthogonalized columns).

## Acceptance

- The three contributions numerically separated on Ti (and spot-checked
  on Sc + one closed-shell atom as the negative control).
- The convention part FIXED in core.py/operators.py (our fold matches
  PySCF's ROHF convention) or precisely documented as a convention
  difference with its formula.
- The raw-ruler vs PySCF cross-check reaches the bookkeeping class
  (~0) for identical inputs, or the residual is named.
- README §7.3's cross-check legend rewritten on the outcome; the
  elements T3b cross-check column (ticket 26) uses the resolved
  meaning.
- The ⟨L^2⟩/⟨S^2⟩/term-table path is fold-independent — unaffected
  either way (state this in the resolution).

## Notes

- The campaign's own contract ("bookkeeping must be 0", fact 14) is
  what this ticket enforces on the ROHF atoms.
- The as-run energy column of T3b is PySCF's total and is untouched by
  this ticket; only the ANALYSIS ruler's relation to it is at stake.
