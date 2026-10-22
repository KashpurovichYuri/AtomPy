# Ticket 24: matched-cost window comparison — AO-exact vs raw-ROHF-broken (U + Sc + Ti)

Status: OPEN (execute first among the new wave; light items only)

Origin: user direction 2026-10-21 (post-4f2d discussion). Frame: pure
AO windows make "below ROHF" expectable — enough span trivially does
it; the INFORMATIVE comparison is exact-symmetry vs broken-symmetry AT
MATCHED COST (same determinant count), plus where the two converge as
the window grows.

## The physics questions

1. **Does exactness cost energy at matched cost?** Window A (AO
   primitives, m-complete, ⟨L²⟩ exact) vs window B (raw ROHF columns,
   m-incomplete, ⟨L²⟩ broken) at the SAME det count: energies and
   symmetry fingerprints side by side. Expectation: A ≥ B in energy
   (B's extra freedom = member/asymmetric radials) — measure by how
   much, per atom class.
2. **Where does the fact-16b plateau melt?** (user's full-space
   intuition: exact and broken converge as the window → full space —
   the +24.4/+25.1 mEh radial-vs-escape gap on Sc MUST close once the
   lattice includes the escape's directions m-completely; the escape
   was a window-identity change, and identity is available exactly in a
   wide-enough m-complete window). Measure E_lattice vs the escape
   energy / radial ceiling on Sc (Ti as the second TM point).

## The runs (ALL short; ≤ 30 min/sector, nothing bigger — user rule)

- **U methodology-B scan**: window = the 4 occupied open MOs + up to 4
  character-selected virtuals (S-projector metric — canonical order is
  junk on this basis), admitted one at a time: CAS(4,4) → CAS(4,12)
  class. Per step: E, ⟨L²⟩ off-grid magnitude, ‖[H,Lz]‖, ‖[H,L²]‖ —
  the broken-symmetry-size-vs-extension curve. Minutes per step.
- **Sc/Ti lattices, methodology A**, both basis kinds side by side:
  contracted (shell granularity: {3d}, {4d}, {4s}, … — one converged
  radial per shell; select_valence_ao route) and decontracted
  (primitive granularity; primitive_chain route). Per-block ΔE tables
  (references EXPLICIT: Δf = E(k f) − E(k−1 f) at fixed d; Δd
  analogously) + ⟨L²⟩ columns; read against the escape/ceiling numbers
  of facts 15d/16b.
- **Sc/Ti methodology B at matched det counts**: the question-1
  verdict.

## Deliverables

- Comparison tables under results/ (per atom × methodology × basis
  kind), per-block ΔE matrices with defined references.
- **README §7.x section** (draft reviewed by the user IN-CONVERSATION
  before landing — the standing docs rule): the U window-lattice table
  + the brief per-block ΔE matrix + the representability comment (a
  CASCI bounds the ROHF only when the window span contains the
  occupied MOs; window nesting ≠ determinant nesting; the crossing
  measures how much radial span the mean field uses). Brief.
- §3 entries in ticket 20 for any new friction.

## Out of scope / postponed (user rule: no computationally intensive)

CAS(6,13) on U (even the short (5,1) sector — 16,731 dets — is
postponed, not cancelled); any sector > 30 min; 4f+2d-class runs; the
3rd-d-primitive U point; NO compression (ticket 22); SA wave (23).

## Findings (2026-10-22, in execution)

15. **U methodology-B scan** (`scripts/run_u_methb.py`, ≤ 495 dets,
    seconds/step): (a) the window-B virtual admission MUST rank by
    S-weight on the VALENCE-BAND rows — all-f+d-chain character admits
    the near-core tight MOs (canonical mo123, eps = +38.8 Eh, f-weight
    1.000 on the tight primitives) whose singles sit tens of Eh up and
    never mix (the curve reads flat); (b) NEW FACT: the gated U ROHF
    reference determinant is a numerically pure ⁵S° state — ⟨L²⟩ = 8e-10,
    parity −1 — i.e. the "raw ROHF" ground is an L=0 term 10 mEh above
    the ⁵L°, not an m-broken ⁵L°; (c) the B correlation gain SATURATES
    at −0.24 mEh by k=6→8 virtuals (mean-field-representable but
    correlation-starved); (d) MATCHED-COST VERDICT on U (495 dets):
    B = −50.3825 (0.66 below ROHF, broken: |[H,Lz]| 2.6, ⁵S° content)
    vs A(1f+1d) = −49.4587 (+923 mEh, exact) — exactness costs ~0.9 Eh
    at 495 dets; A catches B only at the 4f+2d corner (73,815 dets,
    −50.3920, exact). The B ceiling (−50.3825) is broken only by the
    saturated lattice — the ⁵L° ground needs radial freedom, not just
    the mean field.
16. **The TM lattices, contracted (shell granularity)** (`scripts/
    run_tm_lattice.py`; sector (2,1)/(3,1); the {d,s}-line with the 4s
    active; pur substrate; every sector ≤ 12k dets both-sides exact):
    (a) the {3d,4s} point REPRODUCES casci_avg_pur bit-wise (+0.00 mEh
    both atoms) — construction validated; (b) **THE PLATEAU MELTS**:
    Ti {3d,4d,4s,4p} (5,096 dets) = −848.438477 = 3.50 mEh BELOW the
    escaped free-CASSCF basin (−848.434973) at exact symmetry
    (|[H,Lz]| 7.2e-13, ⟨L²⟩ = 12.0000); Sc {3d,4d,4s,4p} (1,274 dets)
    = −759.765064 = +0.93 mEh from the escape and descending ~5 mEh
    per axis step — the user's full-space intuition CONFIRMED: the
    escape's correlation becomes available at exact closure once the
    window includes the p directions m-completely; (c) caveat: the
    uniform-nelec {d}/{dd} tiers force d³/d⁴ configurations (+293/+
    173 mEh on Sc) — lattice points on the d axis at fixed nelec, NOT
    the campaign's folded-4s convention; (d) trap fixed in-script:
    build_full_mo returns the FULL matrix — the analysis space is
    mo[:, ncore:ncore+ncas] (a [:, ncore:] slice silently runs a
    CAS(3,34): measured −46.9 mEh against PySCF's correct fold).
17. **Decontracted granularity** (Sc DONE; Ti UNBLOCKED same day): de-
    contract_basis never changes the one-particle span (the ROHF energy
    is IDENTICAL to the contracted value bit-wise) — only the window-
    unit granularity differs. Primitive-band windows (RELATIVE 0.05
    threshold; the absolute one admits the fact-16e needles and crashes
    build_full_mo's Löwdin) reach −34.50 mEh vs the contracted shell
    windows' −36.17 at equal cost: the contracted shell radials are
    near-optimal primitive combos; granularity is a non-axis on Sc.
    **Ti-unc RESOLVED (user direction: try accelerators/guesses)**: the
    whole-l-span forced spec + occ_e_ref = −848.40648875 + NEWTON lands
    the good basin at E = −848.4064887504744 — the pinned contracted
    gate energy to 10 decimals (decontraction preserves the span, so
    the variational minimum transfers). TWO earlier statements
    corrected: the gate energy DOES transfer across basis kinds; and
    the Sc-unc −737 attempt was a DIIS-with-forced-occ convergence
    failure (the U docstring's warning), NOT whole-l smearing — Sc-unc
    plain ROHF was already valid. Operational: the gate only runs on a
    FRESH chk — a cached runaway chk is trusted and short-circuits it
    (cost us one confusing pass).

## Open (next waves)

- **README §7.5**: drafted (lattice table + per-block Δ matrix with
  explicit references + the representability comment) — awaiting the
  user's in-conversation review before landing.
- **Sc/Ti methodology B at matched det counts**: RUNNING (queued chain:
  run_tm_methb.py; matched to ds/dsd/dsdp = 90/605/1274 on Sc,
  120/1815/5096 on Ti). The U matched-cost verdict is in.
- **Sc crossing tiers** (dsdpd5/dsdpf): in the running chain.
- Then ticket 23 (direct CASSCF first), ticket 22, the review wave.
