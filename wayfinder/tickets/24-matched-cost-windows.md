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

18. **Matched-cost comparison, Sc/Ti (run_tm_methb.py; all both-sides
    exact; chk keys renamed mo_norb<ncas> — self-describing)**:
    equal-size pairs (same determinant count), exact m-complete
    primitive window vs raw ROHF-orbital window:
    Sc: 90 dets: A −759.728891 / B −759.736003 (B by 7.1 mEh);
    605: A −759.739490 / B −759.745047 (B by 5.6);
    1,274: A −759.760397 / B −745353 (A by 15.0).
    Ti: 120: A −848.394585 / B −848.405235 (B by 10.7);
    1,815: A −848.408935 / B −848.414096 (B by 5.2);
    5,096: A −848.438477 / B −848.419572 (A by 18.9).
    VERDICT: opposite small-size behavior to U — on d-electrons B
    wins only the ~100-det point (representability) and A overtakes
    by ~600 dets, reaching the escape at 5k; B's size-gain is again
    tiny (Sc 9.4 / Ti 14.3 mEh over the whole ladder) while A's
    radial freedom gains 31.5/43.9. ⟨L²⟩: A exact (6/12); B ≈ 0-0.2
    (the Ti B ground also sits at ⟨L²⟩ = 0.0000 at w=6 — the same
    L=0-det phenomenon as U's ROHF, now on a d-shell atom).
19. **U matched B-grid complete** (k = 8/15/18/20/22 = 495/3,876/
    7,315/10,626/14,950 dets, matched to 1f1d/2f1d/1f3d/2f2d/3f1d):
    B = −50.3825/−50.3827/−50.3839/−50.3840/−50.3844 (−0.42 → −2.57
    mEh vs the ⁵S° ref det; ⟨L²⟩ 0.13-0.26; ‖[H,Lz]‖ up to 44) vs
    A = −49.4587/−49.7506/−49.7333/−49.8943/−49.9928: B wins every
    pair below 15k (representability), A only at the 4f+2d corner
    (73,815) — where it is exact. The 1f+3d A point: ⁵I° f¹d³
    −49.733255 ×13 (the d-heavy window favors f¹d³, like 1f+1d).
20. **TM crossing tiers**: Sc {3d,4d,5d,4s,4p} (3,249) +0.73 vs
    escape — 5d is flat; Sc {3d,4d,4s,4p,4f} (4,410) −759.769055 =
    −3.06 mEh BELOW the escape: **Sc crosses via the 4f radial**
    (the p-mixing escape is expressed exactly once f-room exists).
    Ti {3d,4d,5d,4s,4p} (18,411) −848.439334 = −4.36 (deeper than
    dsdp's −3.50). BOTH TMs now cross below their broken-symmetry
    basins at exact symmetry — the fact-16b ceiling statement is
    revised to a WINDOW statement with measured crossing points.
21. **Ti-decontracted lattice** (the gated reference, finding 17):
    dalls (2d+2s, 245) −848.400893; dallsp (2d+2s+2p, 3,718)
    −848.432303 = +2.67 vs escape (approaches, not crossed — the
    contracted 4d-inclusive set crosses). Granularity verdict
    consistent with Sc: shell radials ≈ near-optimal primitive
    combos; window-unit granularity is a weak axis.

## Open (next waves)

- **README §7.5**: LANDED (lattice + Δ matrix + representability +
  B-curve + ⁵S° fact + TM plateau table + planned-extension lines);
  the matched Sc/Ti table goes in with the physical naming at the
  next README pass (user: "we will return to README a bit later").
- **Ticket 23** (next): direct single-root CASSCF FIRST, then SA-17,
  then relax/pur-relax.
- **Ticket 22**: NO column of the matched-cost tables (window-units
  registry: shell / primitive / SCF-MO / NO).
- **Ticket 20 entries**: full-matrix slice guard; ReferenceSpec
  exponent-band pattern vocabulary (basis-kind-aware forcing); the
  needle linear-dependence helper; --fresh / key-versioning for chks.
- Deferred: 3rd-d U point, CAS(6,13), NIST.
