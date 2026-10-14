# Ticket 09: Si full analysis — diagnose the escaped basin with the ladder

Status: resolved

## Question

Ticket 07 found Si's CASSCF/SA escaping the 3p shell into an ³S-like well.
With the ticket-08 machinery (window l-weights, per-root energies,
configuration labels, virtmix rung) we can now ANALYZE that escape instead
of skipping it: which configuration does each rung relax into, at what
energy cost, and is the intended ³P basin reachable?

## Resolution (to implement)

1. Generalize `ladder.py` beyond carbon: `run_ladder(cfg, shell_patterns,
   nelec_act, workdir, sa_weights, frac params)`; `run_carbon_ladder` stays.
   Si: patterns `["3p"]`, nelec_act 2, sector (2,0), frac-occ (2,3) —
   n_closed = 6 = [Ne]+3s, fractional 3p block.
2. ΔE column vs reference rung (casci_avg) in the comparison table.
3. virtmix rung: block ONLY active↔inactive rotations (active↔virtual and
   inactive↔virtual allowed) — isolates the valence↔virtual correlation
   energy from identity-changing core↔window exchange.
4. `configuration_summary` (analysis.py): approximate configuration label
   per root — S-metric l-weights of the window columns mapped through the
   dominant CI determinants (exact for raw-AO windows, approximate for
   optimized ones). Method-agnostic.
5. Si ladder rerun: all rungs, basin diagnosis recorded. Per-root CASCI
   extraction shows whether ³P appears among the roots of the escaped
   solution; if not, census via natural-orbital start / mcscf.newton —
   TIMEBOXED: if ³P is not found with standard tools, that is the finding.
6. Checks: canonical core [:6] = [Ne]+3s with no 3p charge; chaotic-band
   test assertions (loose); chkfile store-names distinct per rung.

## Acceptance

- Si ladder table + configuration labels per rung; AGENTS.md records the
  diagnosis; Sc (ticket 10) reuses the machinery unchanged.

## Resolution (closed)

**Headline: the ticket-07 "³S multi-solution well" was a construction bug,
now fixed.** Chain of findings:

1. New ortho_err detector (max|C†SC − I| per rung) showed the shell-start
   CASSCF/SA/virtmix rungs returning MO matrices with orthonormality
   error ~1 — the "escaped ³S well at −303 Eh" (14 Eh BELOW the variational
   bound) was a diverged optimizer fed a broken MO matrix, not physics.
2. Root cause in `build_full_mo`: on Si the raw 3p AOs overlap the
   canonical 2p core MOs by ~0.3 (adjacent radial shells, same l), so the
   complement projector `I − P` was inexact, its spurious singular values
   (0.107) passed the 1e-8 filter, and the returned matrix had shape
   (nao, nao+6) with junk columns.
3. Block-safe construction (fact 14) fixed it. Si now behaves exactly like
   carbon: CASCI textbook-exact (³P/¹D/¹S, raw-core split 3.1e-4 in the
   {0}/{±1} pattern, comm 6e-4); free CASSCF converges in the ³P basin
   from BOTH starts (−0.127 Eh correlation, ⟨L²⟩ 1.58–1.96 by gauge);
   SA restores ⟨L²⟩ = 2.000000 (comm 4e-9). Si CASES row re-enabled
   (run_casscf=True, sa_weights restored) and green.
4. `casscf_relax` on Si: converges now (window p^2.00, ⟨L²⟩ = 2.000,
   comm 6.6e-5). Its dE vs CASCI (+9.1e-2) is dominated by the
   window-definition offset: Si's raw 3p AOs vs PySCF's core-orthogonalized
   window differ physically (carbon's offset ~1e-16). Documented in the
   test and fact 14: max|dE| has up to three contributions (fold
   anisotropy, window definition, bookkeeping-must-be-0).
5. Machinery generalized: `run_ladder(cfg, patterns, nelec_act, ...)` (any
   atom), virtmix rung, ΔE column, `configuration_summary` labels
   (p^2.00 read off determinants exactly for shell windows), sa_pyscf
   rung, per-rung chkfile store-names.

Tests: 5 Si ladder tests + carbon ladder re-validated with the fixed
constructor (45 total green). Ticket 07's "CASSCF/SA disabled for Si"
is reverted; the escape narrative is retired with a pointer to fact 13.
