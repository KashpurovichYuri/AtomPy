# Ticket 21: U (u_14e) — reference, engine, and the ⁵L° ground by window convergence

Status: PART 1 RESOLVED (CASCI story: reference + engine + ground-state
identification). SA/ladder = the next session's wave.

Origin: user request 2026-10-20 ("identify its ground state"; Stuttgart
RLC ECP78 + cut-Dyall with high harmonics cut). All numbers measured;
full narrative in results/u_14e/U_casci_report.md.

## Resolved findings

1. **Reference**: Newton + forced whole-l spans ("l0","l1" closed,
   "l2","l3" open) + occ_e_ref gate = −50.38185233 Eh, 5f³6d¹7s², from
   4/4 guesses (the most reproducible reference in the project). Plain
   DIIS with forced occupations diverges from every start (E > +10³);
   the 1:1 open ranking was the poison (whole-l spans let closed
   leftovers tie) — now 2·w_open + w_closed, verified bit-identical on
   the Sc/Ti-ecp pinned basins.
2. **The label trap, quantified**: uncontracted labels number from the
   TIGHT end ("5f" = exp 36.8). The 5f is compact: open-MO weight on
   kin 3.0/6.5; selection by measurement (fact-1 amendment).
3. **Three-strike canonical-order bug**: range(ncore)/default c_core
   sliced canonical columns = junk virtuals (closed block at 22/61/88/
   120/143); 14.5 and 13.8 Eh discrepancies + a pur crash. All library
   defaults now route through occupied_core_idx.
4. **Fact 4 triggered** (first multi-radial window): overlapping f
   primitives are not S-orthonormal; raw-column analysis gives L²
   eigenvalues 172.98, non-quantized. Analysis uses the Löwdin-orthon-
   alized window (l-blocks survive exactly; ⟨N_f⟩/⟨N_d⟩ consistent).
5. **Pur is mandatory on U**: raw substrate [H,L²] = 2.2 (6p adjacent
   to 5f); pur closes at 5.1e-13 with exact multiplets.
6. **The window-convergence series** (the campaign's core finding):
   1f: f¹d³ bottom → 2f: f²d² (even ⁵L ×17, f³d¹ +52 mEh) → 3f:
   **⁵L° ×17 exactly degenerate at −3.745754, ⟨L²⟩ = 72.00, odd,
   ⟨N_f⟩ 2.33** — the physical ⁵L°₆ parent. Configuration ordering is
   window-sensitive; radial flexibility of f decides it. Parity is the
   exact odd/even separator (pair transfer mixes only within a parity
   manifold — the even f²d² is exactly excluded from the odd ground).
7. **Engine**: slater_condon_matrix vectorized (identity-pinned to the
   kept reference; 18,411 dets hours → 21.7 s; real operators real
   matrices). Regressions exact (2f block grounds to six decimals).
8. **Operational**: 255-thread BLAS auto-pools livelock on our sizes
   (the "silent deaths" were spin-wait contention at 36 GB RSS);
   default 128 threads (user rule). Full-spectrum expectation einsums
   are the O(n³) trap — truncate to the roots read. U timing at 128t:
   95 s ops + 370 s eigh at 14,950 dets.

## Open (next session)

4f stability point + 2d variant; per-configuration term tables as the
standard report; run_u.py upgrade; the SA wave (17 roots over ⁵L° M_L;
does SA radial relaxation replace window growth); ticket 20 execution;
sparse O·v mode for 7s-active windows.
