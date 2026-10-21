# Ticket 21: U (u_14e) — reference, engine, and the ⁵L° ground by window convergence

Status: PART 2 RESOLVED (bands + tables + stability verdicts; the
4f+2d completion point launched same day). SA/ladder = ticket 23.

Origin: user request 2026-10-20 ("identify its ground state"; Stuttgart
RLC ECP78 + cut-Dyall with high harmonics cut). All numbers measured;
full narrative in results/u_14e/U_casci_report.md; serial tables in
results/u_14e/u_run_part2.log (10 sectors, 2h17m).

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

## Part-2 findings (2026-10-21)

9. **The selector rework**: `primitive_chain` + `select_primitives`
   (exponent walk, band/count semantics, contracted shells raise).
   Bands are contiguous sets around the MEASURED open-MO Mulliken
   peaks (f: 1.44 w1.12 / 0.674 w0.91 / 2.92 w0.53 / 0.294 w0.41;
   d: 0.172 w0.41 / 0.372 w0.30 / 0.076 w0.22) — the kin proxy retired.
10. **Standard tables**: term × parity × ⟨N_f⟩/⟨N_d⟩ × dominant config
   × E per root + per-block max|dE| (population_operator = exact
   S-metric projector; integer per-det N_f = exact config weights).
   CORRECTION of finding 6: the ⁵L° ground is configuration-PURE f³d¹
   (⟨N_f⟩ = 3.000, w 1.00) — part 1's ⟨N_f⟩ 2.33 "f¹d³ mix" was an
   artifact; L = 8 is angular-momentum-unreachable from f¹d³ (L ≤ 5).
11. **The two-axis race, RESOLVED at the symmetric corner** (Ms=2
   block grounds, Eh): f-side growth favors f³d¹ (3f −49.9928 → 4f
   −50.3320, +339 mEh; f⁴ surfaces at +25 mEh); d-side growth favors
   f²d² (2f1d −49.7506 → 2f2d −49.8943 → 3f2d −50.0809). 3f+2d flips
   to even ⁵L by 19.1 mEh — a TRANSIENT of the asymmetric lattice
   point: 4f+2d (73,815 dets, 2h50m) restores ⁵L° f³d¹ by > 37 mEh
   (f²d² ground above the 48th root −50.3547) AND crosses below the
   ROHF: −50.3920 < −50.38185 (measured-gradient extrapolation
   predicted −50.40). GROUND TERM: ⁵L° (5f³6d¹, L=8, S=2, odd), stable
   at every saturated corner. Last open neighbor: the 3rd d primitive
   (0.076, w 0.216); flip-back needs f²d² to gain > 37 mEh MORE than
   f³d¹ from it (the 2nd d gave +119 mEh relative at 3f — not closed,
   but the burden is now on the d side).
12. **PySCF Davidson straggles inside degenerate multiplets**: up to
   38 mEh within exact 13–17-fold clusters (nroots=48, measured vs the
   exact spectrum). Fix: run_casci(fci_opts={"pspace_size": ndet}) =
   their exact branch; both-sides-exact max|dE| 7e-15..8.6e-13.
13. **Dense LAPACK eig is the wall at >= 18k dets** (memory-bound
   ~40-80 GFLOPs effective regardless of threads; 18,411 dets spent
   >45 min incomplete — killed). ARPACK Lanczos (low_roots: k = nroots
   + 16 for the clusters, fixed v0) does 40,920 dets in ~35 min/sector;
   validated against the exact 2f2d spectrum.
14. **255-thread pools are two pools**: OpenBLAS 128 + PySCF's own OMP
   (PYSCF_NUM_THREADS) — set BOTH explicitly at launch; utime-advancing
   threads are NOT proof of progress (spin-wait burns utime too): track
   RSS motion + log checkpoints, kill on frozen-RSS windows.

## Open (next waves — separated into tickets)

- **4f+2d**: RESOLVED (2026-10-21 22:15): ⁵L° f³d¹ restored by > 37
  mEh; ground −50.3920 CROSSES BELOW the ROHF (−50.38185). The ground
  term is stable at every saturated lattice corner; last open neighbor
  = the 3rd d primitive (0.076, w 0.216) — POSTPONED (no heavy runs,
  user rule 2026-10-21).
- **Ticket 24** (next): the matched-cost comparison (AO-exact vs
  raw-ROHF-broken at equal det counts, U + Sc/Ti, contracted AND
  decontracted, short windows only; the Sc plateau-vs-escape question;
  the README §7.x section with the lattice table + per-block ΔE
  matrix). The raw-ROHF CASCI baseline and CAS(6,13) items moved
  there / postponed.
- **Ticket 22**: window methodologies AOs→NOs (the comparison frame:
  primitive bands / replicated SCF radials incl. virtuals / NO
  occupation selection / hybrid NO compression; core/window split as
  a parallel axis; orthonormality + m-complete replication the only
  hard constraints). The U lattice is methodology A's reference.
- **Ticket 23**: direct single-root CASSCF FIRST (the multi-minima/
  term-basin phenomenology), then the SA wave (17-root SA over the
  ⁵L° M_L set); relax/pur-relax; relaxation-vs-growth per axis.
- Then: the review wave (two-axis code-review + codebase-design over
  4f3cf27...HEAD), the architecture discussion, and ticket 20
  (ergonomics) last, on reviewed foundations.
- Deferred: sparse O·v mode, NIST comparison.
