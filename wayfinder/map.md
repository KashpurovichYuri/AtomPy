# AtomTerm — the effort map

The navigator: what is live, what is open, where the worked examples
live. History is git's job — resolved tickets are deleted, survey
dumps are not kept.

## Tickets

The wave spec: [`specs/pt2-perturbation.md`](specs/pt2-perturbation.md)
— the PT2 analyzation's contract, REWRITTEN PySCF-mrpt-primary
(2026-11-19): the window-ladder diagnostic; the symmetry probes; the
delivery chain.

- [45 — the PT2 runs + PERTURBATION_TABLE.md]: COMPLETE 2026-11-20
  and RETIRED to git history (see the resolved record below); the
  declared-out classes' gated recoveries moved to the Disposition
  candidates; the 2026-11-20 review + the V&V instrument recorded in
  Issues.
- [46 — the TABLE.md physics pass](tickets/46-table-physics-pass.md):
  DESIGN-FIRST, user-gated — the render-level significance policy
  (the stores keep everything; the table keeps the physics).

Resolved (git history the archive): 45 the PT2 wave, COMPLETE
2026-11-20 (2026-11-19..20): the lane seam (`cells.readout_lane` --
the store's readout twin, kernel-ready per-block `mc` out;
`perturb.nevpt2_rows_on` -- the pricing split; `casscf.irrep_solvers`
-- the pin's one construction; `cases.lane_substrate` moved in) +
step 5's ladder on ALL SEVEN atoms (C 29, Si 37, Sc 52, Ti 46, the
uncs 26/26, U14e 32 cells; the combined gate logs/validate_pt2_all.log
-- 1259 roots certified) + the flagship CAS-size read (the corrected
ladder PREFERS (4,12): the (4,15) promotion's CAS gain is dynamic
double-counting, E(2) −341.7 → −336.0 mEh, the corrected total
worsens −50.7655 → −50.7624) + the seam's measured bands
(`readout_mc(ci0=)`, e_noise 1e-6, p_tol 1e-4 — the per-atom tail
sweep Si 2.2e-7 / Ti 3.9e-6 / Ti-unc 5.0e-5; the meters matrix-free
— the dense build hit ~1.35 TB on the M-shell pol twins) + the
wave's reading discipline in README · Analysis.  THE WAVE'S ANSWER:
the optimized ladder SATURATES (the optimizer repairs the poverty
class; the thin spans fuel and PT2 does not repair them — the
remainder is static); E(2) is orbital-sensitive,
rotation-insensitive.  CLOSED OUT with the U14e V&V instrument
(`scripts/u14e_vv.py`, the run instrument; logs/u14e_vv.log green
6/6 — the wave's answers re-derived standalone at machine precision:
E(2) delta 0.0000 mEh, ⟨L²⟩ 72 at 8.9e-15 free / 1.9e-07 pinned,
route delta 5.3e-09 Eh) and the 2026-11-20 two-axis review + design
survey (the fixes + declinations in Issues).  45 step 3 — THE SYMMETRY
PROBES (2026-11-19, `scripts/probes/probe_pt2_block.py`, the certified pur
cell `c_2p_2-0_sa-casscf_ao_-_pur_avg`): E(2) = −58.1340 mEh
(−1.5819 eV) IDENTICAL across the ³P Sz=+1 block; probe A (gauge
redraw, 8 unitaries × 3 roots): controls 1.2e-14 (⟨L²⟩ = 2 exact,
vHv = e_cas exact), spread 2.7e-11 mEh = **7.3e-13 eV** — the
engine's wobble floor; probe B: pure M_L states (−1/0/+1) exact
through the complex-Hermitian Lz block (Lz is anti-Hermitian in the
real p gauge — the .real symmetrization zeroes it; the machine
lesson), the REAL ML=0 component reads the same constant, ±1 are
complex-gauge (covered by A's transitivity).  VERDICT: clean —
per-root rows stand as stored; the degenerate-block gauge class is
CLOSED for SC-NEVPT2 at L=1.  THE NUMBERS: E(2) is ~7.5× the
EN-era preview (−7.7 mEh) — the H₀-choice axis is a big effect
(EN's exact-configuration denominators vs Dyall's orbital gaps);
the class breakdown: Sir (0)' −27.4 mEh leads, then Srsi −10.7,
Sijrs −7.2, Sij −4.9, Srs −4.5, Sijr −3.4; Sr/Si read EXACTLY zero
(perturber norms 0 — no N±1 coupling at this 2e/3o CAS on pur); the
mechanism (structural vs selection-rule) rides the ladder read
across cells/atoms.  The stored block sits 26.6 µEh
below ROHF at the same partition (fact 13's near-equality + the SA
orbital gain); the guard's first real bite: store-flat vs
reader-column CI shapes (fixed: amplitudes are the assertion).
44 the PySCF-mrpt adapter
(2026-11-19): `perturb.py` rebuilt at the PySCF boundary —
`readout_mc` (the run_casci twin without the checkpoint) +
`assert_stored_roots` (fail-fast; energies strict, CI vectors
gauge-aware through the degenerate block's joint projector) +
`nevpt2_rows` (the kernel's canonicalization prelude + the 8
Koopmans classes via mrpt's public functions; Sijrs/eris and the
per-class call-shape dispatch mirrored).  Certified: the fixture's
exact ³P block on the pur core; the RAW-core split pinned at unit
scale (4.4e-5 Eh, the anisotropy detector class); the gauge-tolerant
refusal (a random unitary in the block passes, a zapped root
refuses); the plumbing-vs-kernel twin — 1e-10 on the small CAS, and
AT CAMPAIGN SIZE THE ENGINE'S OWN WOBBLE FLOOR ~3e-7 Eh RUN-TO-RUN
(PySCF's kernel three-run record on O2: −0.16978547 / −0.16978542 /
−0.16978512 — thread-parallel reductions in the big contractions;
the pins sit at that floor); the O2/H14 anchors at 1e-6 of the
examples' retained values; two-root H2 root selection; the SA guard
(kernelled).  13 tests; counts 365/271.  Two unit-tier facts the
build caught: PySCF returns CI vectors as (n,1) column arrays (the
store convention handles it), and the kernel loop WRITES the
canonicalized CI into the shared mc.ci — our plumbing detaches the
list, the twin pins run on fresh readouts.  the PT wave REPLAN — DELETION
(2026-11-19): the spec's "external-anchor-absent" premise was FALSE —
pyscf 2.13.1 ships `pyscf/mrpt` natively (`NEVPT` = sc_nevpt,
strongly contracted; `dfnevpt2`; NO caspt2 in 2.13.1, master-only —
the env fact the original spec missed). The home-grown ladder DELETED
in one commit: `perturb.py` (the FOIS/EN rungs 0–1, the Dyall rung,
the density layer), `test_perturb.py`, `validate_perturb.py`, the
committed `logs/validate_perturb.log` — git history the archive. THE
BUG that forced honesty: `dyall_h0`'s virtual-range `eigh` mixed α/β
virtuals through degenerate eigenspaces (coupling norm 0.016 →
5.4e-11 measured on the carbon fixture) — the committed Dyall script
numbers were INVALID; the 44b density certificates stood but their
consumer was gone (the deletion test passes every rung). EN's one
future perspective — the EXPLICIT first-order wavefunction
(c_I = V_I0/(E0−H_II)) — archived with it; revival = the
wavefunction-level analyzation if the probes demand it. Architecture
now: PySCF owns the PT2 algebra (`mrpt.NEVPT` per root on the stored
roots — the bit-exact readout precedent); we own the adapter + the
meters + the probes + the render (README's Analysis chapter, 1adb82a,
is the reading discipline: the gauge ledger, the two probes, the
per-L expectations). Step 0 anchors: O2/6-31G CASCI(6,8) NEVPT
−0.169785475 (vs the example's retained −0.169785157, 3.2e-7), C2
SA-CASSCF(8,8) → 4-root CASCI → NEVPT(root=1) green in 0.3 s.
Counts 372/278 → 352/258.  43 the validation script (2026-11-18;
the artifacts deleted at the 2026-11-19 replan — this record
describes the DELETED chain, git history the archive; its three
findings carry: the ruler + the gauge class (now formalized in
README · Analysis) + the reach class):
`scripts/validate_perturb.py` — the run instrument (the run.py
precedent), the three-layer chain printed as certificate lines on the
declared C/Si SA scope (grid + planned, ao start, free-route: 14+18
runnable cells, 3 roots each): (a) twins at campaign size (production
coupling vs `pyscf_apply`, EN vs the dense Slater–Condon twin), (b)
the e0 external anchor (the store's PySCF CI energy; the E(2) anchor
absent — the pyscf-forge note prints), (c) whole-space exactness +
completeness (production == dense on a synthetic full-support CI) +
reach containment.  GREEN: C 197 / Si 253 PASS, 0 FAIL
(logs/validate_perturb.log committed).  THREE FINDINGS the
certification caught: (1) THE RULER — the analyzation must fold the
RUN'S OWN core field RAW (`spherical_core=False`, `c_core =
mo[:, :ncore]`): the stored energies are the kernel's values on that
field, and the optimized core's gauge p-admixture moves e0 at the
5e-5 Eh scale between the raw and averaged folds (fact 14 applied to
the post-hoc analyzation; with it the anchor sits at 1e-15).  (2) THE
GAUGE CLASS — the Sz=+1 p² CAS block is exactly degenerate: the
stored SA roots are arbitrary block gauges (support 1/3 dets on the
pur substrates — Davidson's draw), disclosed per root by the script's
gauge line; the meters read the stored draw.  (3) THE REACH CLASS —
on symmetry-adapted tensors (pur) selection rules zero couplings
INSIDE the brute reach (43 of 63): completeness is certified as
production == dense (independent constructions on the same tensors),
the brute count is the disclosed upper bound.  THE LADDER READ
(preview of 45): ground-root E(2) is FLAT across the ladder — C
−7.77..−7.79 mEh (compact/@r1/@s1), Si −6.9 → −5.5 (the @s2 decline)
— against the +250 mEh CAS spread: the external space supplies the
same small dynamic share on every window (cc-pVDZ's valence is thin —
the reading key must carry the basis-thinness caveat); 42 the module born
(2026-11-18; DELETED at the 2026-11-19 replan, git history the
archive):
`perturb.py` rungs 0–1 — the FOIS (the exact apply-pattern
enumeration beside `fois_indices`, the brute spin-preserving
completeness certificate) and the EN E(2) in two independent
constructions (`epstein_nesbet` matrix-free production vs
`epstein_nesbet_dense`, the dense Slater-Condon twin) + the
occupation-sum `diagonal_energies` (the antisymmetrized pair sum —
the symmetric direct form is WRONG on general arrays; the random
non-symmetric integrals caught it, the verification rule's
purpose). Meters: floor + w_ref in `PT2Result`. 7 tests: whole-space
→ 0, one-electron closed answer, the zero-denominator fail-fast, the
twins elementwise (complex non-symmetric, 3 shapes), the diagonal
twin, the cardinality meter, the carbon pin (frozen fixture, NO SCF:
e0 −1.35139 valence-space ruler, e2 −7.689 mEh, floor 0.054,
w_ref 0.996 at n_fois 140); 33 the symmetry-anatomy
deliverable — RETIRED UNEXECUTED 2026-11-18 (the user's call); the
idea stays specified: per-atom Readings sections with (1) the
m-symmetric θ-grid (θ | ΔE | ‖[H,Lz]‖), (2) the member-tilt ladder
exposing the ±pair blind spots (descend-with-zero-commutator = the
directions the sector commutator cannot see), (3) one README reading
key; everything regenerable from `scans.py`, verification = renderer
pins. Would have delivered: the explanatory record behind the
disclosure columns. Unneeded now: the verdicts it would visualize
already live in the map's settled facts, the element files, and the
README rules. Revival = re-open a ticket citing this entry (shape
decisions recorded: descent ladder ~10 rows over the full 31-subset
map; per-atom Readings; the Readings grammar home); 29 the Ti-first protocol; 30 the
module walk; 31 the ECP wave; 32 the symmetry lane (C/Si/Sc/Ti one
table, the route axis, the UHF anchors; closed 2026-11-05); 34 the
folded reference (steps 1–5 + the gap closures; the After's T-C done
as 36A/38, T-D as the 2026-11-16 docs round); 35 the
basis-preparation layer (the design landed, the honest prior table,
the census + the Newton-disease certificate); 36 the runs wave
(batches B/C/D/E ran — the stores carry 66–208 prior-treatment legs
per atom; batch A closed as ticket 38); 38 the U14e close-out (T1 +
T3a/b/c + T4, 2026-11-14..16); 37 the basis-repair bench, the @s
wave, the subspace-ROHF (Stages 0–3, 2026-11-12..16: the seeded-MOM
repair certified + integrated, the twins' series + @s legs 36/36 per
twin, `_unc_s_block` retired, the UHF anchor re-landed, the frac
class measured-refused, `subspace_rohf` integrated with the battery
green and the per-contraction ΔE certificate measured — the record in
the ticket's git history + elements/<twin>/logs/subspace.log); 39 the
U14e prior axis (2026-11-17: the series census measured-REFUSED the
@s half — the ECP substrate holds no member basin, every converged
draw is the head f³d¹; the @r1 legs ran and CURED the band-identity
poverty — +146.9/+149.3 vs the incumbent's +923.1, the exact ⁵L
(72.0000) at machine-zero closure on pur; the optimizer follow-on is
the disposition candidate below); 40 the U14e series REVIVAL
(2026-11-17: the refusal re-attributed — the weight hook's
character-only premium, not the landscape; the members landed under
the self-seeded MOM (`run_rohf(mom_seed=)` + `SeriesSpec
seed="self"`/`open_counts`, certified at ~4e-9, f²d² under
diis+newton); @s1 ran — +141.3/+143.7, the series average BEATS the
head's own block by ~5.6 mEh, exact ⁵L on pur; @s2 MEASURED-REFUSED
at construction — no second f radial, the 5f radial
configuration-inert).

## Issues

### Open

- **The u14e V&V instrument rebuilt and reviewed — 2026-11-27** (the
  code-review skill over 396ddd0..7ff6e47 + the codebase-design
  vocabulary; commits 330b807, 92c5aac, 7274ba5, b64e359, 1f71ee9).
  WHAT CHANGED: the table repaired (one route token `free@pur`, one
  block per method, undeclared keys refused); `_vv_common.py` -> the
  library `scripts/u14e/vv` (a pure move, AST-certified) and then
  restructured around ONE rule -- a State (orbitals + CI of one
  solve) is born only in `vv.routes` and certified at birth
  (`assert_consistent`); the legs became `u14e_vv_ladder.py` (the
  table's one writer over `vv.cases.CASES`), `u14e_vv_uhf.py` and
  the driver; 61 pins in `tests/unit/test_u14e_vv_*.py` (fast tier
  345).  The eight per-CAS legs named in the entries below are
  deleted -- their logs are the records, their code is in git
  history.  For the (4,16) family the table carries the RECORDED run
  (`logs/u14e_vv_4_16.log`): the later duplicate rows had no
  committed log.  FINDINGS, measured, FOR THE USER'S CALL:
  F1 PySCF's CASSCF does not freeze the core (frozen=None rotates
     it with the active and virtual space); the V&V free routes and
     src's `run_irrep_sa` rely on the default.  Core leak of the
     pinned SA: 1.0e-2 ((4,12) class) .. 1.7e-2 ((6,13)); the src
     runner's internal readout (the frozen-core fold with the
     returned active block) sits 1.6e-2..4.1e-2 Eh from the
     consistent one.  "@pur" names the START core.  Freezing is a
     physics decision; for `run_irrep_sa` a src edit.  NOT checked
     whether the production cells' default mask freezes the core.
  F2 A polarized core makes the raw and the averaged core folds
     differ (2.9e-5 Eh on the frozen carbon reference's canonical
     core): compare folds only on the same core columns.
  F3 The (6,13)@r1 "s pair forms" gate reads RED (occupation 0.53 vs
     >= 1.5); the re-scoping proposal stays parked.
  F4 (4,21)@r1 as declared (6p at k=3) is refused structurally: its
     four table keys can never fill -- remove, or keep as a refused
     class.
  F5 E(2) conventions differ: a free row is root 0, a pinned row the
     component mean (3.4 mEh apart at (6,13) for that reason alone).
  F6 (FIXED) The free routes exist only on C1: on a symmetry-tagged
     substrate PySCF's symmetry-adapted classes pin a plain CASCI to
     one irrep (2e-2 Eh on the (4,12) window).  The free routes
     refuse it; the wrapper gate compares the ERIs.
  F7 (FIXED) PySCF's symmetry-adapted CASSCF returns orbitals
     orthonormal only to ~4e-5 (the virtual block; draw-dependent);
     every old pinned row was priced on them unchecked.
     `windows.orthonormalize` runs before the readout.
  F8 (VISIBLE in the V&V library; the src side is the user's call)
     The NEVPT2 prelude -- PySCF's canonicalize(cas_natorb=True),
     which src/atomterm/perturb.py and the V&V library both run --
     RE-PARTITIONS empty active orbitals: a natural orbital occupied
     below 1e-06 is canonicalized together with the virtuals, sorted
     by orbital energy, and the lowest virtuals take its active slots;
     E(2) then depends on the number of empty slots, not on the
     window.  Measured on the frozen carbon 3P ground (p_z empty): a
     non-coupling s orbital moves E(2) by 100 mEh (-265.8 -> -165.1);
     with the window priced as solved, 0.57 mEh.  On U14e 0.17 mEh
     per three slots -- it is what the table's frozen "14p"/"7s" rows
     measure ((4,15) and (4,16) frozen CASCI both read -1056.2x;
     (4,15)@r1/(4,16)@r1 -445.17/-445.13 against (4,12)@r1's
     -445.00).  The V&V library counts it (resorted) and can price
     the window as solved (partition_as_solved); its default stays
     PySCF's convention (the production twin).
  F9 (MEASURED 2026-11-28; what follows is the user's call) The
     FREE state average does not keep a window's composition: in C1
     the optimizer rotates empty non-f orbitals into f.  An s and a p
     channel beside the (4,12)@r1 base converge with f-weight
     0.999..1.000, the term's roots splitting by 0.4..1.3 mEh (<L2>
     71.97 / 71.79); the table's "(4,15)@r1" window (5f6d14p@r1)
     likewise -- its three p columns end at f-weight 1.000.  Those SA
     rows are m-incomplete pieces of a second f shell, not p/s
     windows, and are not comparable at their nominal size.  On
     carbon the free route keeps compositions.
  F10 (MEASURED 2026-11-28) The frozen (4,12)@r1 window's lowest 17
     roots are the 5L term of 5f^2 6d^2 (2.0 electrons in f, 2.0 in
     d), 10.83 mEh below a 5I level: the table's frozen casci row of
     that window -- and the channel screen below -- read a different
     configuration of the same (L, S).  The free state average's
     CONVERGED window holds the ground term (5f^3 6d^1: 3.0 / 1.0;
     its next level is 5K at +9.8 mEh).
  F11 (MEASURED 2026-11-28) F1 priced on the base: with the core
     frozen on the purified span the free state average reads
     -50.4237465 (core leak 5e-16, E(2) -342.686 mEh) against the
     table's -50.4238100 (leak 1.0e-02, -341.670 mEh).  The core
     relaxation buys 0.064 mEh of static energy and costs 1.02 mEh of
     dynamic correlation: 0.95 mEh worse in the total.
  THE INSTRUMENT (4daecb6): `routes.radial_sa` -- the state average
     CONSTRAINED to the window's composition (the production driver
     `constrained_casscf` under the production radial mask, on an
     lm-pure sector basis `windows.lm_basis`; the core frozen on a
     rotation-closed span), seeded with the term's components
     (`term_seed`: a cold Davidson misses components of an exactly
     degenerate term -- 14 of 17 here; its P-space guess does not
     cover every M_L block -- and the cold solves are outright
     diagonalizations up to 20000 determinants), and checked at
     convergence by a fresh exact solve (`_lowest_level`: a
     warm-started solve keeps tracking the set it is given while the
     window holds lower states -- a 5f^2 6d^2 term was tracked 30 mEh
     above lower states of its own window).  Route token
     `radial@pur`: declared in the vocabulary, refused by the store
     until a table block carries it.  11 pins on the frozen carbon
     reference.
  THE WINDOW-VARIANT SWEEP, sizes 12-26 -- MEASURED 2026-11-28
  (`u14e_vv_variants.py` + logs/u14e_vv_variants.log: 105 variants =
  every subset of the declared channels beside the base, all
  converged, all checks passed; by `radial_sa`; core frozen on the
  purified span; mEh against ROHF; E_SA = the state average's
  ensemble energy; sizes 20-26 carry no E(2) -- by F12 it does not
  rank).  The minimal window of each size, by E_SA:
       size  minimal window                 static   gain vs base   E(2) as solved
        12   base [5f, 6d]                  -41.894        --        -342.72
        13   base + s                       -41.894     0.000        -342.69
        14   base + 2s                      -41.894     0.000        -342.68
        15   base + a p radial              -42.543    -0.648        -343.85
        16   base + p + s                   -42.543    -0.648        -343.8
        17   base + a second d radial       -43.827    -1.932        -317.50
        18   base + second d + s            -43.827    -1.932        -317.5
        19   base + a second f radial       -47.827    -5.933        -330.36
        20   base + second f + s            -47.829    -5.935        (not priced)
        21   base + second f + 2s           -47.829    -5.935
        22   base + second f + p            -48.555    -6.661
        23   base + second f + p + s        -48.556    -6.662
        24   base + second f + second d     -52.979   -11.085
        25   base + second f + second d + s -52.981   -11.087
        26   base + second f + second d + 2s -52.981  -11.087
  Every variant of a size converges to ONE window (7p = 8p, 7d = 8d,
  5f = 6f in E_SA to 1e-08; an s channel never moves E_SA): the
  construction of the channel does not matter, only its l.  The added
  span holds 0.002..0.005 electrons: the gains are a static radial /
  angular correlation of the valence shell, small.  ADDITIVITY
  (mEh): p -0.648 alone, -0.728 beside the second f radial (additive);
  the second d radial -1.932 alone but -5.152 beside the second f
  radial (2.7x: synergy -3.2); a THIRD f radial -0.432, a third d
  radial -0.174.  Runners-up: size 22, two second d radials -2.107
  (vs second f + p -6.661); size 26, three f radials -6.365 (vs
  second f + second d + 2s -11.087).
  F12 (MEASURED 2026-11-28; probe_u14e_e2_radials + log) E(2) of the
     (4,N) family is 98% the FROZEN 6s^2 6p^6 7s^2 shell's
     correlation: six of the eight Koopmans classes carry a core hole
     (-336.4 of -342.7 mEh); the valence shell's own correlation into
     virtuals (Sr, Srs) is -8.4 mEh.  It is PARTITION-DEPENDENT: with
     an EMPTY active orbital the core's coupling to that orbital moves
     from the virtual-target classes to the active-target ones (the
     total perturber norm conserved, 1.32 -> 1.31; -95 mEh out, +73
     mEh in), and the net (+25 mEh for a second d radial, +12 for a
     second f radial) is what makes those totals look worse although
     the state is 99.8% the same.  The valence part is consistent: the
     base's PT2 attributes -5.0 mEh to valence -> f'-type and -3.0 to
     valence -> d'-type excitations, which become the static gains
     (-5.93, -1.93 mEh).  CONSEQUENCES (the user's call): NEVPT2
     totals cannot rank windows that differ by an empty active orbital
     -- the state-averaged energy can; the table's "static/dynamic
     trade" across windows (incl. "the corrected ladder prefers
     (4,12)" over the f-enriched (4,15)/(4,16)) is this redistribution
     of the frozen core's correlation; an E(2) column that is 98%
     frozen core says little about the valence shell's dynamical
     correlation.
  THE HF-VIRTUAL CHANNEL SCREEN -- MEASURED 2026-11-27
  (`u14e_vv_screen.py` + logs/u14e_vv_screen.log,
  logs/u14e_vv_screen_spectrum.log; the inventory `vv.channels`):
  frozen CASCI + NEVPT2 per channel beside the (4,12)@r1 base, mEh
  dE_CAS / dE(2) / dE_tot: 5f (the compact 5f radial partner)
  -43.09 / +20.70 / -22.40; 6f -4.58 / +2.04 / -2.54; 7d
  -0.29 / -0.44 / -0.72; 8d -0.12 / +0.71 / +0.58; the p and s
  shells (7p, 8p, 8s, 9s) are SPECTATORS at the frozen level
  (dE_CAS = 0; their E(2) shifts are F8's re-partition).  THE
  READINGS: (i) the p/s channels carry nothing at the frozen level
  -- the proposed p/s rungs are void; (ii) a coupling channel's
  frozen gain is orbital RELAXATION of the occupied orbitals toward
  it (the channel's natural occupations stay ~1e-04), which the
  optimizer recovers on its own -- the screen ranks channels, it
  does not size the rungs; (iii) the base's 5f block already holds
  ~70% of the compact 5f virtuals and its 6d block ~55-60% of four
  "7d" virtuals.  The basis has no g functions (l = s, p, d, f).
  CAVEAT (2026-11-28, F10): the screen's frozen readings are those
  of the 5f^2 6d^2 term of the frozen window, not of the ground
  configuration; its ranking of channels is a frozen-level ranking
  only.  The "(4,19)@r1" rung proposal is superseded by the sweep.
  DECLINED (review): the lane-core route's pinned rows (outside the
  keep-set) are print-only; the "unification" leg's C4 lane-window
  disclosure is dropped (record: logs/u14e_vv_unification.log).
  NEXT (the user's decisions, 2026-11-28): (a) table rows for the
  minimal windows need a `radial@pur` block (vv/table.py BLOCKS,
  GROUPS) + cases in vv/cases.py; (b) the converted (4,15)/(4,16)
  free rows; (c) the frozen 7s^2 pair and the 98%-frozen-core E(2)
  (F12) -- a (6,N) family or an unfrozen-core treatment is a physics
  decision.  Pricing the sizes >= 20 waits on (c).
- **The u14e V&V wave continued — the (6,13) fix + the unification
  ladder — COMPLETE (2026-11-23)** (the lx4 plan: the optimizer
  budgets, the (6,13) @r1 fix, the pin-on-pur unification, the
  backfills).  THE OPTIMIZER BUDGET: `free_kernel` carries
  `max_cycle_macro=200` (the campaign's effective cap was PySCF's
  default 50 — its budgets carry only `sa_weights`) + explicit
  `conv_tol=1e-7`; the conv-column tracker's two bugs fixed (PySCF's
  `imacro` is 1-based — the old `max+1` overcounted; `final_de`
  captured a stale one-cycle-late delta; verified on a probe: the
  old rule read 3 cycles / 1.5e-06 for a 2-cycle / 2.1e-10 run);
  `cells._cycle_trace` counts DISTINCT `imacro`s now (the store's
  callback-count records predating 2026-11-23 inflate by the micro
  overhead — the campaign cells' "200" = 50 macro + ~150 micro; the
  `cycles == max_cycle` cap-hitting identity starts working);
  `run_irrep_sa` gained `mc_opts` (default None = the historical
  path).  THE (6,13) FIX (leg `u14e_vv_6_13_r1.py`, the @r1+MOs
  window): the as-built window's 7s capture is declaration-bound —
  the pair MO (occupied col 143, the promoted 5th pair, occ 2.0) is
  s-diffuse across >= 6 primitives, the declared 3-band span's
  radial captures cos^2 0.71, the single-band fallback ~0.27
  (structurally insufficient — DECLINED); the fix: the s channel =
  the pair MO projected onto its PURE-s sector (the C1 gauge's
  m-mixed d tail 1.3% removed — the raw MO breaks closure at
  8.4e-02, caught by the certificate; restored to 3.3e-13) —
  frozen ground = the 5L ITSELF at +140.5 mEh (occ(7s) 0.53; the
  occ~2 control re-scoping is the user's pending call — the leg
  reads RED on that one check as committed).  THE UNIFICATION
  LADDER (`u14e_vv_unification.py`, green 5/5;
  `logs/u14e_vv_unification.log`): C0 the wrapper probe (the D2h
  clone's integrals bit-identical, ovlp/hcore dev 0.0; CASCI
  agreement 2.95e-08 reproducing the campaign's frozen @r1 record
  -50.2325; the distilled columns irrep-pure at 2.2e-16); C1 THE
  GATE HELD BOTH HALVES (17 x 72.000000 frozen, max dev 1.1e-13 +
  the SA 17/17 at 72.000000, 7 macro cycles) — the certified lane
  record reproduced with NO D2h SCF; C2 pin-on-pur (4,15): the
  frozen block grounds 0..42 (the p-promotion configurations f3d1p1
  cheaper frozen than the 5L) + the SA near-5L (71.79-71.96, 25
  macro cycles) — the FIRST pinned-(4,15) SA; C3 the class trilogy:
  as-built pinpur 12..42 / lane 16x90+72 / @r1+MOs pinpur
  17 x 72.000000 (dev 3.3e-08) — the declared 5L the exact ground
  of every irrep block on the honest window; C4 the lane window's
  pinned SA REFUSED at the labeling seam (promoted columns
  core-span content, measured 3.3e-06 raw / 2.7e-07 post-cleanup —
  the certified no-window class, witnessed at both seams;
  disclosed).  NEW INSTRUMENT: `orthonormal_per_irrep` in
  `_vv_common` (the pre-assembly treatment for blocks whose
  replaced columns overlap the pure ones — the assembly's global
  Loewdin would mix the pure columns across irrep blocks; per-group
  treatment preserves purity at 2.2e-16).  FINDING: the (6,13)
  class is lane-construction-dependent (the counterweight in the
  lane-retirement decision).  README carries the
  window-composition review table + the <L2>-expectation ladder
  (closure buys exactness; the window's economics buy the value).
  THE CONVERGENCE ROUND (the user's review, 2026-11-23): the table's
  method blocks demarcated by === rule rows; the label
  "(6,13)@r1+MOs" -> "(6,13)@r1"; the as-built (4,15)/(6,13) legs
  re-run under the 200-macro budget — the capped "51*" rows
  resolve: the (4,15) CASSCF converges at 150 macro, the (6,13) at
  192 (right at the edge — the flat manifold genuinely needs the
  budget), and the SA<CASSCF inversion the user flagged is
  RESOLVED (the capped draw was the artifact; the converged
  (6,13) CASSCF reads 5.3 mEh deeper and restores the normal
  ordering — 1.9 mEh below the SA ground).  THE @r1 LADDER
  (`u14e_vv_r1_ladder.py`, green 4/4): the matched prior-axis
  source across the CAS axis — (4,12)@r1 reproduces the campaign's
  certified frozen record EXACTLY (-50.23250753, +149.3 mEh,
  <L2> 72.000); (4,15)@r1's frozen total is IDENTICAL to (4,12)@r1
  to 8 decimals — the promoted 14p column is a frozen-level
  spectator (its +3.09 mEh cost appears only in the optimized
  totals); the converged CASSCF/SAs land the incumbent optima on
  every window (the optimizer erases the radial signal's history
  at (4,12)/(4,15); at (6,13) the @r1 draw descends 7.4 mEh
  deeper than the as-built's — the better s channel).  The V&V
  table reads in === blocks; the conv column publishes true macro
  counts everywhere.
- **The u14e V&V family + the scripts regroup — COMPLETE
  (2026-11-21/22)** (the user's directive: several INDEPENDENT
  scripts, one per CAS, each readable top-to-bottom, PySCF-direct;
  metadata from one shared place only): `scripts/u14e/` —
  `_vv_common.py` (metadata, the trusted-rung loader, anchors,
  meters, the derivations' thin helpers, the corrected
  `symmetrize_columns`, the ledger, the row store + the one `VV.md`
  table), the per-CAS legs running the FIVE-route matrix (free
  CASCI — the frozen-window pricing class, free CASSCF, free SA,
  pinned CI, pinned SA), `u14e_vv_uhf.py` (the reference rows +
  (e)), the thin driver.  The free routes are fully independent of
  the production construction chain (the window written out,
  certified elementwise against the production builder at 2.2e-16);
  the meters are the kept analysis; the pinned route and the 6_13
  distilled radial are marked production calls.  THE TABLE (the
  user's schema): the ROHF row leads (the dE zero), the UHF row
  beside it; per row E_method / dE_method / E(2) / E_tot / dE_tot,
  each delta in the `mEh / eV` slash form, <= 4 decimals, no date
  column; E_tot = E_method + E(2) everywhere; the footers carry the
  E_method semantics (the pinned CI value is the TERM energy) + the
  ⟨L²⟩/⟨S²⟩ reading key.  THE SYMMETRIZER (certified standalone;
  the certification caught two bugs before any run: the
  coordinate-only projection overcounts on non-S-orthonormal D2h
  blocks — pyscf keeps raw AO selectors in single-irrep blocks —
  and the projector's transcribe dropped the S): the per-irrep
  METRIC projector `P_k = D_k(D_k^T S D_k)^{-1}D_k^T S`, impurity
  in [0,1] whatever the blocks' gauge; pure → 0 exact, dirtied →
  the analytic value at machine precision, a 50/50 mix raises.
  THE PINNED (4,15) — ATTEMPTED AND MEASURED (the user's
  re-attempt, under the corrected symmetrizer); the mechanism
  CORRECTED AND VALIDATED (the 2026-11-22 second validation, the
  user's challenge — the block-content ranks measured): the
  declared core span of `@6s6p7s` realizes the chemical 6p as
  THREE bands (13p/14p/15p — `shell_bands`), so the promoted
  "14p" is the OCCUPIED SHELL'S OWN MIDDLE RADIAL, not a diffuse
  function (its columns overlap a single 6p MO at 0.87/0.86/0.94).
  The pinned route's frozen core is the C1 reference's m-mixed
  occupied 6p S-projected onto that declared span — it
  contributes THREE independent directions per p block
  (instrumented: B2u rank 3, B3u rank 3 — the committed
  probe `scripts/probes/probe_u14e_vv_span.py`, log
  `logs/probe_u14e_vv_span.log`, carries the ranks, the MO
  overlaps, and the m-pure-feed vanish reproducibly), SPANNING
  THE 14p SLOTS
  ENTIRELY: the promoted columns are core-span content BY
  CONSTRUCTION for the pinned route, the pinned SA has no window
  (`purify_window` refuses every run — the refusal's form is
  draw-dependent: the vanish at the 1e-09..1e-08
  scale or the closed-form tolerance miss ~5e-05, the verdict
  constant), and the pinned
  (4,15) IS the pinned (4,12) at this partition.  The earlier
  mechanisms both RETRACTED (the projector non-commutation; the
  cross-block-cancellation story — the vanishing does not depend
  on the feed's gauge: the M-PURE feed vanishes identically).
  The pinned casci RUNS (the solvers' window keeps the 14p with a
  5.2e-04 core lean — the slots are inside its core's span, the
  block projection cannot remove them) and its block grounds come
  out NOT the 5L — L classes 0..42 (the same route at (4,12)
  reads 17 × 72.000000): the promoted direction DESTABILIZES the
  term identity in the pinned route, the (6,13)
  block-ground-shift class.  The FREE route's core (the pur
  substrate) does NOT contain the slots (survival 0.28) and
  hosts the promotion — its measured cost (+3.09 mEh corrected,
  closure 71.78) is the experiment's answer.  VERDICT: the
  project-then-validate promotion into `lane_core` DECLINED (src
  unchanged); no (6,21)-class partition is needed for now (the
  user's call).  THE
  INPUT-WINDOW GROUND: the frozen CASCI at (4,12)/(4,15) reads the
  certified ⁵I (42.000000, definite L) — the SA promotes it to the
  ⁵L (the record's "~1 Eh gain IS the ⁵L relaxation"); the (6,13)
  frozen CASCI reads the ⁷F (⟨L²⟩ = ⟨S²⟩ = 12.000000, definite
  L = 3, S = 3 — the d⁵ high-spin class, the 7s column
  unoccupied; the leg's earlier "18.0085, NOT definite-L" was the
  shifted-slice meter artifact, RETRACTED above).  THE SPIN
  CROSS-CHECK: the (6,13) ⟨S²⟩ = 12 class verified through
  pyscf's own spin_square (an independent spin construction):
  (S², 2S+1) = (12.00000000, 7) → S = 3, agreeing with atomterm's
  meters at machine precision — the all-parallel d⁵+f¹ septet
  class, the electron count right.  Green: (a0) the ⁵I
  re-derived, (a) −0.0029 mEh, (b) 3.0e-06, (d) +1.5e-09 Eh, (a15)
  +0.0032 within the ±2 mEh band, (c) fresh-vs-fresh +3.091
  (stored +3.086), 6_13 2/2, uhf green (the raw-start draw
  RECOVERED this draw, dev 2.8e-05 — the draw-dependent class
  disclosed).  THE REGROUP: `scripts/probes/` (`probe_pt2_block`
  stays the gate's live dependency), `scripts/generators/`,
  `complete_u14e_table.sh` retired; the moved instruments' usage
  paths lockstep; `collect_*.tmp.log` renamed.  THE DECLINATIONS
  (the 2026-11-22 two-axis review): the shared derivation helpers
  in `_vv_common` (`free_kernel`/`price` — the independence the
  directive wants is vs src/, not between sibling legs); the
  (6,13) converged-False disclosed in the log, not the table (the
  provenance convention); c4ad3e6's `--pinned-only` message
  overclaim stands as history (superseded).  THE DRAW-DEPENDENCE
  DISCLOSURES carried: the free single-root CASSCF rows sit in the
  draw class by physics (⟨L²⟩ 71.82/71.60/22.30 observed); the
  (4,15) free SA hit the cycle cap in one draw (converged False,
  numbers within bands — the log discloses); the
  UHF-seeded SA from a RAW start is a coin toss (one draw descended
  broken, ⟨L²⟩ dev 7.2e+01, git history at 824b5d8; this wave's
  draw recovered); the (6,13) free-SA meters partition the flat
  near-degenerate manifold draw-dependently; the (6,13)
  "assembly defect" is RETRACTED — MEASURED OUT (the 2026-11-22
  triangulation second pass, two sources): the earlier span
  diagnostic's reads were ARTIFACTS OF ITS OWN SHIFTED SLICE (it
  and the 6_13 V&V leg hardcoded NCORE=5 while the key's core is
  ncore_of(14, 6) = 4 pairs — the leg's meter slices ran one
  column off; the KERNELS were never affected: PySCF derives
  ncore=4 from the sector).  The assembled span IS the projected
  distilled span exactly (principal cosines all 1 against
  (I−P)act; the one 0.8536 tilt = the honest 7s projection
  survival), the closure meters ~1e-13 on act/projected/assembled
  AND on c_pur, and the CAMPAIGN CELL agreed all along
  (`cells.json` (6,13) casci: ground_l2 12.00000003, win_lz/lperp
  ~1e-13).  The TRUE frozen (6,13) ground: the 7F — L = 3, S = 3
  definite (⟨L²⟩ = ⟨S²⟩ = 12.000000; the occupations measured:
  five d columns at ~0.99 + one f at 0.998, the 7s column
  UNOCCUPIED at 3e-03): the as-built window's variational ground
  abandons the intended 5f³6d¹7s² configuration — driven by the
  WEAK DISTILLATION CAPTURES (10f 0.82, 13d 0.92, 7s 0.27–0.71),
  not by any assembly corruption.  THE WINDOW REWORK the user
  intuits lands HERE: a larger s k in the (6,13) distillation
  declaration (k_s = 2 or 3 of the 3 s bands — the user's
  declaration change), a new-case decision.  THE LEG FIXED:
  NCORE = 4 + the key-identity assert; the SGS replacement of
  build_full_mo's Loewdin DECLINED — no measured defect (the
  gauge never corrupted a span; the windows are healthy at
  machine precision).  The (4,15) lean/verdict paragraphs above
  are untouched (their mechanisms are lane-side and stand).
- **The 2026-11-20 wave review + design survey** (ticket 45's close;
  the two axes + the codebase-design vocabulary over `src/`, fixed
  point a0fe3fe): fixed = the renderer's ground-block spread got one
  home (`_ground_spread_meh`); the dead `size` branch retired (no
  manifest declared it); the lane block-support `assert` → a
  fail-fast raise (rule 5); `key_core`/`lane_core`/`window_columns`
  promoted public (the store-free V&V instrument is the third
  consumer; the docstring refs already pointed there; `__all__` +
  the layout row + `kernel_cell_payload`'s row gap, one lockstep
  commit); test_perturb's band labels run (a)-(d); the committed
  PERTURBATION_TABLEs were one blank line short of the renderer's
  own bytes (pre-existing drift, surfaced by the re-render, lockstep
  restored); AGENTS/README test counts lockstep (376/282 — the
  handoff's 374 was stale at 67ab1da).  DECLINED (judgement calls,
  recorded): `lane_substrate`'s docstring carries measured instances
  (the repo's exemplars land measured values in Notes-style prose —
  `run_rohf`'s own Notes cite Ti numbers; the slot rule reads those
  as Notes, not contract); `lane_substrate`'s basin-steer retry
  branch untested (campaign chemistry, run-verified — a pin would be
  an hours-class test costume); `_meters`/`EH_EV` duplicated across
  the standalone instruments (the run.py precedent: scripts are
  self-contained run instruments); `_lane_block`'s 13-parameter
  clump (a block-context bundle type — a later refactor, not this
  close); `validate_pt2._probe`'s string-parse coupling to
  probe_pt2_block.py's stdout (gate + probe are committed together,
  one home); `validate_pt2`'s `FAILS` module global (script style);
  the (6,13) leg print-only in the V&V (the handoff's assertion list
  is (a)-(e); the (6,13) meters carry the certified L=9 disclosure);
  the UHF-seeded start MO is the symmetric window (the handoff's
  literal wiring — the reference feeds PySCF's kernel; a
  broken-orbital-start variant is a follow-up if the spin-broken
  START class is ever wanted).  THE V&V RUN: logs/u14e_vv.log green
  6/6 — E(2) delta 0.0000 mEh, ⟨L²⟩ 72 at 8.9e-15 free / 1.9e-07
  pinned, the (4,15) worsening +3.081 mEh (stored +3.086), route
  delta 5.3e-09 Eh, the UHF sibling printed (E −6.961 mEh, ⟨S²⟩
  6.229, all five occupied pairs broken); the (6,13) meters
  reproduce the certified L=9 disclosure.
- **The 2026-11-19 review's declinations** (the two-axis review over
  the PT2 wave, fixed point 7ab52df): fixed = probe_tmp.py scratch
  deleted; the AGENTS cases.py layout row (lockstep); the measured
  values moved out of assert_stored_roots' contract slot; the
  shape-mismatch refusal message (a numpy broadcast error used to
  replace the declared text); e_tol/c_tol force-pinned; the
  compression meter (ClassRow.e_per_norm) landed; the absolute
  fixture E(2) pin (observed −0.26575851 Eh — which reproduced the
  campaign poverty-class cell's E(2) from an independent substrate
  rebuild, ~1e-7 Eh); H2 build one home; probe-A onto _expect; the
  spec letter amended to the declared stored-roots form + the mrpt
  internal-prelude honesty.  DECLINED (judgement calls): the probe's
  _root_row private reach (one instrument; a promoted seam would
  widen nevpt2_rows' interface for a single consumer — the deletion
  test fails it); cell_substrate returning a bundle type (fights the
  declared substrate-signature convention); the _VARIANTS guard
  (the registry deferral is documented); readout_mc's unused `mol`
  (the substrate-signature convention).
- **The future-wave queue** (ticket 34's After + the standing
  candidates): the contraction-recipe provenance review (our radials
  vs the basis authors' coefficients); then spin–orbit, transition
  moments, g-factors.  (The term-ID round — splittings vs NIST —
  DECLINED 2026-11-18, the user's call.)
- **The C/Si pinned-lane legs — APPROVED 2026-11-18** (the 2+2
  pending rows, the last non-twin pendings): `--group lane` on C and
  Si; minutes-class; closes the lane families completely.
- **The queued CASCI review — EXECUTED 2026-11-17** (the three
  classes classified; no analysis bugs found): (1) the window-poverty
  class (+923 mEh) is the correct declared reading — U14e's
  band-identity 5f6d holds the compact radial, poverty is radial not
  symmetry (pur payload: exact closure + exact ⟨L²⟩ at +923.1), and
  the class is the normal CASCI shape campaign-wide (C/Si/Sc rows at
  +180…+615 with exact certificates); fixed by the README ΔE rule
  (the partition's mean-field price) + the D-window band clause.
  (2) the reference-echo roots (ΔE ≡ 0: the reference determinant's
  singles Brillouin-dead, the CI reads the mean field back) are
  presented honestly — Ti-unc's (2,5) rows carry them, Sc-unc has
  none (d¹ lacks the same-l double channel); noted in the twin's
  element file.  (3) the raw-fold ⟨L²⟩ mixing is the declared
  (avg, pol)-ruler anisotropy — U14e's 41.9649 is the honest ground
  root (⁵I on all three fields, no root swap; pur exact 42.0000);
  noted in U14e's element file.  Candidate left open below (the
  renderer's note-class).
- **The folded-field campaign identity** (ticket 34's open
  decision): the substrate tag in the key grammar, the free-route
  dominance declaration, the table re-anchoring — user-visible
  before any table carries folded rows.
- **The pinned lane payload seam — RESOLVED 2026-11-18** (7e4ab14):
  the renderer's row gate reads all three payload shapes (the lane
  payloads' ``total`` beside the free cells' ``e_ground``/
  ``energies_pyscf``); the stored pinned legs render with the
  bookkeeping-form disclosure.  Ti/Sc 4→0 pending each, the twins
  54→51, U14e 26→24.
- **test_symm_tm's SA legs** re-scope to store-anchored verification
  (the AGENTS Tests boundary, user-directed 2026-11-05): fresh CI
  legs stay, the SA verdicts load the committed route-batch stores —
  lands at the file's next touch.
- **The char-selected fold seam** (the m-complete 6p promotion,
  CAS(10,15)): the fold-tag grammar cannot express multi-band folds
  and `closed_electrons` assumes one band per shell; the fix is
  character-selected fold spans with reference-occupation counts —
  a real seam extension when wanted.
- **Review declinations, ticket-40 wave** (2026-11-17): (a)
  `SeriesSpec.seed` stays `str | None` + runtime refusal (NOT
  `Literal["self"]`) — consistency with the declaration-style
  validation (`ReferenceSpec.accel`'s own pattern), the battery pins
  the vocabulary; (b) the `_cap` log-line closure stays nested in
  `_window_columns` — display shaping at its one call site, a module
  level name would widen the interface for a print; (c) the probes'
  pick twins (`probe_u14e_mom.seed_columns`,
  `probe_u14e_series.seed_density`) stay un-retired to the
  production `mom_self_seed` — committed historical instruments whose
  logs are the provenance of what ran (their duplicate BLOCK, a
  real bug, was fixed); (d) the seeded-refusal's surviving
  `<chk>.<seed_basis>.source.chk` cache is CORRECT, not an
  asymmetry to fix — the source is the small-basis run certified
  independently at its own creation (its own post-check); the
  target's refusal does not impugn it, the cache convention applies.

### Disposition candidates (the user's future call)

- **The (6,13) active-space identity — MEASURED 2026-11-25**
  (`probe_u14e_span_compare.py` + log): the three optimized runs —
  the free SA (C1), the pinned@pur on the as-built window, the
  pinned@pur on the @r1 window — converge to the SAME active
  subspace: the span principal cosines 13/13 at 1.0000 across ALL
  three pairs; the pinned ensemble energies equal to 1e-8
  (-50.45150714/3); the ensemble densities agree (the sorted
  profiles <=0.046).  THE READINGS: (i) the optimized SA is window-
  and route-INDEPENDENT — the attractor is the ensemble density's
  dominant 13-NO span (the complete-multiplet ensemble has a unique
  eigenbasis — the SA-anchors-the-gauge fact, now measured at the
  span level); (ii) the window economics are a FROZEN-window story
  — the trilogy's frozen classes rank the windows' CAPTURE, while
  the optimized rows all land the same space and (pinned) the same
  e_tot; (iii) the optimizer MOVES the span massively from the
  as-built input (all 13 principal cosines 0.85-0.96 — no direction
  shared at >0.99), while the @r1 input starts 6/13-in (the
  prior-axis device captures the attractor better — the +140.5
  frozen vs +1065 reading, confirmed); (iv) the optimized partition
  = the L=9 class dominates: the free 17x90, the pinned 16x90+72
  (the Au block's 5th root exactly 72.0 — the block multiplicity's
  boundary; the counts lines' 'other' buckets are the 1e-3 snap
  tolerance losing the 90.000x draw noise, not physics).  THE
  LANE-RETIREMENT INPUT: the routes differ in GAUGE only — the span
  is the invariant; the lane's remaining value is the frozen-class
  readout, not a different optimized space.  (The clause that made
  PT2 gauge-visible in the column gauge is retracted, 2026-11-27:
  the large pinned E(2) was a readout defect -- the autopsy entry.)
- **The static/dynamic correlation trade across the (4,0) ladder —
  MEASURED 2026-11-25** (the V&V table's rows; the frozen dE / the
  SA dE / the static gain / the E(2) / the NET E_tot move): (4,12)
  +923/-42/+965/-342/-253; (4,12)@r1 +149/-42/+191/-342/-88; (4,15)
  +923/-45/+968/-336/-247; (4,16) +923/-45/+968/-334/-247; (6,13)
  +1065/-71/+1136/-282/-337; (6,13)@r1 +140/-71/+211/-282/-56.  THE
  READINGS: (i) the E(2) remainder is nearly constant (~-282..-342
  across ALL the windows) — the trade's dynamic leg barely moves
  (free-route rows; the pinned twins now agree to <= 0.7 mEh where
  the term is closed, 2026-11-27);
  (ii) the STATIC gain carries the window economics (the @r1
  windows pre-capture it frozen: +149 vs +923); (iii) the NET E_tot
  move = the window's missing correlation — the SMALLEST for the
  @r1 windows ((6,13)@r1: -56 vs the (6,13) as-built: -337) — the
  @r1 windows are the most complete; (iv) the ANTICORRELATION (the
  static up => the dynamic down) does NOT hold across windows (the
  E(2) is stable while the static moves) — the trade is a
  within-window conversion (the frozen -> the optimized), not an
  across-window law.

- **The CAS(4,16)/(4,18)/(4,21) correlation-channel ladder — MEASURED
  2026-11-24** (`u14e_vv_4_16_21.py` + log; the compositions as
  addresses/declarations, never code): the (4,16) pair (5f+6d+14p+7s)
  — the frozen classes unchanged (the raw window freezes the 5I at
  +923, the @r1 freezes the 5L at +149), the CASSCF descends ~5 mEh
  deeper than (4,12) (the p/s correlation channels pay), the SAs land
  the (4,12) optima — 8 rows.  THE (4,18)/(4,21) COMPOSITIONS:
  REFUSED BY THE TAIL CERTIFICATES — the 6p channel's radial signal
  is RANK-1 in every available reference (the incumbent Gram's 2nd
  eigenvalue 8.8e-7, the @r density's 2.2e-5 — both tail): the 6p⁶
  trio's MOs are m-separated and the device's m-summing collapses
  their density to rank-1 — a CLOSED-shell channel carries one radial
  by construction; multi-radial 6p needs an OPEN-6p configuration in
  the series (a different reference, not a window change).  The (4,21)
  composition is certificate-measured-out.  16 rows -> 8 + the two
  certificate disclosures.
- **The irrep-occupation approach on uranium — MEASURED OUT
  2026-11-24** (`probe_u14e_irrep_scf.py` + log): PySCF honors
  `irrep_nelec` for ROHF (the per-irrep occupations match the dict
  exactly), but the filling comes out WRONG for uranium's valence —
  the Ag block conflates 6s + 7s + two d functions, and
  block-aufbau drops the diffuse 7s² pair in favor of the valence d
  columns (the occupied characters read s 2.0 / d 3.0 against the
  intended s 4 / d 1; +155 mEh above the trusted rung; converged
  False; three runs reproducible).  The D2h l-span arm (our hook on
  the D2h mol) ALSO fails — the C1-gauge hook mis-fills (p 8 / s 2 /
  d 4) and diverges (E ~2400-2800).  Downstream: the irrep
  reference's MO-direct window reads <L2> 0.005 (closure destroyed
  5.75/8.04 — the distillation necessary, decisively) and its @r1
  device has no 7s pair to distill (the s-dominated occupied list is
  empty — the pair was never formed), reading 55.8 / +261.8 mEh with
  broken closure.  VERDICT: the l-span character-forced C1 rung is
  the necessary reference; the irrep vocabulary structurally cannot
  declare "6s2 6p6 7s2 closed, 5f3 6d1 open" for this valence.  The
  D2h framework keeps its value as the pinned routes' readout
  wrapper only.  The irrep roles themselves (the four jobs + the
  decoded PySCF example) are documented in README.
- **The pinned@pur E(2) "blow-up" — RESOLVED 2026-11-27: an instrument
  defect, not physics** (`probe_u14e_pt2_autopsy.py` + log;
  `logs/u14e_vv_ladder_pinned.log`).  Six pinned@pur SA rows once
  read E(2) -1270..-1580 mEh where the free route reads -282..-342
  ((4,12), (4,15)@r1, (4,16), (4,16)@r1, (6,13), (6,13)@r1).  CAUSE:
  their readout paired the CI vectors of the CONVERGED window with
  the START window's columns (`make_blocks(c_read=)`); the <L2> meter
  was built on the same wrong columns.  THE AUTOPSY (one pinned SA at
  (6,13); E(2) per component, mEh): as the defect built it -1580 (Sr
  norm 15.2, its one-body term 14.3, effective denominator 12 Eh --
  not small denominators); the SAME foreign CI on an ORTHONORMAL set
  -1939 (Sr 15.4) -- orthonormality is not the driver, the foreign CI
  vector is; the same orbitals with the CI re-solved in them -1098
  (Sr 3.2, two virtuals at denominators 1.0 Eh) -- the
  frozen-wrong-term class, i.e. the plain-window frozen casci rows
  (-1054..-1081) are real numbers for a bad active space; the
  corrected object -282.5 per component (certificate 1e-14) = the
  free route.  THE FIX (7274ba5): a state carries the orbitals and
  the CI of ONE solve and is certified at birth; `c_read` no longer
  exists.  THE RE-PRICE (1f71ee9; 8/8 assertions): -341.73, -336.72,
  -334.30, -334.29, -285.50, -285.52 -- each at its free twin ((6,13):
  the pinned row is the component mean incl. the one 5L component at
  -333; its sixteen L=9 components read -282.5).  The two rows the
  old code priced consistently were re-priced on the certified path
  (`logs/u14e_vv_ladder_pinned2.log`; 5/5): (4,12)@r1 -341.74
  (unchanged at the table's precision), (4,15) -336.72 (-336.71
  before).  RETRACTED: the
  earlier mechanism texts ("the Sr class carries a misplaced valence
  orbital"; "gauge-visible at the Eh scale"; "~790 mEh pure gauge").
  The retired probe's log stays as the record of the defect.
- **The (6,13) occ control re-scoping — proposed 2026-11-23** (the
  user's call): the frozen ⁵L ground on the @r1+MOs window uses the
  pair column at **0.53** — the L=8 root mixes 7s²/7s¹/7s⁰
  configurations internally (the class is right; the
  pure-configuration expectation isn't).  The control as written
  (occ ≈ 2) reads RED; re-scope to a disclosure — the column-usage
  jump 3e-03 → 5.3e-01 plus the ⁵L class IS the fix's deliverable —
  and keep closure / definite-L / ΔE as the gates.
- **The lane retirement — PARKED 2026-11-23** (the user's call; the
  gate HELD): the pin-on-pur route reproduces the (4,12) lane
  record at machine precision with no D2h SCF anywhere, and hosts
  what the lane refuses (the (4,15) SA).  The counterweight: the
  (6,13) class is lane-construction-dependent (12..42 pur-as-built
  vs 16×90+72 lane vs 17×72.000000 @r1+MOs) — the lane's (6,13)
  record is a property of the lane construction, worth knowing
  before retiring.  The campaign lane cells stay as the certified
  record either way; any regen is a separate call.
- **The V&V error-detection instruments — proposed 2026-11-22,
  FIRST TWO LANDED 2026-11-22** (the wave's caught-error base: four
  certification-caught math bugs, two latent crashes, two
  presentation bugs, one wrong assertion, one construction defect
  — and the (6,13) shifted-slice meter artifact, caught by the
  two-source comparison against the campaign cell): (a) the
  window-closure preamble — LANDED in the 6_13 leg
  (`window_closure` printed per run; the other legs get it at
  their next legitimate re-run, keeping their committed logs
  byte-reproducible from HEAD); (b) the render-time consistency
  check — LANDED in `render_table` (raises on any row violating
  E_tot = E_method + E(2)); (c) the SLICE-IDENTITY assert — LANDED
  (the leg asserts `key_core`'s ncore == its NCORE constant: the
  shifted-slice class dies at the seam); (d) the pre-flight gate —
  the driver's `--ground-only` pass (minutes) before any full
  chain, catching the latent-crash class at 1% cost (exists);
  (e) the span-health tool — COMMITTED as
  `scripts/probes/probe_u14e_vv_span.py` (the (6,13)
  triangulation + the (4,15) pinned-mechanism instrument, one
  probe, its log committed): the principal-angles comparison of
  the intended vs assembled/projected spans, the closure meters,
  the block-content ranks, and the m-pure-feed refusal.  The user's call on further scope.
- **The (6,13) window rework candidate — the user's declaration
  decision (2026-11-22)**: the as-built (6,13) window's frozen
  ground is the ⁷F (the d⁵ high-spin class; the 7s column
  unoccupied) — the 7s distillation's k=1-of-3 capture (0.27–0.71)
  is the weakest radial in the campaign.  The candidate: declare
  k_s = 2 or 3 for the (6,13) s chain (a CaseSpec/`shell_radials`
  change — new windows, new cells), THEN re-read the frozen
  ground's term.  No code — a declaration change.
- **The NO probe for (4,15) — POSTPONED (the user's debt,
  2026-11-22)**: occupations + NO-group selection
  (`natural.no_group_windows`) from the term-averaged 17-root SA
  density — does a p-type direction survive the selection, with
  what occupation?  A window-quality datum, not a pinned-route
  rescue (the block-ground phenomenon is span physics).
- **The ⁵L-selected-density → NO-window construction — POSTPONED
  (the user's debt, 2026-11-22)**: the root-extraction-by-identity
  machinery (extra roots per pinned block; the meters classify
  ⟨L²⟩≈72, ⟨S²⟩=6) as the tool toward a better initial window — the
  ⁵L-selected roots' density feeding the NO window.  The same
  experiment quantifies the ⁵L's push-down under the block grounds
  (the destabilization measure) when wanted.
- **The (4,15) disjoint partition — the (6,21)-class follow-up —
  NOT NEEDED for now (the user's call, 2026-11-22)**: the second
  validation showed the promoted 14p band is the occupied 6p
  shell's own middle radial (`shell_bands`: 6p = 13p/14p/15p), and
  the pinned route's declared core span contains those slots by
  construction — the pinned (4,15) is the pinned (4,12).  Dropping
  6p from the core would be a different chemistry (the semi-class
  move), not a repair of this experiment.

- **The PT2 wave's declared-out classes** (gated recoveries; each
  carries its measured reason in ticket 45's git history): the uncs'
  lane legs (the lane substrate's basin steer diverges on today's
  box); Ti's 3d4d@ar / 3d4s4p4d@ar lanes (the same-l distilled wall
  — `window_distilled`'s future seam — plus schema-stale stores) and
  the semi lanes (the Löwdin needle on fresh construction); Sc's
  semi lanes (the same Löwdin class) + 5 substrate drifts (the regen
  loop, gated); U14e's 11 refusals (the raw-AO windows' re-solve
  does not reconstruct; the 7s SA cells past e_noise; the lane legs
  on real reconstruction failures).

- **The U14e prior optimizer follow-on — EXECUTED 2026-11-18**
  (bb39a70 + the close commit; the prior-anchor-first scenario's
  dated claim re-tested by `--full`, both prior windows, 8 legs
  each): the escape class is UNIVERSAL — the free single-root
  CASSCF descends ~190 mEh below its CI anchor BY BREAKING ⟨L²⟩
  (65.0–65.2, [H,Lz] 5e+00–1.3e+01, cycle-capped 1e-5/200); SA
  holds the cure: −42.0 m at ⟨L²⟩ 72.0000 (4e-7), commutators
  1e-4, converged, on BOTH priors; the frozen pur class certifies
  machine-zero exact ⁵L under optimization (1e-13/1e-14, both
  objectives, both windows).  THE DEEP FINDING: the optimizer's
  landscape is window-INDEPENDENT — @r1 and @s1 free SA read
  −50.4238/−1.1418 identically (the escape floors coincide at
  −50.4298 too): the 5.6 mEh CI-level @s1 advantage was a
  span-capture artifact of the frozen-orbital read; the relaxation
  channel dominates.  The basis decision's evidence is complete:
  both priors optimization-stable under SA, neither prior's
  identity surviving to the optimizer level.  The twins' adoption
  decision carries this evidence (the same un-gating).
- **The swap-basin optimizer class re-run gauge-anchored** (the
  Stage-2 leftover): if the campaign wants in-window solutions on the
  twins, the optimizer legs re-run seeded -- rides the Stage-3
  adoption decision that un-gates `run_casscf`.
- **`_COMM_CAP` (cells) = 6000** skips the dense [H,L] meters past
  6000 dets (the U14e (6,13) sector's em-dashes): this box could
  afford the raise + an ~80-min targeted re-run of the 24 cached
  legs if the certification is ever wanted there.
- The U14e off-basin `moc pur` starts (the multi-solution class,
  honest meters in the payloads) + the `d14p sa moc_frozen pol`
  27.7-mEh reconstruction outlier (the frozen-moc gauge's
  analysis-ruler artifact class).
- The (4,19)/(4,21)-class lattice rows: the ndet² analysis ceiling —
  the same class as `_COMM_CAP`.
- **A data-derived note-class for the poverty rows** (the CASCI
  review's leftover, 2026-11-17): C/Si's D-window CASCI rows read
  +221…+615 mEh with exact closure certificates but an EMPTY note
  column — the poverty signal is invisible there.  A renderer
  note-class (magnitude-derived, keyed on CASCI + ΔE threshold) would
  surface it; the notes are currently identity-derived only.  Rides
  `build_methods_tables.py`'s next touch; the user's call.
- **The Stage-3 adoption decision** (the user's call, the measured
  basis in hand 2026-11-16; RE-SCOPED PROBE-FIRST 2026-11-18 — the
  U14e optimizer wave's evidence: the escape class is universal and
  the optimizer landscape is window-independent, so the twins' full
  51+51-row grids would mostly re-measure a characterized class):
  run a scoped probe batch first (@s2 pur, ~4–6 legs per twin:
  casscf/sa × maskless/frozen), then scale or decline the full grid
  ON THE EVIDENCE.  `subspace_rohf`'s certificate says the
  @s2 windows are MEAN-FIELD-SOUND everywhere (the confinement costs
  1.3–2.2 mEh at the ROHF level, every span, every filling, both
  twins), @s1 is filling-dependent (3.3–38 mEh), @r1 is
  head-filling-only (+1.4/+1.9 for its own head, up to +199 mEh for
  the d-richer members).  Un-gating `run_casscf` on the @s2 windows
  is certified mean-field ground; the queued CASCI review EXECUTED
  2026-11-17 (no analysis bugs found) — the execution waits on the
  probe batch's go.  The swap-basin re-run (below) bundles with this
  decision (the same un-gating).

### Watchlist

- The latent refusal on `--fresh`: the incumbent distilled FRAC
  cells predate the distillation declaration — their caches mask the
  construction refusal.
- The adjudication watch: the pinned solver's roots are proposals
  our exact H corrects (self-verification, not solver independence)
  — a shared-convention error is covered by the dense commutator
  certificates + the ion anchor; revisit if a block class appears
  where adjudication is unaffordable AND the solver skips.
- Revisit conditions: the lane-family extraction from cells (if the
  radial legs grow it); the m-replication/assembly duplication
  (windows vs basis — extract at the next legitimate change, bitwise
  A/B); the per-atom prepared-basis bundle (a second consumer or a
  fifth field); the prior_bands-defaults rule hoist from cells (when
  next touched); the projected-Gram alternative (if the per-radial
  spectrum shows truncated-ranking distortions); the D2h determinism
  `--verify` mode (a second regen cycle); the D2h mol shared builder
  (a third frozen reference); `run_frac_occ(dm0=)`'s disposition —
  kept as the reference family's symmetric explicit-seed form, its
  one consumer was the 2d measurement attempt, revisit if no
  frac-bearing atom appears by the next reference-layer wave.

### Settled (standing facts, one line each; git history the archive)

- **The U14e band-identity poverty is the compact radial, and the
  prior cures it** (tickets 39 + 40, 2026-11-17): the incumbent
  10f+13d window captures only 82-83% of the head's f-valence density
  (the full 9f-11f family 96-97%; 13d 92% vs 98.5%); the @r1 prior
  recovers 84% of the +923 mEh poverty at CASCI (+147) with the exact
  ⁵L on pur, and the @s1 series average is the BEST prior
  (+141.3/+143.7 — ~5.6 mEh below @r1, ~85% cured, exact ⁵L on pur).
  The census's member refusal was the WEIGHT HOOK's vocabulary, not
  the landscape: the character-only premium ties ~12 pure d/f valence
  virtuals + the tight tail at premium ≈ 2.000 (a combinatorial tie
  broken by column index); the self-seeded MOM lands both members
  (f⁴ −50.31571501 l3:4.000; f²d² −50.34788201 under two
  accelerators — plain newton is draw-sensitive there, diis+newton
  lands the anchor in three fresh draws; on f⁴ diis+newton RELAXES
  INTO THE HEAD).  The seed pick rule is ε-first within the majority
  floor, never weight-ranked (the tight-primitive trap).  @s2 is
  construction-refused: the ROHF 5f radial is configuration-inert —
  the series average carries no second radial (k=2 eigenvalue 5.9e-5
  relative), the d shell never asked (k uniform by the tag grammar).
  The basis carries NO near-dependence tail (minima 2.2e-6 s … 7.8e-3
  f; the conditioning cut NO-GO by the twins' precedent; the f block
  is the BEST-conditioned, κ 428).  `run_rohf` deletes a refused
  draw's chkfile (the gate's or the post-check's) — a refused
  draw never poisons the trusted cache.

- **The U14e campaign COMPLETE** (tickets 36A + 38): the free-mask
  SA legs certify the ground term END-TO-END — ⁵L, ⟨L²⟩ = 72.00
  exact, closure + commutators at machine zero, bookkeeping ~1e-3
  mEh; single-root legs break symmetry as declared; the input CI
  grounds are ⁵I/⁵F (the ~1 Eh gain IS the ⁵L relaxation); best
  energies (6,13) −50.456 < (4,15) −50.434 < (4,12) −50.431; the
  retired tier windows' −50.386 validates the promotion; 8
  unconverged legs ACCEPTED (the user, 2026-11-16 — the no-free-mask
  crawl; the physics line rests on the certified legs).
- **The honest prior table** (Sc/Ti/C/Si, the @r/@s windows): NO
  universal k=1 winner (atom-dependent, ~11–100 mEh spread) — the
  SA-over-ladder optimization is the only road to a dominant k=1;
  @s2's dominance is variational headroom at double budget (every
  comparison row carries the CAS shape); C's @s2 REFUSES (the
  rank-1 series average) while Si's stays (a real second radial);
  the committed probe logs are the flipped-metric era — the module's
  numbers are the honest successors.
- **The uncontracted-twins reference disease + the seeded repair**:
  fresh SCF paths diverge on the unc bases — the CHARACTER HOOK, not
  the SCF (metric-amplified virtuals carrying span weights ~1 +
  shell-counter labels selecting AO slivers); the seeded-MOM repair
  (the projected 6-31G source, source-vocabulary hooks, the
  character-of-occupied post-check) lands the declared anchors
  EXACTLY — the product path since Stage 1c; both twins' trusted
  records verified to machine epsilon.
- **Stage 2 completed the twins' prior axis** (2026-11-16): the
  members' recipe is `diis+newton` — Newton alone crosses into the
  deeper sibling basin (the d3→d2s1 fall, −128 mEh, reproduced
  exactly), plain DIIS holds the configuration but never converges;
  `projected_source` must thread the declaring ref's accelerator
  (the member sources are the DIIS-oscillation class). The @s legs
  ran 36/36 per twin, zero refusals; the contracted prior verdict
  holds at the basis limit (the series mean loses to the head's own
  block, @s2 buys the double-budget headroom). The UHF anchor row
  re-landed the censused basin through `uhf_seed_density` (the
  pre-repair record was 11.4 mEh shallow). The frac class
  MEASURED-REFUSED on the twins (both accelerators diverge seeded —
  the 43-AO block feeds the amplified subspace); the averaged core
  is their spherical field.
- **Stage 3 measured the mean-field certificate** (2026-11-16):
  `subspace_rohf` confined-valence solves on both twins, 12 rows
  each — the whole-space legs exact to ≤ 2e-6 mEh (the degenerate
  limit reproduces the anchors), @s2's confinement cost 1.3–2.2 mEh
  everywhere (mean-field-sound), @s1 3.3–38 (filling-dependent),
  @r1 1.4–199 (head-only — excellent for its own head, catastrophic
  for the d-richer members); within-@s monotone everywhere, the
  cross-family r1/s1 nesting measured NOT nested (different
  eigenproblems).
- **The subspace-ROHF's five measured defects** (the bench's
  evidence chain, 2026-11-16): the seed density must reach the
  kernel SPIN-SEPARATED (a summed dm is split (d/2, d/2) — the
  open-shell asymmetry lost, the filling off-basin from cycle one);
  the veff shim returns the TRUE per-spin pair through the substrate
  (a symmetric pair from the summed density carries the
  closed-shell exchange formula, Eh-scale wrong on opens); fills
  size COLUMNS not electrons; NO block-index pools (the subspace
  Fock's eigvectors interleave across the family/window layout — the
  4s ranks 9th globally, inside the family's index range); the MOM
  assignment runs on the CLEAN spin-summed Fock with the anchors in
  the clean gauge (the kernel's Roothaan effective Fock is
  MIRROR-DEGENERATE in the low spectrum — the pairs split at 1e-7
  and flip gauge every cycle; identical Fock, the picks flip from
  exact to garbage).
- **The twins' confined_core_orbitals floor**: the closed-core
  Roothaan on the 65-column family span floors at ~5e-10 past its
  own 100-cycle budget (the conditioning class) — the ion
  cross-check runs on the carbon-class span (exact), the twins'
  refusal recorded.
- **The folded-field records are draw-documented single draws** (no
  replay; the post-kernel converged-filling certificate guards
  production now); the +400 mEh folded offset is basis-stable
  construction physics (the missing core-valence self-consistency),
  not span bloat; the certified numbers stand as measurements.
- **The D2h substrate**: NEVER symmetrize a degenerate-gauge seed
  (`symmetrize_orb` amputates the open MO, `symmetrize_space`
  scrambles the pairing); the solver's irrep blocks erase the gauge;
  pinned-cell windows are D2h-born sources only (any future intake
  needs an irrep-purity check at the seam).
- **The basin zoo is census-resolved** (the twins + U14e): the
  deepest basins are the adopted anchors; plain aufbau can land 10+
  mEh high at the right filling — the seed, not the gate, is the
  cure when gate cycling diverges; the Newton path freezes
  occupations (the occ_ok certificate refuses converged wrong
  configurations).

## Examples (worked records in `elements/`)

- `elements/C/`, `elements/Si/` — the p² campaign.
- `elements/Sc/`, `elements/Ti/` — the transition-metal campaign
  (the radial axis, the d² lattice, the prior windows).
- `elements/Sc-unc/`, `elements/Ti-unc/` — the uncontracted twins
  (the folded/lane records, the prior axis, the seeded references).
- `elements/U14e/` — the ECP atom: the certified ⁵L record.
- `elements/U/` — the frozen pre-campaign ECP workdir;
  `elements/ecp/` — the ECP seam's certification workdir.

## Architecture (edge-verified)

    pure math      system · detbasis · secondquant · natural
    substrate      windows (the hub)
    1-particle     angular (uses windows' shell_slots)
    assembly       operators (both PySCF boundaries) · irreps (the
                   D2h dictionary; the symm boundary)
    engines        reference · core · casscf
    analysis       analysis · stages · scans
    declaration    cases -> cells

`natural`'s mathematics is pure at module level (numpy + secondquant
only); its one driver (`no_group_windows`) lazily imports the
assembly layers. No cycles, no upward edges. Division of labor over
windows: declaration NAMES (labels), windows BUILDS (the kinds +
substrate builders), cells BINDS (the start axis picks the builder).
