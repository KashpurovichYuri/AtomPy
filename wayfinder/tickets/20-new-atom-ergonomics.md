# Ticket 20: new-atom ergonomics — setup, run, read

Status: OPEN (write-only for now). EXTENDED DURING THE URANIUM CAMPAIGN:
every pressure point found while running U is appended to §3 the moment
it recurs twice. RESOLVED AFTER U: one deepening pass over the full
inventory. Do not execute during the campaign except items the campaign
literally cannot proceed without (scope split in §4).

Origin: user request 2026-10-20, after three new atoms in two days
(Sc-ecp, Ti-ecp, U-preflight) each re-hand-rolled the same rituals.

## 1. Destination

Adding and running an atom is: **one data row of setup → one serial run
→ readable results**, and problems surface loudly at the step where they
occur (not five columns later, not in a /tmp script). "We often would
need some corrections anyway" — accepted: the run must make the
corrections *findable*.

## 2. Items (each with the session evidence)

- **A. CaseSpec + run_case**: the atom definition is scattered (cfg +
  ECP/basis strings + spans + ReferenceSpec + sa_weights + frac params
  across six call sites). One dataclass row per atom; one serial
  orchestrator: reference → sector plan → census stages → ladder → all
  artifacts under one dir (`results/u_14e/` = a workdir tag, not a
  directory scheme). Precedent: `ladder._registry` (new rung = new
  entry, report cannot drop a rung).
- **B. basin_census()**: the basin hunt was hand-rolled for Ti-AE,
  Sc-ecp, Ti-ecp, and U (guesses x forcing x fingerprints x window
  defect, in /tmp scripts). One driver: E / converged / occupation
  characters / window-defect per guess; the SCF-methodology matrix
  (DIIS+damp+level-shift, Newton, two-stage, dm0-from-CASSCF round-trip)
  belongs in it.
- **C. occupation_characters()**: span-weight fingerprinting re-written
  three times this session. One helper, shared by censuses and reports.
- **D. select_diffuse()**: RESOLVED by ticket 21 part 2 (2026-10-21):
  `primitive_chain` + `select_primitives` (exponent walk, band/count
  semantics, contracted shells raise with the label/MO pointer; the kin
  proxy retired). Remaining ergonomics gap for this item: nothing.
- **E. sector_plan()**: ndet = C(norb,na)*C(norb,nb) and GB = ndet^2*8*nops
  computed by hand all session. Print the plan (sectors, dets, dense-op
  GB estimate + threshold warning) before any 30-minute burn.
- **F. per-rung health line**: the known failure classes (fact-13 ortho,
  fact-17 input-window baseline, escape/swap windows, fold bookkeeping)
  require holding five facts in head while reading five columns. One
  verdict line per rung with fact pointers (summary line, not a
  framework).
- **G. small enablers**: custom basis/ECP strings through AtomConfig
  (+ label tag); parity operator + degree mark + L-letters to L<=11;
  open-count forcing in ReferenceSpec if the U census shows f3d1-vs-f4
  steering is needed (measure first).

## 3. Pressure points found during U (append-only)

- Census scripts re-implemented in /tmp for Ti-AE, Sc-ecp, Ti-ecp, U
  (guesses x forcing x methods x fingerprints) -> item B, confirmed.
- **Canonical-order assumption, three strikes in one campaign**: (1) the
  run script's `range(ncore)` core (max|dE| 14.5 Eh), (2) the analysis
  call relying on `active_space_integrals`' default c_core (13.8 Eh),
  (3) `purified_core_orbitals`' own slice (the pur construction crashed).
  All library defaults now route through `occupied_core_idx`; the
  pressure point: "canonical order = occupied block" was never a checked
  invariant anywhere. CaseSpec should own the occupied-block semantics.
- **Non-orthonormal windows**: two overlapping f primitives are not
  S-orthonormal; raw-column analysis silently breaks every product
  operator (L^2 eigenvalues 172.98, non-quantized -- fact 4). The
  analysis must use the build_full_mo-orthonormalized window;
  `build_operator_parts` could GUARD (assert S-orthonormality of c_act).
- **Comparison arithmetic traps**: active-space eigenvalues vs PySCF
  totals (needs + core_energy); sorted-vs-per-root alignment. Two bugs
  in one campaign script. sector_plan (item E) should print the
  bookkeeping contract next to the sizes.
- **slater_condon_matrix is an O(ndet^2) pure-Python double loop**:
  18,411 dets (tier-1b Ms=1) is hours-to-days; 29k impossible. The
  single blocker for sector depth on U. Vectorization (grouped by
  excitation class) is now a campaign-critical item; the secondquant
  witness + existing pins make it safe to build.
- **Window convergence is itself an experiment axis on U**: the
  f3d1/f4/f2d2 ordering flips with the primitive selection (1f window:
  f4 5I ground; 2f window: f2d2 5L ground by 6.3 mEh; physical answer:
  f3d1 5L). Near-degenerate manifolds need a systematic window scan --
  a new item: window_scan (add primitives by measured open-MO weight
  until the low-spectrum ordering stabilizes).
- **Parity's l-purity invariant needed a percent-class tolerance**: the
  core-orthogonalized 6d window column carries 0.8% of 7s (fact-17
  class); <P> reads 0.98 and the pur substrate restores 1.0.
- **PySCF's Davidson straggles inside degenerate multiplets** (found
  2026-10-21, part 2): nroots=48 on 10,626 dets, up to 38 mEh errors
  WITHIN exact 13-17-fold clusters (single-root tolerance tests pass
  while cluster interiors drift). The exact fix: `pspace_size = ndet`
  (pyscf's full-diagonalization branch, needs norb < 64) -- now
  `run_casci(fci_opts=...)`; both-sides-exact max|dE| = 7e-15. Health
  lines (item F) must distinguish "fold difference" from "solver
  residual" before quoting any dE.


## 4. Scope split with the U campaign

The U campaign (u_14e) implements its MINIMAL enablers itself (G's
custom-basis wiring, select_diffuse, parity, L-letters) because it
cannot run without them. Ticket 20 afterwards absorbs and generalizes:
promotes the campaign's inline hacks into the CaseSpec vocabulary,
lands A/B/C/E/F, and deduplicates whatever §3 recorded. Tests per the
working agreement; the campaign's artifacts (scripts/run_u.py,
results/u_14e/) become the regression examples.

## 5. Acceptance (post-U resolution)

- A new atom = one CaseSpec row + `run_case`; no new orchestration code.
- The serial log answers, in order: did the reference land (basin,
  gate), is the sector plan sane (dets/GB), are the stages healthy
  (ortho/purity/closure/commutators/fold classes with fact pointers),
  where are the artifacts.
- §3 empty-or-resolved: every pressure point recorded during U maps to
  an item in §2 or a new one, and none remain unhandled.


## 6. Axes-systematization acceptance (added 2026-10-24, from the
ticket-22/23 waves + the user's window x treatment direction)

- `Window.columns(mol, mf)` resolver: one first-class object per window
  spec (AO shell, primitive band, MO span, NO group); NO windows carry
  their dependency (a parent CASCI for the term-averaged density).
- Axes runner: window x treatment x start, per-cell budget knobs
  (mc_opts pass-through exists; the user wants longer runs later) --
  fills each atom's elements T3 product table.
- Completion runs (the measured gaps): U SA-24, relax-24, SA-(6,13);
  Sc/Ti NO-start ladders (SA + relax); then the full product.
- `--fresh` / chkfile key versioning (cached artifacts trusted today).
- Intent manifest (precursor: the six13 printed intent header).
- Full-matrix slice guard at library level (ncas from mc.ncas;
  [ncore:ncore+ncas] on outputs too).
- Needle linear-dependence helper (the relative-threshold path threw
  one unreproduced ValueError session; map.md survey 2026-10-24).
- Promote `closed_shell_core_energy` from the scan scripts to the
  library.
