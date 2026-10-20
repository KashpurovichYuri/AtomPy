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
- **D. select_diffuse()**: labels lie on uncontracted bases (measured:
  "5f" = the exp-36.8 primitive; valence f = the diffuse end). First-
  class exponent-aware selector + fact-1 amendment.
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

## 3. Pressure points found during U (append-only; empty now)

-

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
