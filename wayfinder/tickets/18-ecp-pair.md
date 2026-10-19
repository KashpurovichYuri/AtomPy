# Ticket 18: ECP support — ccECP pairing, wiring, p² validation (C, Si)

Status: OPEN
Origin: map.md "Post-wave survey" ECP design (pre-verified seam) + the
2026-10-19 pre-flight (all numbers below are MEASURED in this env, PySCF
2.13 bundled sets — no new dependencies).

## Pre-flight facts (measured; do not re-derive)

1. PAIRING: the only bundled PP family covering our atoms is **ccECP**
   with its cc-pVDZ-quality basis: `basis="ccecpccpvdz"`, `ecp="ccecp"`.
   Measured rejections: cc-pVDZ-PP (bundled .dat skips Si, Sc, Ti
   entirely — Ga onwards plus Cu/Zn/Y/Zr+ only); bfd (`bfdvdz` lacks
   Sc/Ti basis); def2-ecp does not parse in 2.13 and def2 is
   all-electron for H–Kr anyway; basis-set-exchange is NOT installed.
   ccECP is correlation-consistent (constructed to track all-electron
   cc-pVnZ) — the right physics pairing for term structure.
2. STRUCTURE (mol.nelectron is ECP-honest automatically):
   - C: 4 e (He-core in ECP), ncore=1 (the explicit 2s²), window {2p},
     CAS(2,3). ROHF converges; E = −5.31427464 Eh.
   - Si: 4 e (Ne-core), ncore=0 (!), window {3p}, CAS(2,3). The first
     ncore=0 case — the fold has an empty core (veff_core must be 0;
     verify, don't assume).
   - Sc/Ti: 11/12 e (Ne-core; 3s3p EXPLICIT valence), ncore=3,
     window {3d,4s} — same active space as all-electron, the [Ar] core
     replaced by ECP + explicit 3s²3p⁶. (TM execution = ticket 19.)
3. AO labels truthful and spherical under ccECP (C: 2s,3s,2px,2py,2pz,
   3p...; Si: 3s,4s,3p,4p...) — fact-1 prefix matching works; m-sets
   complete; lm_sectors decomposition expected to transfer (verify
   fused shells via bas_nctr as usual).
4. SMOKE (C/ccECP, Ms=1, raw 2p window, averaged fold): ⟨S²⟩=⟨L²⟩=2.0000
   exactly, ³P ×3, OUR-fold split 0.00e+00, ||[H,Lz]|| = 6.06e-17 —
   machine-exact, identical to the all-electron control (2.67e-17).
   RAW fold: split 8.02e-05 Eh, ||[H,Lz]|| = 1.09e-04 (all-electron
   raw: 4.37e-05 / 6.06e-06) — the polarized-core detector fires
   LOUDER under ECP (the explicit 2s² sits closer to the open shell);
   record as the new fingerprint scale. ||[H,L²]|| trivially ~0 in the
   complete-irrep Ms=1 sector (fact 12).
5. L operators are pure AO objects — ECP enters only through h1
   (get_hcore() picks up int1e_ecp automatically at BOTH operators.py
   sites: `active_space_integrals` :151 and `full_space_integrals`
   :277). No angular code changes needed.

## Scope (this ticket)

A. `AtomConfig.ecp: str | None = None` + `build_mol` wiring
   (`mol.ecp = cfg.ecp` before build; extend the docstring: ECP honesty
   — mol.nelectron is valence-only, NEVER recompute from Z (the
   ncore_of contract already enforces it)). Update `label` to include
   the ecp tag when set.
B. ncore=0 path (Si): `active_space_integrals`/fold with an empty core
   — measure veff_core == 0 (or skip the fold entirely) and that
   `build_full_mo(mf, [], act)` works; fix whatever the empty edge
   breaks (defensively: CASCI with ncore=0 is already exercised
   all-electron-wise by the core-explicit test? — measure).
C. Pipeline + ladder on C and Si with ECP (CAS(2,3), both Ms sectors):
   - term tables match all-electron ORDERING (³P < ¹D < ¹S) with exact
     ⟨L²⟩ = 2/6/0 and machine-exact commutators (averaged fold);
   - term SPACINGS vs all-electron cc-pVDZ: measure and report (ccECP
     should track cc-pVDZ at the mEh level — record the actual
     deviations, don't assert them);
   - ladder smoke: casci/sa rows on Si-ncore=0 (the fold-free mode).
D. Tests (`tests/test_ecp.py`, slow tier): wiring (nelectron, labels,
   spherical); C-ecp ³P exact closure + detector scale (raw 1.1e-4
   class, branch not pin — chaos discipline); Si-ecp ncore=0 fold
   identity + ³P/¹D/¹S table; spacings-vs-all-electron bands with the
   measured values in comments.
E. Docs: AGENTS fact (ECP: the one measured fingerprint scale + the
   ncore=0 mode + pointers), README section (the ECP reference mode:
   what moves [absolute energies, E_ref] and what doesn't [term
   structure, closure classes]), map update.

## Out of scope (ticket 19)

Sc/Ti ECP: ReferenceSpec re-measurement (E_ref basins are new under
ECP), the 3s3p explicit-core purified story (fact 20 with ncore=3),
TM ladder fingerprints, escape/swap basins under the ECP field.

## Acceptance

- Full suite green; `pytest -m fast` unaffected (ECP tests are slow
  tier). Zero interface change outside AtomConfig.ecp (additive).
- Term tables: ordering + L eigenvalues exact; spacings reported with
  measured numbers.
- No new dependency (bundled ccECP only).
