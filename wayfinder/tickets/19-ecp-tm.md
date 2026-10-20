# Ticket 19: ECP Sc/Ti — gated references, free-ion tables, the pur family at ncore=4

Status: RESOLVED (code zero-change; tests + measurements + this record)

Origin: ticket 18's deferred TM half. All numbers below measured in this
env (PySCF 2.13, ccECP + ccecpccpvdz, chkfiles in results/).

## The references (basin census, fact-21a under ECP)

- **Ti** (12 e, spin=2): NO aufbau runaway (no [Ar] to hide in; the ECP
  core is 3s3p-only), but the member-gauge basin zoo transfers with the
  default REVERSED: forced minao/atom land the good basin
  E = −57.60649802 (build_full_mo window defect 1.4e-8); huckel AND the
  unforced default land the metastable +11.1 mEh basin (defect 1.6e-2 —
  the same numbers as all-electron Ti's +11.2 mEh / 1.6e-2). Gate with
  occ_e_ref = −57.60649802.
- **Sc** (11 e, spin=1): THREE classes. Unforced minao/vsap/dk: physical
  3d¹4s² at −46.12083302. Unforced atom/huckel: a 4s²4p¹ runaway
  (+104.7 mEh — the open electron goes to 4p). FORCED from every guess
  (minao/atom/huckel/vsap): −46.12052897, window defect ~1e-2. The
  forcing hook itself shifts the basin (+0.30 mEh vs the unforced
  minimum, same occupations) — the forced basin is chosen as THE
  reference because it is reachable from every guess (reproducibility
  beats 0.3 mEh; the shift moves all rungs together). occ_e_ref =
  −46.12052897.
- Sc's ~1e-2 input-window defect is INTRINSIC to this reference (every
  basin carries it): the fact-17 mechanism amplified — the explicit
  3s3p core sits directly adjacent to the {3d,4s} window (no deep core
  to dilute; ⟨4s|3s-MO⟩ large). Ti's good basin: 1.4e-8.

## The tables (raw shell window, averaged fold; our folded H)

- Sc CAS(3,6) Ms=½, 90 dets: ²D×5 ground, then ⁴F×7, ⁴P×3, ²F×7, ²D×5,
  ²G×9, ²P×3, ²S … — textbook, exactly degenerate within each term,
  integer ⟨L²⟩/⟨S²⟩; ||[H,Lz]|| = 1.5e-12.
- Ti CAS(4,6) Ms=1, 120 dets: ³F×7 ground (⁵F×7 at +4.2 mEh), ³P, ⁵P,
  ³G, ³H, ³D, ⁵D, ³H, … the full free-ion table; ||[H,Lz]|| = 5.7e-13.
- Fold bookkeeping: max|dE| vs PySCF ~1e-3 (Sc) / 2.3e-3 (Ti) on CASCI
  rows — the avg-vs-raw fold difference on a percent-polarized ECP core
  (casci_raw ground shifts −1.0/−2.3 mEh vs casci_avg).

## The ladder (fingerprints; chaos discipline — branch, don't pin)

- Escape present on both: Sc virtmix/shell/pyscf −51.6 mEh vs avg
  (⟨L²⟩ → 0.05, window s/p/d-mixed); Ti −57.3 mEh. NEW FLAVOR on Ti:
  casscf_shell grabbed 4f (window p²d¹·³f⁰·⁷, −0.88 Eh, ||[H,L²]|| = 11)
  — the "complete window still escapes" attractor reaching into the
  ccECP 4f functions.
- Sc sa = the 4s↔3p swap basin (window p¹d², ⟨L²⟩ = 6.002, fold
  illusion −0.83 Eh) — survives ECP exactly as on AE. Ti sa landed the
  HEALTHY window this draw (comm 1.1e-5, max|dE| 4.5e-8) — free SA is
  bimodal under ECP too; sa_pur took the swapped/escaped basin (−0.88
  Eh, p¹d³ window). Branch on window weights, never pin.
- casscf_radial / casscf_radial_pur: single-root gauge chaos as on AE
  (⟨L²⟩ 4.08 / 7.34 on Sc/Ti radial_pur).
- Ti casscf_relax = 6e-8 (the good-basin substrate class, mirrors
  all-electron Ti's 2.3e-8); Sc relax carries its intrinsic 1e-2.

## The pur family at ncore=4 (the fact-20 deliverable under ECP)

- purified_core_orbitals on {3s²3p⁶}: contamination max 5.7e-5 (Sc) /
  7.0e-5 (Ti) — foreign-l tails only.
- Closure certificate: both input windows (build_full_mo, radial_basis)
  close at 4.9e-15 (Sc, vs 1.4e-2 on the canonical core) and 7.2e-15
  (Ti, vs 9.5e-8).
- relax_pur / sa_relax_pur: MACHINE-EXACT on both atoms (1.5e-12 /
  5.3e-13; sa max|dE| 2.8e-7 / 3.6e-8; the pur-fold removes even the
  1e-3 fold-bookkeeping class — the purified density is exactly
  m-scalar, so avg fold == raw fold). sa_radial_pur at its floor
  (2.2e-6 / 2.0e-6); free sa_pur bimodal (above).

## Test lessons (recorded for the next atom)

- The per-root labeling MUST go through analysis.analyze / V†OV with the
  ROOT index on the columns. An einsum with the root index on the wrong
  axis returns garbage that is REPRODUCIBLE and plausible-looking on
  exactly-degenerate multiplets (Sc's ²D×5 read as ⟨L²⟩ ≈ 13.9, ⟨S²⟩ ≈
  2.0); it cost a full diagnosis loop before the subscripts were read
  carefully. Tests now pin the library path.
- Tests added: references/basins, multiplet tables, pur closure
  (tests/test_ecp.py; suite 99). Ladder fingerprints live here + the
  results/ artifacts (re-running optimizers in tests buys no
  determinism).
