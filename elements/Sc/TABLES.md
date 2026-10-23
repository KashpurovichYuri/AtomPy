# Sc — scandium (d¹)

cc-pVDZ (Balabanov–Peterson); reference **ROHF −759.736255 Eh**
(3d¹4s²; M_S = 1/2; S can be half-integer — snap to the half-grid).
pur substrate: the raw [Ar] core polarizes at the percent level
(anisotropy 1.3e-2 rel, [H,Lz] 2.7e-2); pur closes it.

Provenance: `results/sc_ladder.md` + `results/*_ladder.md`, the Sc
NO pilot `results/tm_no.log`, `tests/test_ladder.py` Sc block.
Cross-atom verdicts: README §7.3–7.5.

## Term content

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(1,5) / {3d} | Ms=½ | 5 | ²D: 5 |
| CAS(3,6) / {3d,4s} | Ms=½ | 90 | ²D: 20, ⁴F: 14, ²F: 14, ²G: 18, ⁴P: 6, ²P: 6, ²S: 1, ²H: 11 |

Sc CAS(1,5) is the Sc⁺-like d¹ space on the [Ar]+4s² closed core.

## Multiplet spectrum (CAS(1,5), our fold, exact ≤1e-9)

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ²D | 5 | −0.342549 | — | −2…2 |

## Window lattice — exact symmetry vs raw MOs at matched cost

| window | CAS | dets | AO, Eh | ⟨L²⟩ | ROHF MO, Eh | ⟨L²⟩ | NO |
|---|---|---|---|---|---|---|---|
| 3d | (3,5) | 50 | −759.436 | 12.000 | — | — | — |
| 3d4s | (3,6) | 90 | −759.729 | 6.000 | −759.736 | 0.000 | −759.731 |
| 3d4d | (3,10) | 450 | −759.555 | 12.000 | — | — | — |
| 3d4s4p | (3,9) | 324 | −759.760 | 6.000 | — | — | −759.762 |
| 3d4s4d | (3,11) | 605 | −759.739 | 6.000 | −759.745 | 0.008 | — |
| 3d4s4p4d | (3,14) | 1,274 | **−759.765** | 6.000 | −759.745 | 0.143 | **−759.765** |
| 3d4s4p4d5d | (3,19) | 3,249 | −759.765 | 6.000 | — | — | — |
| 3d4s4p4d4f | (3,21) | 4,410 | **−759.769** | 6.000 | — | — | — |

## Treatment ladder (ΔE vs ROHF, mEh)

| rung | ΔE | note |
|---|---|---|
| CASCI, averaged fold | +6.7 | exact |
| CASCI, raw core fold | +6.7 | comm 2.8e-2 — the raw-core detector |
| CASCI, frac reference | +36.8 | frac septet quasi-degeneracy — never pin the ordering |
| CASSCF, free | **−31.5** | escape to a ²S-like basin ({s,p,d}-mixed window) |
| CASSCF, window frozen | +6.6 | |
| CASSCF, radial restriction | −6.6 | single root, chaotic; ²D→²P relabeling across draws |
| SA-CASSCF, free | −13.3 | swap basin (totals only; the fold reads −730) |
| SA-CASSCF, frozen, purified core | +7.3 | exact (pur machine-exact) |
| SA-CASSCF, radial, purified core | +0.0 | flat-gauge floor |

Scan anatomy: the m-symmetric radial curve only RISES (3d→4d, 4s→5s/6s
all harmful); member tilts descend (−0.3/−0.9 mEh) — the JT class.
Ticket-24 amendment: the +25 mEh radial ceiling is a WINDOW statement —
the grown lattice crosses below the escaped CASSCF value (4f radial
−3.06 mEh) at exact symmetry.

## NO groups (ticket 22)

Term-average of the ⁵L... (Sc: the largest-window ground term) → groups
4s(1) / 3d(5) / 4p(3) / [3d+4d](5); NO beats the shell lattice at every
equal-det count measured (90/324 dets) and the ROHF-MO column at the
top pairs; every group window at ~1e-13 certificate, exact ⟨L²⟩.

## Open questions

- Semicore (3p) treatment — user-deprioritized.
- The frac-occ septet ordering (SCF-chaos-dependent — never pin).
- Spin–orbit splitting of the ²D multiplet.
