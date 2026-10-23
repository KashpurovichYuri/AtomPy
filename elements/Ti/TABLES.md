# Ti — titanium (d²)

cc-pVDZ (Balabanov–Peterson); reference **ROHF −848.406489 Eh** — the
gated good-basin reference: default ROHF is an aufbau runaway (+123.07
mEh) and the SCF is multi-basin even with character-forced occupations
(ReferenceSpec gate; ticket 12/21a). Standard sector M_S = 1. Ground
term ³F (first L = 3 ground; first ¹G in the Ms=0 spectrum).

Provenance: `results/ti_ladder.md`, `results/*_ladder.md` (full 17-rung
registry), `tests/test_ladder.py` Ti block, `TI_EXACT_TERMS` in
`tests/conftest.py`. Cross-atom verdicts: README §7.3–7.5.

## Term content (what each sector must contain)

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(2,5) / {3d} | Ms=1 | 10 | ³F: 7, ³P: 3 |
| | Ms=0 | 25 | ³F: 7, ¹D: 5, ³P: 3, ¹G: 9, ¹S: 1 (NIST order ³F<¹D<³P<¹G<¹S) |
| CAS(4,6) / {3d,4s} | Ms=1 | 120 | ³D: 15, ³F: 35, ³G: 18, ³H: 22, ³P: 15, ⁵D: 5, ⁵F: 7, ⁵P: 3 (TI_EXACT_TERMS, S ≥ 1) |

## Multiplet spectrum (CAS(2,5), our fold, exact ≤1e-9)

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ³F | 7 | −1.468663 | — | −3…3 |
| ¹D | 5 | −1.423747 | +0.0449 | −2…2 |
| ³P | 3 | −1.413969 | +0.0547 | −1…1 |
| ¹G | 9 | −1.398223 | +0.0704 | −4…4 |
| ¹S | 1 | −1.295050 | +0.1736 | 0 |

## Window lattice — exact symmetry vs raw MOs at matched cost

ROHF reference −848.406 Eh; M_S = 1:

| window | CAS | dets | AO, Eh | ⟨L²⟩ | ROHF MO, Eh | ⟨L²⟩ | NO |
|---|---|---|---|---|---|---|---|
| 3d | (4,5) | 50 | −848.054 | 6.000 | — | — | — |
| 3d4s | (4,6) | 120 | −848.395 | 12.000 | −848.405 | 0.000 | −848.399 |
| 3d4d | (4,10) | 1,200 | −848.232 | 6.000 | — | — | — |
| 3d4s4p | (4,9) | 756 | −848.429 | 12.000 | — | — | −848.432 |
| 3d4s4d | (4,11) | 1,815 | −848.409 | 12.000 | −848.414 | 0.010 | — |
| 3d4s4p4d | (4,14) | 5,096 | **−848.438** | 12.000 | −848.420 | 0.201 | **−848.438** |
| 3d4s4p4d5d | (4,19) | 18,411 | −848.439 | 12.000 | — | — | — |
| 3d4s4p4d4f | (4,21) | 27,930 | **−848.446** | 12.000 | — | — | — |

Raw-core caveat: the ground ⟨L²⟩ is UNSHIFTED at percent-class
commutators (4.9e-2 / 2.6e-1) — ⟨L²⟩ alone is not a raw-core detector;
commutators are.

## Treatment ladder (ΔE vs ROHF, mEh)

| rung | ΔE | note |
|---|---|---|
| CASCI, averaged fold | +10.6 | exact |
| CASCI, raw core fold | +10.6 | comm 4.9e-2 — the raw-core detector |
| CASCI, frac reference | +48.0 | |
| CASSCF, free | **−31.3** | escape (window-identity change; terms scramble) |
| CASSCF, window frozen | +10.4 | gated-basin relax honesty (2.3e-8 input defect) |
| CASSCF, radial restriction | −2.1 | single root, chaotic; d² tilts descend (JT class) |
| SA-CASSCF, free | −10.1 | swap basin (totals only; the fold reads −786) |
| SA-CASSCF, frozen, purified core | +11.7 | exact (pur machine-exactness 6.2e-14) |
| SA-CASSCF, radial, purified core | +1.0 | flat-gauge floor (sa_radial_pur 1.4e-4/5.3e-4 above) |

Scan anatomy: raw valence rotations all RISE (+9.8 mEh @ θ = 0.05);
member tilts DESCEND (−163 µEh @ δ = 0.1; the Lz-blind pair pattern) —
the JT drive behind the chaotically broken single roots; the ~11 mEh
radial gain lives on the ORTHOGONIZED substrate's s slot. Ticket-24
amendment: growing the m-complete lattice crosses below the escaped
CASSCF values (5d −4.36, 4f −11.46) — the ceilings are window
statements, not physical ones.

## Open questions

- The ⟨L²⟩ = 0 MO-window phenomenon (same class as the U ⁵S°
  reference; mechanism open, AOC-smearing hypothesis).
- Semicore (3p) treatment — user-deprioritized.
- Spin–orbit splitting of the ³F/¹D multiplets.
