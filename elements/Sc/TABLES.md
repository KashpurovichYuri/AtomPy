# Sc — scandium (d¹)

cc-pVDZ. Reference **ROHF −759.7362 Eh** (3d¹4s²; M_S = ½; S
half-integer). Section grammar and column keys — README · Results.

## Term table

| CAS / window | Sector | Dets | Term → states |
|---|---|---|---|
| CAS(1,5) / {3d} | Ms=½ | 5 | ²D: 5 |
| CAS(3,6) / {3d,4s} | Ms=½ | 90 | ²D: 20, ⁴F: 14, ²F: 14, ²G: 18, ⁴P: 6, ²P: 6, ²S: 1, ²H: 11 |

## Multiplet energies

Route free · partition averaged fold · source `report.md` (`scripts/run.py Sc`).

| Term | States | Energy (Eh) | Gap to ground | M_L set |
|---|---|---|---|---|
| ²D | 5 | −0.342549 | — | −2…2 |

## CASCI — window survey

CI level · ΔE anchored on 3d4s pinned AO (the certified row) · source `cells.json` (as-run) · pins `test_symm_tm.py::test_cas_ladder_*`.

<!-- survey:begin -->
| window | CAS | fold | start | route | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3d | (3,5) | Ar | AO | free | 50 | −759.4359 | +291.6 | −20665.3136 | +7.9335 | 11.9999 | 7.0e-04 | 3.1e-02 |
| 3d | (3,5) | Ar | MO | free | 50 | −759.7467 | −19.3 | −20673.7714 | −0.5242 | 0.0065 | 1.1e-01 | 2.2e-01 |
| 3d4s | (3,6) | Ar | AO | free | 90 | −759.7296 | −2.1 | −20673.3051 | −0.0579 | 5.9998 | 8.9e-04 | 4.0e-02 |
| 3d4s | (3,6) | Ar | MO | free | 90 | −759.7470 | −19.6 | −20673.7795 | −0.5324 | 0.0062 | 1.5e-01 | 7.2e-01 |
| 3d4s | (3,6) | Ar·raw | AO | free | 90 | −759.7296 | −2.1 | −20673.3051 | −0.0579 | 6.0003 | 2.8e-02 | 1.3e-01 |
| 3d4s | (3,6) | Ar·pur | AO | free | 90 | −759.7289 | −1.5 | −20673.2867 | −0.0396 | 6.0000 | 1.6e-13 | 5.6e-13 |
| 3d4s | (3,6) | Ar·frac | AO | free | 90 | −759.6995 | +27.9 | −20672.4868 | +0.7604 | 6.0000 | 1.3e-13 | 5.4e-13 |
| 3d4s | (3,6) | Ar | AO | pinned | 90 | −759.7274 | 0 | −20673.2471 | +0.0000 | 6.0000 | 3.4e-14 | 3.7e-13 |
| 3d4s | (3,6) | Ar | AO | free | 90 | −759.7274 | +0.0 | −20673.2471 | +0.0000 | 6.0000 | 6.5e-14 | 2.4e-13 |
| 3d4d | (3,10) | Ar | AO | free | 450 | −759.5553 | +172.1 | −20668.5639 | +4.6833 | 11.9999 | 2.2e-03 | 1.2e-01 |
| 3d4d | (3,10) | Ar | MO | free | 450 | −759.7482 | −20.8 | −20673.8133 | −0.5662 | 1.5780 | 4.7e+00 | 1.3e+01 |
| 3d4s4p | (3,9) | Ar | AO | free | 324 | −759.7611 | −33.7 | −20674.1629 | −0.9157 | 5.9998 | 1.4e-03 | 5.8e-02 |
| 3d4s4p | (3,9) | Ar | MO | free | 324 | −759.7481 | −20.7 | −20673.8103 | −0.5632 | 0.6272 | 3.7e-01 | 8.2e+00 |
| 3d4s4d | (3,11) | Ar | AO | free | 605 | −759.7402 | −12.8 | −20673.5943 | −0.3472 | 5.9998 | 2.5e-03 | 1.3e-01 |
| 3d4s4d | (3,11) | Ar | MO | free | 605 | −759.7484 | −21.0 | −20673.8183 | −0.5711 | 1.6033 | 5.8e+00 | 2.6e+01 |
| 3d4s4p4d | (3,14) | Ar | AO | free | 1,274 | −759.7658 | −38.4 | −20674.2909 | −1.0437 | 5.9998 | 3.3e-03 | 1.6e-01 |
| 3d4s4p4d | (3,14) | Ar | MO | free | 1,274 | −759.7510 | −23.5 | −20673.8875 | −0.6404 | 2.8678 | 6.9e+00 | 2.8e+01 |
| 3d4s4p4d5d | (3,19) | Ar | AO | free | 3,249 | −759.7660 | −38.5 | −20674.2961 | −1.0489 | 5.9998 | 4.2e-03 | 2.1e-01 |
| 3d4s4p4d5d | (3,19) | Ar | MO | free | 3,249 | −759.7696 | −42.2 | −20674.3942 | −1.1471 | 5.5905 | 1.1e+01 | 4.9e+01 |
| 3d4s4p4d4f | (3,21) | Ar | AO | free | 4,410 | −759.7697 | −42.3 | −20674.3985 | −1.1513 | 5.9998 | 1.4e-01 | 8.0e-01 |
| 3d4s4p4d4f | (3,21) | Ar | MO | free | 4,410 | −759.7720 | −44.6 | −20674.4604 | −1.2133 | 5.6353 | 1.2e+01 | 7.2e+01 |
| 3d3p | (7,8) | semi | AO | pinned | 3,920 | −759.7871 | −59.7 | −20674.8712 | −1.6241 | 6.0000 | 7.1e-13 | 4.9e-12 |
| 3d3p | (7,8) | semi | AO | free | 3,920 | −759.7871 | −59.7 | −20674.8712 | −1.6241 | 6.0000 | 6.4e-13 | 5.2e-12 |
| 3d3p4d | (7,13) | semi | AO | pinned | 204,490 | −759.8409 | −113.4 | −20676.3338 | −3.0867 | 6.0000 | 3.6e-14ψ | 9.6e-14ψ |
| 3d3p4d | (7,13) | semi | AO | free | 204,490 | −759.8409 | −113.4 | −20676.3338 | −3.0867 | 6.0000 | 1.5e-14ψ | 3.8e-14ψ |
| 3d3p4s | (9,9) | {1s,2s,2p,3s} | AO | pinned | 15,876 | −759.7899 | −62.4 | −20674.9462 | −1.6990 | 6.0000 | 7.8e-15ψ | 6.1e-14ψ |
| 3d3p4s | (9,9) | {1s,2s,2p,3s} | AO | free | 15,876 | −759.7899 | −62.4 | −20674.9462 | −1.6990 | 6.0000 | 2.7e-15ψ | 3.8e-14ψ |
| 3s3p3d | (9,9) | {1s,2s,2p,4s} | AO | pinned | 15,876 | −759.7900 | −62.5 | −20674.9492 | −1.7021 | 6.0000 | 4.3e-14ψ | 1.2e-13ψ |
| 3s3p3d | (9,9) | {1s,2s,2p,4s} | AO | free | 15,876 | −759.7900 | −62.5 | −20674.9492 | −1.7021 | 6.0000 | 1.7e-14ψ | 7.2e-14ψ |
<!-- survey:end -->

The semi-core window pays −59.7 mEh over the declared window at
frozen shells — the 3p⁶ correlation — with ²D hosted exactly
(components degenerate, every meter exact).  The extension windows place
the pair-correlating 4d shell and the 4s/3s pairs against it.

## Method ladder, CAS(3,6)/{3d,4s}

Optimizer rows across CAS windows · ΔE vs the ROHF anchor (-759.7363 Eh) · routes in rows · source `cells.json` (as-run).

| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| — | ROHF (anchor) | — | — | — | — | −759.7363 | 0 | −20673.4871 | +0.0000 | — | — | — | — | — |
| free | CASSCF, free | 3d4s | (3,6) | Ar | 90 | −759.7680 | −31.7 | −20674.3496 | −0.8625 | 1e-8/54 | 0.0469 | — | 5.7e-01 | 1.2e+00 |
| free | CASSCF, free, MO | 3d4s | (3,6) | Ar | 90 | −759.7675 | −31.3 | −20674.3377 | −0.8506 | 1e-8/24 | 0.0533 | — | 3.8e-01 | 7.6e-01 |
| free | CASSCF, free, NO | 3d4sno2 | (3,2) | Ar | 2 | −759.7366 | −0.3 | −20673.4952 | −0.0081 | 0/152 | 0.0000 | — | 3.2e-06 | 3.6e-21 |
| free | CASSCF, virtual-mix | 3d4s | (3,6) | Ar | 90 | −759.7676 | −31.3 | −20674.3397 | −0.8525 | 1e-8/73 | 0.0465 | — | 5.9e-01 | 1.2e+00 |
| free | CASSCF, free, pur (+no) | 3d4s | (3,6) | Ar·pur | 90 | −759.7678 | −31.5 | −20674.3447 | −0.8575 | 1e-8/95 | 0.0526 | — | 3.5e-01 | 7.2e-01 |
| free | CASSCF, free, pur, MO | 3d4s | (3,6) | Ar·pur | 90 | −759.7675 | −31.3 | −20674.3377 | −0.8506 | 1e-8/24 | 0.0533 | — | 3.8e-01 | 7.6e-01 |
| free | CASSCF, radial | 3d4s | (3,6) | Ar | 90 | −759.7431 | −6.9 | −20673.6741 | −0.1870 | 1e-12/15 | 3.5331 | — | 1.6e+00 | 4.0e+00 |
| free | CASSCF, radial, pur | 3d4s | (3,6) | Ar·pur | 90 | −759.7425 | −6.2 | −20673.6558 | −0.1686 | 1e-10/15 | 3.5384 | — | 1.6e+00 | 4.0e+00 |
| free | CASSCF, frozen | 3d4s | (3,6) | Ar | 90 | −759.7297 | +6.6 | −20673.3075 | +0.1796 | 1e-11/5 | 5.9998 | — | 8.9e-04 | 4.0e-02 |
| free | CASSCF, frozen, MO | 3d4s | (3,6) | Ar | 90 | −759.7470 | −10.7 | −20673.7795 | −0.2924 | 1e-11/4 | 0.0062 | — | 1.5e-01 | 7.2e-01 |
| free | CASSCF, frozen, NO | 3d4sno2 | (3,2) | Ar | 2 | −759.7266 | +9.6 | −20673.2251 | +0.2620 | 0/6 | 0.0000 | — | 5.9e-18 | 1.1e-32 |
| free | CASSCF, frozen, pur (+no) | 3d4s | (3,6) | Ar·pur | 90 | −759.7295 | +6.7 | −20673.3042 | +0.1830 | 0/6 | 6.0000 | — | 3.4e-14 | 1.7e-13 |
| free | CASSCF, frozen, pur, MO | 3d4s | (3,6) | Ar·pur | 90 | −759.7470 | −10.7 | −20673.7786 | −0.2915 | 1e-10/5 | 0.0062 | — | 1.5e-01 | 7.2e-01 |
| free | SA-CASSCF, free | 3d4s | (3,6) | Ar | 90 | −759.7495 | −13.3 | −20673.8478 | −0.3607 | 1e-12/129 | 6.0020 | — | 7.9e-01 | 3.8e+00 |
| free | SA-CASSCF, free, MO | 3d4s | (3,6) | Ar | 90 | −759.7626 | −26.3 | −20674.2032 | −0.7161 | 1e-11/14 | 0.0309 | — | 3.1e-01 | 1.2e+00 |
| free | SA-CASSCF, free, pur | 3d4s | (3,6) | Ar·pur | 90 | −759.7397 | −3.5 | −20673.5819 | −0.0948 | — | 6.0000 | — | 4.6e-05 | 3.0e-04 |
| free | SA-CASSCF, radial | 3d4s | (3,6) | Ar | 90 | −759.7370 | −0.8 | −20673.5085 | −0.0214 | — | 5.9997 | — | 1.2e-03 | 5.8e-02 |
| free | SA-CASSCF, radial, pur | 3d4s | (3,6) | Ar·pur | 90 | −759.7362 | +0.0 | −20673.4869 | +0.0002 | — | 6.0000 | — | 7.0e-08 | 2.0e-05 |
| free | SA-CASSCF, frozen, pur | 3d4s | (3,6) | Ar·pur | 90 | −759.7290 | +7.3 | −20673.2891 | +0.1980 | — | 6.0000 | — | 6.5e-14 | 2.1e-13 |


### CAS(3,6)/{3d,4s} · Ar

Optimizer rows across CAS windows · ΔE vs the ROHF anchor (-759.7363 Eh) · routes in rows · source `cells.json` (as-run).

| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | CASSCF | 3d4s | (3,6) | Ar | 90 | −759.7436 | −7.4 | −20673.6877 | −0.2005 | 1e-9/18 | 3.5457 | — | 1.5e+00 | 3.7e+00 |
| free | SA-CASSCF | 3d4s | (3,6) | Ar | 90 | −759.7495 | −13.3 | −20673.8479 | −0.3607 | 1e-9/75 | 6.0020 | — | 2.3e-01 | 3.8e+00 |
| pinned | SA-CASSCF | 3d4s | (3,6) | Ar | 90 | −759.7397 | −3.5 | −20673.5819 | −0.0948 | 1e-13/32 | 6.0000 | 2e-06 | 1.8e-05 | 1.0e-04 |


### CAS(7,8)/{3d,3p} · semi

Optimizer rows across CAS windows · ΔE vs the ROHF anchor (-759.7363 Eh) · routes in rows · source `cells.json` (as-run).

| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | CASSCF | 3d3p | (7,8) | semi | 3,920 | −759.8454 | −109.1 | −20676.4572 | −2.9701 | 1e-9/13 | 5.0700 | — | 6.2e+00 | 6.0e+01 |
| free | SA-CASSCF | 3d3p | (7,8) | semi | 3,920 | −759.8218 | −85.5 | −20675.8147 | −2.3276 | 1e-11/12 | 6.0000 | — | 5.5e-05 | 1.7e-03 |
| pinned | SA-CASSCF | 3d3p | (7,8) | semi | 3,920 | −759.8218 | −85.5 | −20675.8147 | −2.3276 | 1e-9/10 | 6.0000 | 7e-07 | 2.2e-05 | 3.8e-04 |


### CAS(7,13)/{3d,3p,4d} · semi

Optimizer rows across CAS windows · ΔE vs the ROHF anchor (-759.7363 Eh) · routes in rows · source `cells.json` (as-run).

| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | CASSCF | 3d3p4d | (7,13) | semi | 204,490 | −759.8668 | −130.6 | −20677.0398 | −3.5527 | 1e-6/200 | 4.8569 | — | 2.5e-01ψ | 1.3e+00ψ |
| free | SA-CASSCF | 3d3p4d | (7,13) | semi | 204,490 | −759.8523 | −116.0 | −20676.6438 | −3.1567 | 1e-9/17 | 6.0000 | — | 3.7e-05ψ | 2.5e-04ψ |
| pinned | SA-CASSCF | 3d3p4d | (7,13) | semi | 204,490 | −759.8523 | −116.0 | −20676.6438 | −3.1567 | 1e-9/16 | 6.0000 | 5e-08 | 3.4e-07ψ | 2.2e-06ψ |


### CAS(9,9)/{3d,3p,4s} · {1s,2s,2p,3s}

Optimizer rows across CAS windows · ΔE vs the ROHF anchor (-759.7363 Eh) · routes in rows · source `cells.json` (as-run).

| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | CASSCF | 3d3p4s | (9,9) | {1s,2s,2p,3s} | 15,876 | −759.8511 | −114.8 | −20676.6116 | −3.1245 | — | 5.0601 | — | 6.3e-05ψ | 2.8e-01ψ |
| free | SA-CASSCF | 3d3p4s | (9,9) | {1s,2s,2p,3s} | 15,876 | −759.8266 | −90.4 | −20675.9462 | −2.4591 | — | 6.0000 | — | 7.2e-07ψ | 2.4e-06ψ |
| pinned | SA-CASSCF | 3d3p4s | (9,9) | {1s,2s,2p,3s} | 15,876 | −759.8266 | −90.4 | −20675.9462 | −2.4591 | — | 6.0000 | 8e-07 | 1.1e-06ψ | 1.2e-06ψ |


### CAS(9,9)/{3s,3p,3d} · {1s,2s,2p,4s}

Optimizer rows across CAS windows · ΔE vs the ROHF anchor (-759.7363 Eh) · routes in rows · source `cells.json` (as-run).

| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | CASSCF | 3s3p3d | (9,9) | {1s,2s,2p,4s} | 15,876 | −759.8509 | −114.7 | −20676.6077 | −3.1206 | 1e-11/24 | 5.0674 | — | 2.1e-01ψ | 3.2e-01ψ |
| free | SA-CASSCF | 3s3p3d | (9,9) | {1s,2s,2p,4s} | 15,876 | −759.8266 | −90.4 | −20675.9462 | −2.4591 | 1e-8/18 | 6.0000 | — | 4.9e-07ψ | 2.8e-06ψ |
| pinned | SA-CASSCF | 3s3p3d | (9,9) | {1s,2s,2p,4s} | 15,876 | −759.8266 | −90.4 | −20675.9462 | −2.4591 | 0/55 | 6.0000 | 2e-07 | 2.4e-07ψ | 2.3e-06ψ |

## Method ladder, CAS(3,6)/{3d,4s}

Optimizer rows across CAS windows · ΔE vs the ROHF anchor (-759.7363 Eh) · routes in rows · source `cells.json` (as-run).

| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| — | ROHF (anchor) | — | — | — | — | −759.7363 | 0 | −20673.4871 | +0.0000 | — | — | — | — | — |
| free | CASSCF, free | 3d4s | (3,6) | averaged | 90 | −759.7680 | −31.7 | −20674.3496 | −0.8625 | 1e-8/54 | 0.0469 | — | 5.7e-01 | 1.2e+00 |
| free | CASSCF, free, MO | 3d4s | (3,6) | averaged | 90 | −759.7675 | −31.3 | −20674.3377 | −0.8506 | 1e-8/24 | 0.0533 | — | 3.8e-01 | 7.6e-01 |
| free | CASSCF, free, NO | 3d4sno2 | (3,6) | averaged | 2 | −759.7366 | −0.3 | −20673.4952 | −0.0081 | 0/152 | 0.0000 | — | 3.2e-06 | 3.6e-21 |
| free | CASSCF, virtual−mix | 3d4s | (3,6) | averaged | 90 | −759.7676 | −31.3 | −20674.3397 | −0.8525 | 1e-8/73 | 0.0465 | — | 5.9e-01 | 1.2e+00 |
| free | CASSCF, free, pur (+no) | 3d4s | (3,6) | averaged | 90 | −759.7678 | −31.5 | −20674.3447 | −0.8575 | 1e-8/95 | 0.0526 | — | 3.5e-01 | 7.2e-01 |
| free | CASSCF, free, pur, MO | 3d4s | (3,6) | averaged | 90 | −759.7675 | −31.3 | −20674.3377 | −0.8506 | 1e-8/24 | 0.0533 | — | 3.8e-01 | 7.6e-01 |
| free | CASSCF, radial | 3d4s | (3,6) | averaged | 90 | −759.7431 | −6.9 | −20673.6741 | −0.1870 | 1e-12/15 | 3.5331 | — | 1.6e+00 | 4.0e+00 |
| free | CASSCF, radial, pur | 3d4s | (3,6) | averaged | 90 | −759.7425 | −6.2 | −20673.6558 | −0.1686 | 1e-10/15 | 3.5384 | — | 1.6e+00 | 4.0e+00 |
| free | CASSCF, frozen | 3d4s | (3,6) | averaged | 90 | −759.7297 | +6.6 | −20673.3075 | +0.1796 | 1e-11/5 | 5.9998 | — | 8.9e-04 | 4.0e-02 |
| free | CASSCF, frozen, MO | 3d4s | (3,6) | averaged | 90 | −759.7470 | −10.7 | −20673.7795 | −0.2924 | 1e-11/4 | 0.0062 | — | 1.5e-01 | 7.2e-01 |
| free | CASSCF, frozen, NO | 3d4sno2 | (3,6) | averaged | 2 | −759.7266 | +9.6 | −20673.2251 | +0.2620 | 0/6 | 0.0000 | — | 5.9e-18 | 1.1e-32 |
| free | CASSCF, frozen, pur (+no) | 3d4s | (3,6) | averaged | 90 | −759.7295 | +6.7 | −20673.3042 | +0.1830 | 0/6 | 6.0000 | — | 3.4e-14 | 1.7e-13 |
| free | CASSCF, frozen, pur, MO | 3d4s | (3,6) | averaged | 90 | −759.7470 | −10.7 | −20673.7786 | −0.2915 | 1e-10/5 | 0.0062 | — | 1.5e-01 | 7.2e-01 |
| free | SA−CASSCF, free | 3d4s | (3,6) | averaged | 90 | −759.7495 | −13.3 | −20673.8478 | −0.3607 | 1e-12/129 | 6.0020 | — | 7.9e-01 | 3.8e+00 |
| free | SA−CASSCF, free, MO | 3d4s | (3,6) | averaged | 90 | −759.7626 | −26.3 | −20674.2032 | −0.7161 | 1e-11/14 | 0.0309 | — | 3.1e-01 | 1.2e+00 |
| free | SA−CASSCF, free, pur | 3d4s | (3,6) | averaged | 90 | −759.7397 | −3.5 | −20673.5819 | −0.0948 | — | 6.0000 | — | 4.6e-05 | 3.0e-04 |
| free | SA−CASSCF, radial | 3d4s | (3,6) | averaged | 90 | −759.7370 | −0.8 | −20673.5085 | −0.0214 | — | 5.9997 | — | 1.2e-03 | 5.8e-02 |
| free | SA−CASSCF, radial, pur | 3d4s | (3,6) | averaged | 90 | −759.7362 | +0.0 | −20673.4869 | +0.0002 | — | 6.0000 | — | 7.0e-08 | 2.0e-05 |
| free | SA−CASSCF, frozen, pur | 3d4s | (3,6) | averaged | 90 | −759.7290 | +7.3 | −20673.2891 | +0.1980 | — | 6.0000 | — | 6.5e-14 | 2.1e-13 |


### CAS(3,6)/{3d,4s}

Optimizer rows across CAS windows · ΔE vs the ROHF anchor (-759.7363 Eh) · routes in rows · source `cells.json` (as-run).

| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | CASSCF | 3d4s | (3,6) | Ar | 90 | −759.7436 | −7.4 | −20673.6877 | −0.2005 | 1e-9/18 | 3.5457 | — | 1.5e+00 | 3.7e+00 |
| free | SA−CASSCF | 3d4s | (3,6) | Ar | 90 | −759.7495 | −13.3 | −20673.8479 | −0.3607 | 1e-9/75 | 6.0020 | — | 2.3e-01 | 3.8e+00 |
| pinned | SA−CASSCF | 3d4s | (3,6) | Ar | 90 | −759.7397 | −3.5 | −20673.5819 | −0.0948 | 1e-13/32 | 6.0000 | 2e-06 | 1.8e-05 | 1.0e-04 |


### CAS(7,8)/{3d,3p} · semi

Optimizer rows across CAS windows · ΔE vs the ROHF anchor (-759.7363 Eh) · routes in rows · source `cells.json` (as-run).

| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | CASSCF | 3d3p | (7,8) | semi | 3,920 | −759.8454 | −109.1 | −20676.4572 | −2.9701 | 1e-9/13 | 5.0700 | — | 6.2e+00 | 6.0e+01 |
| free | SA−CASSCF | 3d3p | (7,8) | semi | 3,920 | −759.8218 | −85.5 | −20675.8147 | −2.3276 | 1e-11/12 | 6.0000 | — | 5.5e-05 | 1.7e-03 |
| pinned | SA−CASSCF | 3d3p | (7,8) | semi | 3,920 | −759.8218 | −85.5 | −20675.8147 | −2.3276 | 1e-9/10 | 6.0000 | 7e-07 | 2.2e-05 | 3.8e-04 |


### CAS(7,13)/{3d,3p,4d} · semi

Optimizer rows across CAS windows · ΔE vs the ROHF anchor (-759.7363 Eh) · routes in rows · source `cells.json` (as-run).

| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | CASSCF | 3d3p4d | (7,13) | semi | 204,490 | −759.8668 | −130.6 | −20677.0398 | −3.5527 | 1e-6/200 | 4.8569 | — | 2.5e-01ψ | 1.3e+00ψ |
| free | SA−CASSCF | 3d3p4d | (7,13) | semi | 204,490 | −759.8523 | −116.0 | −20676.6438 | −3.1567 | 1e-9/17 | 6.0000 | — | 3.7e-05ψ | 2.5e-04ψ |
| pinned | SA−CASSCF | 3d3p4d | (7,13) | semi | 204,490 | −759.8523 | −116.0 | −20676.6438 | −3.1567 | 1e-9/16 | 6.0000 | 5e-08 | 3.4e-07ψ | 2.2e-06ψ |


### CAS(9,9)/{3d,3p,4s} · {1s,2s,2p,3s}

Optimizer rows across CAS windows · ΔE vs the ROHF anchor (-759.7363 Eh) · routes in rows · source `cells.json` (as-run).

| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | CASSCF | 3d3p4s | (9,9) | {1s,2s,2p,3s} | 15,876 | −759.8511 | −114.8 | −20676.6116 | −3.1245 | — | 5.0601 | — | 6.3e-05ψ | 2.8e-01ψ |
| free | SA−CASSCF | 3d3p4s | (9,9) | {1s,2s,2p,3s} | 15,876 | −759.8266 | −90.4 | −20675.9462 | −2.4591 | — | 6.0000 | — | 7.2e-07ψ | 2.4e-06ψ |
| pinned | SA−CASSCF | 3d3p4s | (9,9) | {1s,2s,2p,3s} | 15,876 | −759.8266 | −90.4 | −20675.9462 | −2.4591 | — | 6.0000 | 8e-07 | 1.1e-06ψ | 1.2e-06ψ |


### CAS(9,9)/{3s,3p,3d} · {1s,2s,2p,4s}

Optimizer rows across CAS windows · ΔE vs the ROHF anchor (-759.7363 Eh) · routes in rows · source `cells.json` (as-run).

| route | method | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[H,Lz]‖ | ‖[H,L²]‖ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | CASSCF | 3s3p3d | (9,9) | {1s,2s,2p,4s} | 15,876 | −759.8509 | −114.7 | −20676.6077 | −3.1206 | 1e-11/24 | 5.0674 | — | 2.1e-01ψ | 3.2e-01ψ |
| free | SA−CASSCF | 3s3p3d | (9,9) | {1s,2s,2p,4s} | 15,876 | −759.8266 | −90.4 | −20675.9462 | −2.4591 | 1e-8/18 | 6.0000 | — | 4.9e-07ψ | 2.8e-06ψ |
| pinned | SA−CASSCF | 3s3p3d | (9,9) | {1s,2s,2p,4s} | 15,876 | −759.8266 | −90.4 | −20675.9462 | −2.4591 | 0/55 | 6.0000 | 2e-07 | 2.4e-07ψ | 2.3e-06ψ |

## Studies (frozen history; ticket 28)

<!-- study:begin --> (frozen history: the generator (atomterm.study) is retired, commit ad28d3a -- the block stays as the committed record; the lane/lattice/methods blocks are the living generated tables)

| cell | status | E, Eh | ⟨L²⟩ | ‖[P,Lz]‖ | ‖[P,Lp]‖ |
|---|---|---|---|---|---|
| study_sc_no_no2_casci | ok | -759.7306 | 6.0615 | 1.9e-14 | 3.0e-14 |
| study_sc_no_no2_relax | ok | -759.7312 | 6.0602 | 1.9e-14 | 3.0e-14 |
| study_sc_no_no2_sa-casscf | ok | -759.7495 | 8.9253 | 1.9e-14 | 3.0e-14 |
| study_sc_no_no3_casci | ok | -759.7622 | 6.5059 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no3_relax | ok | -759.7627 | 6.5047 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no3_sa-casscf | ok | -759.7657 | 6.3497 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no4_casci | ok | -759.7651 | 6.6032 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no4_relax | ok | -759.7655 | 6.6014 | 2.2e-13 | 3.5e-13 |
| study_sc_no_no4_sa-casscf | ok | -759.7710 | 6.4999 | 2.2e-13 | 3.5e-13 |

<!-- study:end -->
