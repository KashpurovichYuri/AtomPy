# Sc-unc (3d¹ 4s², the uncontracted cc-pVDZ split)

<!-- methods:begin -->
| method | route | window | CAS | core | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −759.7374 | 0 | −20673.5179 | +0.0000 | — | — | — | — | — | — | — | forced filling |
| UHF (anchor) | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === |
| CASCI | free | 3d4s / D | (3,6) | Ar·raw | 90 | −759.7342 | +3.2 | −20673.4319 | +0.0860 | — | 5.9990 | 9e-02 | 3e-02 | 4e-02 | 4.6e-02 | 2.5e-01 |  |
| CASCI | free | 3d4s@r1 / D | (3,6) | Ar·raw | 90 | −759.7342 | +3.1 | −20673.4324 | +0.0855 | — | 5.9990 | 9e-02 | 3e-02 | 4e-02 | 4.6e-02 | 2.5e-01 |  |
| CASCI | free | 3d4s@s1 / D | (3,6) | Ar·raw | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d4s@s2 / D | (3,12) | Ar·raw | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (3,6) | Ar | 90 | −759.7342 | +3.2 | −20673.4319 | +0.0860 | — | 5.9989 | 9e-02 | 3e-02 | 4e-02 | 1.4e-02 | 8.4e-02 |  |
| CASCI | free | 3d4s@r1 / D | (3,6) | Ar | 90 | −759.7342 | +3.1 | −20673.4324 | +0.0855 | — | 5.9989 | 9e-02 | 3e-02 | 4e-02 | 1.4e-02 | 8.4e-02 |  |
| CASCI | free | 3d4s@s1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d4s@s2 / D | (3,12) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7470 | −9.6 | −20673.7794 | −0.2616 | — | 0.0060 | 1e-01 | 2e+00 | 3e+00 | 1.0e-01 | 2.6e-01 | m-incomplete |
| CASCI | free | 3d4s / MOC | (3,6) | Ar | 90 | −759.7395 | −2.1 | −20673.5753 | −0.0574 | — | 1.4464 | 1e-01 | 2e+00 | 4e+00 | 1.5e+00 | 6.5e+00 | m-incomplete |
| CASSCF | free | 3d4s / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, virtual-mix | free | 3d4s / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s@r1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s@r1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, virtual-mix | free | 3d4s@r1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s@r1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s@r1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s@s1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s@s1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, virtual-mix | free | 3d4s@s1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s@s1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s@s1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s@s2 / D | (3,12) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s@s2 / D | (3,12) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, virtual-mix | free | 3d4s@s2 / D | (3,12) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s@s2 / D | (3,12) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s@s2 / D | (3,12) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s / MO | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s / MO | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s / MO | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s / MO | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s / MOC | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s / MOC | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s / MOC | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s / MOC | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7327 | +4.7 | −20673.3903 | +0.1276 | — | 6.0000 | 9e-02 | 2e-14 | 2e-14 | 6.3e-15 | 9.3e-14 |  |
| CASCI | free | 3d4s@r1 / D | (3,6) | Ar·pur | 90 | −759.7327 | +4.7 | −20673.3909 | +0.1270 | — | 6.0000 | 9e-02 | 2e-14 | 2e-14 | 1.0e-14 | 9.2e-14 |  |
| CASCI | free | 3d4s@s1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d4s@s2 / D | (3,12) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7455 | −8.1 | −20673.7395 | −0.2216 | — | 0.0060 | 1e-01 | 2e+00 | 3e+00 | 1.0e-01 | 2.6e-01 | m-incomplete |
| CASCI | free | 3d4s / MOC | (3,6) | Ar·pur | 90 | −759.7380 | −0.6 | −20673.5355 | −0.0176 | — | 1.4467 | 1e-01 | 2e+00 | 4e+00 | 1.5e+00 | 6.5e+00 | m-incomplete |
| CASSCF | free | 3d4s / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s@r1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s@r1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s@r1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s@r1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s@s1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s@s1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s@s1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s@s1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s@s2 / D | (3,12) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s@s2 / D | (3,12) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s@s2 / D | (3,12) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s@s2 / D | (3,12) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s / MO | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s / MO | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s / MO | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s / MO | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s / MOC | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s / MOC | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s / MOC | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s / MOC | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d4s / D | (3,6) | Ar·lane | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d4s / D | (3,6) | Ar·lane | 90 | −759.7323 | +5.0 | −20673.3805 | +0.1374 | — | 6.0000 | 9e-02 | 4e-11 | 5e-11 | 6.7e-15 | 4.4e-14 |  |
| CASSCF | free | 3d4s / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | pinned | 3d4s / D | (3,6) | Ar·lane | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === |
| CASCI | free | 3d / D | (3,5) | Ar·raw | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d@r1 / D | (3,5) | Ar·raw | 50 | −759.3375 | +399.9 | −20662.6368 | +10.8810 | — | 11.9994 | 5e-02 | 3e-02 | 4e-02 | 3.6e-02 | 2.1e-01 |  |
| CASCI | free | 3d@s1 / D | (3,5) | Ar·raw | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d@s2 / D | (3,10) | Ar·raw | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d / D | (3,5) | Ar | 50 | −759.3360 | +401.3 | −20662.5968 | +10.9210 | — | 11.9994 | 5e-02 | 3e-02 | 4e-02 | 1.2e-02 | 7.2e-02 |  |
| CASCI | free | 3d@r1 / D | (3,5) | Ar | 50 | −759.3375 | +399.9 | −20662.6368 | +10.8810 | — | 11.9994 | 5e-02 | 3e-02 | 4e-02 | 1.2e-02 | 7.2e-02 |  |
| CASCI | free | 3d@s1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d@s2 / D | (3,10) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d / MO | (3,5) | Ar | 50 | −759.7466 | −9.2 | −20673.7694 | −0.2515 | — | 0.0057 | 1e-01 | 2e+00 | 3e+00 | 8.4e-02 | 2.2e-01 | m-incomplete |
| CASCI | free | 3d / MOC | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, virtual-mix | free | 3d / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d@r1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d@r1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, virtual-mix | free | 3d@r1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d@r1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d@r1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d@s1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d@s1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, virtual-mix | free | 3d@s1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d@s1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d@s1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d@s2 / D | (3,10) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d@s2 / D | (3,10) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, virtual-mix | free | 3d@s2 / D | (3,10) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d@s2 / D | (3,10) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d@s2 / D | (3,10) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d / MO | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d / MO | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d / MO | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d / MO | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d / MOC | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d / MOC | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d / MOC | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d / MOC | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d@r1 / D | (3,5) | Ar·pur | 50 | −759.3374 | +400.0 | −20662.6332 | +10.8847 | — | 12.0000 | 5e-02 | 2e-14 | 2e-14 | 6.6e-15 | 6.5e-14 |  |
| CASCI | free | 3d@s1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d@s2 / D | (3,10) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d / MO | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d / MOC | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d@r1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d@r1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d@r1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d@r1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d@s1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d@s1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d@s1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d@s1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d@s2 / D | (3,10) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d@s2 / D | (3,10) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d@s2 / D | (3,10) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d@s2 / D | (3,10) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d / MO | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d / MO | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d / MO | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d / MO | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d / MOC | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d / MOC | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d / MOC | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d / MOC | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
<!-- methods:end -->
