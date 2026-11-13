# Sc-unc (3d¹ 4s², the uncontracted cc-pVDZ split)

<!-- methods:begin -->
| method | route | window | CAS | core | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −759.7374 | 0 | −20673.5179 | +0.0000 | — | — | — | — | — | — | — | forced filling |
| UHF (anchor) | — | — | — | — | — | −759.7394 | −2.0 | −20673.5733 | −0.0554 | — | — | — | — | — | — | — | ⟨S²⟩=0.7534, σ_min=0.9989 |
| === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === |
| CASCI | free | 3d4s / D | (3,6) | Ar·raw | 90 | −759.7342 | +3.2 | −20673.4319 | +0.0860 | — | 5.9990 | 9e-02 | 3e-02 | 4e-02 | 4.6e-02 | 2.5e-01 |  |
| CASCI | free | 3d4s@r1 / D | (3,6) | Ar·raw | 90 | −759.7342 | +3.1 | −20673.4324 | +0.0855 | — | 5.9990 | 9e-02 | 3e-02 | 4e-02 | 4.6e-02 | 2.5e-01 | bands 3d4s |
| CASCI | free | 3d4s@s1 / D | (3,6) | Ar·raw | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (3,6) | Ar | 90 | −759.7342 | +3.2 | −20673.4319 | +0.0860 | — | 5.9989 | 9e-02 | 3e-02 | 4e-02 | 1.4e-02 | 8.4e-02 |  |
| CASCI | free | 3d4s@r1 / D | (3,6) | Ar | 90 | −759.7342 | +3.1 | −20673.4324 | +0.0855 | — | 5.9989 | 9e-02 | 3e-02 | 4e-02 | 1.4e-02 | 8.4e-02 | bands 3d4s |
| CASCI | free | 3d4s@s1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7470 | −9.6 | −20673.7794 | −0.2616 | — | 0.0060 | 1e-01 | 2e+00 | 3e+00 | 1.0e-01 | 2.6e-01 | m-incomplete |
| CASCI | free | 3d4s / MOC | (3,6) | Ar | 90 | −759.7395 | −2.1 | −20673.5753 | −0.0574 | — | 1.4464 | 1e-01 | 2e+00 | 4e+00 | 1.5e+00 | 6.5e+00 | m-incomplete |
| CASSCF | free | 3d4s / D | (3,6) | Ar | 90 | −759.7593 | −21.9 | −20674.1143 | −0.5965 | 1e-5/195 | 1.7935 | — | 5e+00 | 5e+00 | 3.5e+00 | 1.2e+01 | swap; ; local minimum |
| CASSCF, frozen | free | 3d4s / D | (3,6) | Ar | 90 | −759.7343 | +3.1 | −20673.4331 | +0.0848 | 0/104 | 5.9989 | — | 3e-02 | 4e-02 | 1.4e-02 | 8.4e-02 |  |
| CASSCF, frozen core | free | 3d4s / D | (3,6) | Ar | 90 | −759.7343 | +3.0 | −20673.4353 | +0.0826 | 1e-5/157 | 5.9987 | — | 4e-02 | 6e-02 | 3.2e-02 | 3.5e-01 | not converged |
| SA-CASSCF | free | 3d4s / D | (3,6) | Ar | 90 | −759.7497 | −12.3 | −20673.8517 | −0.3338 | 1e-6/184 | 6.0035 | 8e-03 | 1e+00 | 2e+00 | 7.3e-01 | 4.0e+00 | swap |
| SA-CASSCF, frozen | free | 3d4s / D | (3,6) | Ar | 90 | −759.7332 | +4.2 | −20673.4043 | +0.1136 | 0/112 | 5.9989 | 2e-03 | 3e-02 | 4e-02 | 1.4e-02 | 8.4e-02 |  |
| SA-CASSCF, frozen core | free | 3d4s / D | (3,6) | Ar | 90 | −759.7280 | +9.4 | −20673.2621 | +0.2558 | 1e-4/116 | 6.0733 | 9e-02 | 5e-02 | 8e-02 | 7.1e+02 | 3.7e+03 | not converged |
| CASSCF | free | 3d4s@r1 / D | (3,6) | Ar | 90 | −759.7604 | −23.1 | −20674.1452 | −0.6274 | 1e-5/196 | 1.9818 | — | 3e+00 | 5e+00 | 4.8e+00 | 2.3e+01 | swap; ; local minimum; bands 3d4s |
| CASSCF, frozen | free | 3d4s@r1 / D | (3,6) | Ar | 90 | −759.7343 | +3.1 | −20673.4336 | +0.0842 | 0/104 | 5.9989 | — | 3e-02 | 4e-02 | 1.4e-02 | 8.4e-02 | bands 3d4s |
| CASSCF, frozen core | free | 3d4s@r1 / D | (3,6) | Ar | 90 | −759.7338 | +3.6 | −20673.4196 | +0.0982 | 1e-4/125 | 6.0081 | — | 3e-02 | 5e-02 | 2.5e+01 | 2.3e+02 | not converged; ; bands 3d4s |
| SA-CASSCF | free | 3d4s@r1 / D | (3,6) | Ar | 90 | −759.7495 | −12.1 | −20673.8477 | −0.3298 | 1e-5/187 | 6.0044 | 8e-03 | 1e+00 | 1e+00 | 8.6e-01 | 4.0e+00 | swap; ; bands 3d4s |
| SA-CASSCF, frozen | free | 3d4s@r1 / D | (3,6) | Ar | 90 | −759.7332 | +4.2 | −20673.4048 | +0.1131 | 0/114 | 5.9989 | 2e-03 | 3e-02 | 4e-02 | 1.4e-02 | 8.4e-02 | bands 3d4s |
| SA-CASSCF, frozen core | free | 3d4s@r1 / D | (3,6) | Ar | 90 | −759.7203 | +17.1 | −20673.0536 | +0.4643 | 1e-3/129 | 6.0295 | 1e-02 | 3e-02 | 4e-02 | 1.3e-01 | 6.8e-01 | not converged; ; bands 3d4s |
| CASSCF | free | 3d4s@s1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s@s1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen core | free | 3d4s@s1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s@s1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s@s1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen core | free | 3d4s@s1 / D | (3,6) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7687 | −31.3 | −20674.3699 | −0.8520 | 1e-8/18 | 0.0541 | — | 2e+00 | 3e+00 | 2.8e-01 | 7.2e-01 | m-incomplete |
| CASSCF, frozen | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7470 | −9.6 | −20673.7795 | −0.2616 | 1e-12/5 | 0.0060 | — | 2e+00 | 3e+00 | 1.0e-01 | 2.6e-01 | m-incomplete |
| CASSCF, frozen core | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7638 | −26.4 | −20674.2362 | −0.7184 | 1e-13/103 | 0.0352 | — | 2e+00 | 3e+00 | 2.2e-01 | 5.9e-01 | not converged; m-incomplete |
| SA-CASSCF | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7636 | −26.3 | −20674.2324 | −0.7146 | 1e-7/112 | 0.0309 | 7e-02 | 3e+00 | 4e+00 | 5.0e-01 | 1.2e+00 | m-incomplete |
| SA-CASSCF, frozen | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7469 | −9.5 | −20673.7771 | −0.2593 | 1e-11/6 | 0.0061 | 8e-02 | 2e+00 | 3e+00 | 1.0e-01 | 2.6e-01 | m-incomplete |
| SA-CASSCF, frozen core | free | 3d4s / MO | (3,6) | Ar | 90 | −759.5728 | +164.6 | −20669.0389 | +4.4789 | 1e-12/100 | 1.9999 | 2e-01 | 2e+00 | 3e+00 | 6.2e+00 | 1.3e+01 | not converged; m-incomplete |
| CASSCF | free | 3d4s / MOC | (3,6) | Ar | 90 | −759.7446 | −7.2 | −20673.7135 | −0.1956 | 1e-8/20 | 3.4898 | — | 2e+00 | 3e+00 | 1.1e+00 | 4.1e+00 | m-incomplete; local minimum |
| CASSCF, frozen | free | 3d4s / MOC | (3,6) | Ar | 90 | −759.7395 | −2.1 | −20673.5753 | −0.0574 | 1e-8/2 | 1.4464 | — | 2e+00 | 4e+00 | 1.5e+00 | 6.5e+00 | m-incomplete |
| CASSCF, frozen core | free | 3d4s / MOC | (3,6) | Ar | 90 | −759.7066 | +30.8 | −20672.6800 | +0.8379 | 0/103 | 3.4454 | — | 2e+00 | 3e+00 | 5.1e+00 | 3.2e+01 | not converged; m-incomplete |
| SA-CASSCF | free | 3d4s / MOC | (3,6) | Ar | 90 | −759.7376 | −0.2 | −20673.5230 | −0.0051 | 0/22 | 5.8000 | 3e-02 | 7e-01 | 1e+00 | 7.2e-01 | 4.6e+00 | m-incomplete |
| SA-CASSCF, frozen | free | 3d4s / MOC | (3,6) | Ar | 90 | −759.7395 | −2.1 | −20673.5743 | −0.0564 | 0/6 | 1.4464 | 1e-01 | 2e+00 | 4e+00 | 1.6e+00 | 6.5e+00 | m-incomplete |
| SA-CASSCF, frozen core | free | 3d4s / MOC | (3,6) | Ar | 90 | −759.5667 | +170.7 | −20668.8728 | +4.6451 | 1e-11/103 | 10.8393 | 8e-03 | 8e-01 | 1e+00 | 3.1e+01 | 1.9e+02 | not converged; m-incomplete |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7327 | +4.7 | −20673.3903 | +0.1276 | — | 6.0000 | 9e-02 | 2e-14 | 2e-14 | 6.3e-15 | 9.3e-14 |  |
| CASCI | free | 3d4s@r1 / D | (3,6) | Ar·pur | 90 | −759.7327 | +4.7 | −20673.3909 | +0.1270 | — | 6.0000 | 9e-02 | 2e-14 | 2e-14 | 1.0e-14 | 9.2e-14 | bands 3d4s |
| CASCI | free | 3d4s@s1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7455 | −8.1 | −20673.7395 | −0.2216 | — | 0.0060 | 1e-01 | 2e+00 | 3e+00 | 1.0e-01 | 2.6e-01 | m-incomplete |
| CASCI | free | 3d4s / MOC | (3,6) | Ar·pur | 90 | −759.7380 | −0.6 | −20673.5355 | −0.0176 | — | 1.4467 | 1e-01 | 2e+00 | 4e+00 | 1.5e+00 | 6.5e+00 | m-incomplete |
| CASSCF | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7606 | −23.2 | −20674.1498 | −0.6319 | 1e-5/189 | 1.9283 | — | 2e+00 | 6e+00 | 1.3e+00 | 9.7e+00 | swap; ; local minimum |
| CASSCF, frozen | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7338 | +3.6 | −20673.4202 | +0.0977 | 0/110 | 6.0000 | — | 2e-14 | 2e-14 | 1.3e-14 | 5.6e-14 |  |
| CASSCF, frozen core | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7326 | +4.8 | −20673.3886 | +0.1293 | 1e-10/100 | 6.0001 | — | 2e-02 | 3e-02 | 9.5e-02 | 5.8e-01 | not converged |
| SA-CASSCF | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7497 | −12.3 | −20673.8533 | −0.3355 | 1e-6/192 | 6.0047 | 8e-03 | 5e-01 | 2e+00 | 3.1e-01 | 4.0e+00 | swap |
| SA-CASSCF, frozen | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7327 | +4.6 | −20673.3915 | +0.1263 | 0/104 | 6.0000 | 1e-06 | 2e-14 | 2e-14 | 8.2e-15 | 6.5e-14 |  |
| SA-CASSCF, frozen core | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7326 | +4.8 | −20673.3884 | +0.1295 | 1e-13/101 | 6.0000 | 4e-06 | 9e-04 | 1e-03 | 4.6e-04 | 3.5e-03 | not converged |
| CASSCF | free | 3d4s@r1 / D | (3,6) | Ar·pur | 90 | −759.7600 | −22.6 | −20674.1321 | −0.6142 | 1e-5/196 | 1.6180 | — | 5e+00 | 4e+00 | 7.8e-01 | 9.1e+00 | swap; ; local minimum; bands 3d4s |
| CASSCF, frozen | free | 3d4s@r1 / D | (3,6) | Ar·pur | 90 | −759.7338 | +3.6 | −20673.4208 | +0.0971 | 0/111 | 6.0000 | — | 2e-14 | 2e-14 | 1.6e-14 | 6.9e-14 | bands 3d4s |
| CASSCF, frozen core | free | 3d4s@r1 / D | (3,6) | Ar·pur | 90 | −759.7313 | +6.1 | −20673.3521 | +0.1658 | 1e-5/155 | 6.0428 | — | 1e-01 | 2e-01 | 1.7e+03 | 7.2e+03 | not converged; ; bands 3d4s |
| SA-CASSCF | free | 3d4s@r1 / D | (3,6) | Ar·pur | 90 | −759.7490 | −11.7 | −20673.8352 | −0.3173 | 1e-5/187 | 6.0037 | 8e-03 | 1e+00 | 1e+00 | 8.5e-01 | 4.0e+00 | swap; ; bands 3d4s |
| SA-CASSCF, frozen | free | 3d4s@r1 / D | (3,6) | Ar·pur | 90 | −759.7328 | +4.6 | −20673.3921 | +0.1257 | 0/105 | 6.0000 | 1e-06 | 2e-14 | 2e-14 | 1.1e-14 | 8.3e-14 | bands 3d4s |
| SA-CASSCF, frozen core | free | 3d4s@r1 / D | (3,6) | Ar·pur | 90 | −759.7328 | +4.6 | −20673.3925 | +0.1254 | 0/100 | 6.0000 | 6e-06 | 2e-03 | 2e-03 | 1.5e-03 | 8.0e-03 | not converged; ; bands 3d4s |
| CASSCF | free | 3d4s@s1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s@s1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen core | free | 3d4s@s1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s@s1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s@s1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen core | free | 3d4s@s1 / D | (3,6) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7687 | −31.3 | −20674.3699 | −0.8520 | 1e-8/18 | 0.0541 | — | 2e+00 | 3e+00 | 2.8e-01 | 7.2e-01 | m-incomplete |
| CASSCF, frozen | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7466 | −9.2 | −20673.7680 | −0.2501 | 1e-12/6 | 0.0060 | — | 2e+00 | 3e+00 | 1.0e-01 | 2.6e-01 | m-incomplete |
| CASSCF, frozen core | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7658 | −28.4 | −20674.2913 | −0.7734 | 0/114 | 0.0515 | — | 2e+00 | 3e+00 | 1.6e+00 | 3.9e+00 | not converged; m-incomplete |
| SA-CASSCF | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7636 | −26.3 | −20674.2323 | −0.7145 | 1e-8/26 | 0.0310 | 7e-02 | 3e+00 | 4e+00 | 5.0e-01 | 1.2e+00 | m-incomplete |
| SA-CASSCF, frozen | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7465 | −9.1 | −20673.7656 | −0.2477 | 1e-11/6 | 0.0061 | 8e-02 | 2e+00 | 3e+00 | 1.0e-01 | 2.6e-01 | m-incomplete |
| SA-CASSCF, frozen core | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.4737 | +263.7 | −20666.3433 | +7.1746 | 1e-2/103 | 2.4610 | 8e-02 | 2e+00 | 3e+00 | 2.8e+01 | 2.0e+02 | not converged; m-incomplete |
| CASSCF | free | 3d4s / MOC | (3,6) | Ar·pur | 90 | −759.7602 | −22.8 | −20674.1379 | −0.6201 | 1e-4/188 | 0.3328 | — | 4e+00 | 5e+00 | 7.9e-01 | 4.3e+00 | m-incomplete; local minimum |
| CASSCF, frozen | free | 3d4s / MOC | (3,6) | Ar·pur | 90 | −759.7391 | −1.7 | −20673.5640 | −0.0461 | 1e-8/111 | 1.4467 | — | 2e+00 | 4e+00 | 1.5e+00 | 6.5e+00 | m-incomplete |
| CASSCF, frozen core | free | 3d4s / MOC | (3,6) | Ar·pur | 90 | −759.7346 | +2.7 | −20673.4433 | +0.0746 | 1e-4/120 | 1.5542 | — | 2e+00 | 4e+00 | 2.6e+03 | 9.6e+03 | not converged; m-incomplete |
| SA-CASSCF | free | 3d4s / MOC | (3,6) | Ar·pur | 90 | −759.7361 | +1.2 | −20673.4842 | +0.0337 | 1e-4/189 | 4.2112 | 7e-02 | 2e+00 | 3e+00 | 1.2e+00 | 4.3e+00 | m-incomplete |
| SA-CASSCF, frozen | free | 3d4s / MOC | (3,6) | Ar·pur | 90 | −759.7390 | −1.7 | −20673.5630 | −0.0452 | 1e-7/112 | 1.4467 | 1e-01 | 2e+00 | 4e+00 | 1.5e+00 | 6.5e+00 | m-incomplete |
| SA-CASSCF, frozen core | free | 3d4s / MOC | (3,6) | Ar·pur | 90 | −759.7322 | +5.2 | −20673.3760 | +0.1418 | 1e-5/143 | 1.6815 | 1e-01 | 2e+00 | 3e+00 | 6.5e+01 | 3.2e+02 | not converged; m-incomplete |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d4s / D | (3,6) | Ar·lane | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d4s / D | (3,6) | Ar·lane | 90 | −759.7323 | +5.0 | −20673.3805 | +0.1374 | — | 6.0000 | 9e-02 | 4e-11 | 5e-11 | 6.7e-15 | 4.4e-14 |  |
| CASSCF | free | 3d4s / D | (3,6) | Ar·lane | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s / D | (3,6) | Ar·lane | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | pinned | 3d4s / D | (3,6) | Ar·lane | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === |
| CASCI | free | 3d4s@s2 / D | (3,12) | Ar·raw | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s@s2 / D | (3,12) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s@s2 / D | (3,12) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s@s2 / D | (3,12) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen core | free | 3d4s@s2 / D | (3,12) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s@s2 / D | (3,12) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s@s2 / D | (3,12) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen core | free | 3d4s@s2 / D | (3,12) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s@s2 / D | (3,12) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d4s@s2 / D | (3,12) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d4s@s2 / D | (3,12) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen core | free | 3d4s@s2 / D | (3,12) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d4s@s2 / D | (3,12) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d4s@s2 / D | (3,12) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen core | free | 3d4s@s2 / D | (3,12) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === |
| CASCI | free | 3d / D | (3,5) | Ar·raw | 50 | −759.3360 | +401.3 | −20662.5968 | +10.9210 | — | 11.9994 | 5e-02 | 3e-02 | 4e-02 | 3.6e-02 | 2.1e-01 |  |
| CASCI | free | 3d@r1 / D | (3,5) | Ar·raw | 50 | −759.3375 | +399.9 | −20662.6368 | +10.8810 | — | 11.9994 | 5e-02 | 3e-02 | 4e-02 | 3.6e-02 | 2.1e-01 | bands 3d |
| CASCI | free | 3d@s1 / D | (3,5) | Ar·raw | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d / D | (3,5) | Ar | 50 | −759.3360 | +401.3 | −20662.5968 | +10.9210 | — | 11.9994 | 5e-02 | 3e-02 | 4e-02 | 1.2e-02 | 7.2e-02 |  |
| CASCI | free | 3d@r1 / D | (3,5) | Ar | 50 | −759.3375 | +399.9 | −20662.6368 | +10.8810 | — | 11.9994 | 5e-02 | 3e-02 | 4e-02 | 1.2e-02 | 7.2e-02 | bands 3d |
| CASCI | free | 3d@s1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d / MO | (3,5) | Ar | 50 | −759.7466 | −9.2 | −20673.7694 | −0.2515 | — | 0.0057 | 1e-01 | 2e+00 | 3e+00 | 8.4e-02 | 2.2e-01 | m-incomplete |
| CASCI | free | 3d / MOC | (3,5) | Ar | 50 | −759.7391 | −1.7 | −20673.5649 | −0.0470 | — | 0.7231 | 1e-01 | 3e+00 | 4e+00 | 9.4e-01 | 2.8e+00 | m-incomplete |
| CASSCF | free | 3d / D | (3,5) | Ar | 50 | −759.5697 | +167.7 | −20668.9548 | +4.5631 | 1e-5/135 | 10.8315 | — | 1e+00 | 2e+00 | 7.5e-01 | 4.3e+00 | local minimum |
| CASSCF, frozen | free | 3d / D | (3,5) | Ar | 50 | −759.4254 | +312.0 | −20665.0275 | +8.4904 | 1e-6/149 | 11.9994 | — | 3e-02 | 4e-02 | 1.2e-02 | 7.2e-02 |  |
| CASSCF, frozen core | free | 3d / D | (3,5) | Ar | 50 | −759.4394 | +298.0 | −20665.4082 | +8.1097 | 1e-4/109 | 11.8869 | — | 3e-01 | 4e-01 | 1.2e+02 | 5.7e+02 | not converged |
| SA-CASSCF | free | 3d / D | (3,5) | Ar | 50 | −759.5672 | +170.2 | −20668.8869 | +4.6309 | 1e-4/119 | 11.7177 | 2e-03 | 6e-01 | 8e-01 | 7.0e-02 | 4.2e-01 |  |
| SA-CASSCF, frozen | free | 3d / D | (3,5) | Ar | 50 | −759.4243 | +313.1 | −20664.9978 | +8.5201 | 1e-7/138 | 11.9994 | 1e-03 | 3e-02 | 4e-02 | 1.2e-02 | 7.2e-02 |  |
| SA-CASSCF, frozen core | free | 3d / D | (3,5) | Ar | 50 | −759.4308 | +306.6 | −20665.1754 | +8.3425 | 1e-3/140 | 11.7866 | 5e-02 | 2e-01 | 3e-01 | 9.8e+00 | 7.8e+01 | not converged |
| CASSCF | free | 3d@r1 / D | (3,5) | Ar | 50 | −759.5609 | +176.4 | −20668.7167 | +4.8012 | 1e-3/117 | 10.6726 | — | 1e+00 | 2e+00 | 5.6e-01 | 3.0e+00 | not converged; ; local minimum; bands 3d |
| CASSCF, frozen | free | 3d@r1 / D | (3,5) | Ar | 50 | −759.4261 | +311.3 | −20665.0468 | +8.4711 | 1e-5/176 | 11.9994 | — | 3e-02 | 4e-02 | 1.2e-02 | 7.2e-02 | not converged; ; bands 3d |
| CASSCF, frozen core | free | 3d@r1 / D | (3,5) | Ar | 50 | −759.4370 | +300.4 | −20665.3430 | +8.1749 | 1e-3/128 | 11.7751 | — | 4e-01 | 6e-01 | 3.3e+02 | 2.6e+03 | not converged; ; bands 3d |
| SA-CASSCF | free | 3d@r1 / D | (3,5) | Ar | 50 | −759.5613 | +176.1 | −20668.7253 | +4.7926 | 1e-3/107 | 11.7419 | 3e-03 | 6e-01 | 8e-01 | 4.8e-02 | 2.5e-01 | not converged; ; bands 3d |
| SA-CASSCF, frozen | free | 3d@r1 / D | (3,5) | Ar | 50 | −759.4251 | +312.3 | −20665.0198 | +8.4981 | 1e-5/175 | 11.9994 | 1e-03 | 3e-02 | 4e-02 | 1.2e-02 | 7.2e-02 | not converged; ; bands 3d |
| SA-CASSCF, frozen core | free | 3d@r1 / D | (3,5) | Ar | 50 | −759.4551 | +282.3 | −20665.8367 | +7.6812 | 1e-2/119 | 11.9330 | 2e-02 | 3e-01 | 4e-01 | 2.9e-01 | 1.4e+00 | not converged; ; bands 3d |
| CASSCF | free | 3d@s1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d@s1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen core | free | 3d@s1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d@s1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d@s1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen core | free | 3d@s1 / D | (3,5) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d / MO | (3,5) | Ar | 50 | −759.7674 | −30.1 | −20674.3359 | −0.8180 | 1e-9/14 | 0.0527 | — | 2e+00 | 3e+00 | 2.3e-01 | 5.9e-01 | m-incomplete |
| CASSCF, frozen | free | 3d / MO | (3,5) | Ar | 50 | −759.7466 | −9.2 | −20673.7694 | −0.2515 | 0/5 | 0.0057 | — | 2e+00 | 3e+00 | 8.4e-02 | 2.2e-01 | m-incomplete |
| CASSCF, frozen core | free | 3d / MO | (3,5) | Ar | 50 | −759.4767 | +260.7 | −20666.4251 | +7.0928 | 0/103 | 1.9796 | — | 2e+00 | 3e+00 | 8.3e-01 | 2.3e+00 | not converged; m-incomplete |
| SA-CASSCF | free | 3d / MO | (3,5) | Ar | 50 | −759.7628 | −25.4 | −20674.2089 | −0.6910 | 1e-8/16 | 0.0332 | 7e-02 | 2e+00 | 3e+00 | 1.6e-01 | 4.2e-01 | m-incomplete |
| SA-CASSCF, frozen | free | 3d / MO | (3,5) | Ar | 50 | −759.7465 | −9.1 | −20673.7664 | −0.2485 | 1e-13/6 | 0.0058 | 8e-02 | 2e+00 | 3e+00 | 8.4e-02 | 2.2e-01 | m-incomplete |
| SA-CASSCF, frozen core | free | 3d / MO | (3,5) | Ar | 50 | −759.7562 | −18.8 | −20674.0307 | −0.5128 | 0/103 | 0.0185 | 7e-02 | 2e+00 | 3e+00 | 1.2e-01 | 3.2e-01 | not converged; m-incomplete |
| CASSCF | free | 3d / MOC | (3,5) | Ar | 50 | −759.7674 | −30.1 | −20674.3359 | −0.8180 | 1e-8/54 | 0.0527 | — | 2e+00 | 3e+00 | 2.4e-01 | 5.9e-01 | m-incomplete |
| CASSCF, frozen | free | 3d / MOC | (3,5) | Ar | 50 | −759.7391 | −1.7 | −20673.5649 | −0.0470 | 1e-8/2 | 0.7231 | — | 3e+00 | 4e+00 | 9.4e-01 | 2.8e+00 | m-incomplete |
| CASSCF, frozen core | free | 3d / MOC | (3,5) | Ar | 50 | −759.7431 | −5.7 | −20673.6734 | −0.1555 | 0/131 | 1.5612 | — | 2e+00 | 3e+00 | 7.4e-01 | 2.3e+00 | not converged; m-incomplete |
| SA-CASSCF | free | 3d / MOC | (3,5) | Ar | 50 | −759.7365 | +0.9 | −20673.4940 | +0.0239 | 1e-12/30 | 5.0033 | 5e-02 | 2e+00 | 3e+00 | 4.4e-01 | 2.2e+00 | m-incomplete |
| SA-CASSCF, frozen | free | 3d / MOC | (3,5) | Ar | 50 | −759.7391 | −1.7 | −20673.5637 | −0.0458 | 1e-13/6 | 0.7231 | 1e-01 | 3e+00 | 4e+00 | 9.4e-01 | 2.8e+00 | m-incomplete |
| SA-CASSCF, frozen core | free | 3d / MOC | (3,5) | Ar | 50 | −759.7331 | +4.3 | −20673.4006 | +0.1172 | 1e-12/103 | 2.5575 | 9e-02 | 2e+00 | 3e+00 | 4.0e+03 | 1.8e+04 | not converged; m-incomplete |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d / D | (3,5) | Ar·pur | 50 | −759.3359 | +401.5 | −20662.5931 | +10.9248 | — | 12.0000 | 5e-02 | 2e-14 | 2e-14 | 7.5e-15 | 7.3e-14 |  |
| CASCI | free | 3d@r1 / D | (3,5) | Ar·pur | 50 | −759.3374 | +400.0 | −20662.6332 | +10.8847 | — | 12.0000 | 5e-02 | 2e-14 | 2e-14 | 6.6e-15 | 6.5e-14 | bands 3d |
| CASCI | free | 3d@s1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d / MO | (3,5) | Ar·pur | 50 | −759.7452 | −7.8 | −20673.7294 | −0.2116 | — | 0.0058 | 1e-01 | 2e+00 | 3e+00 | 8.4e-02 | 2.2e-01 | m-incomplete |
| CASCI | free | 3d / MOC | (3,5) | Ar·pur | 50 | −759.7377 | −0.3 | −20673.5251 | −0.0073 | — | 0.7233 | 1e-01 | 3e+00 | 4e+00 | 9.4e-01 | 2.8e+00 | m-incomplete |
| CASSCF | free | 3d / D | (3,5) | Ar·pur | 50 | −759.5870 | +150.4 | −20669.4258 | +4.0921 | 1e-4/192 | 0.2292 | — | 4e+00 | 5e+00 | 5.9e-01 | 1.5e+00 | local minimum |
| CASSCF, frozen | free | 3d / D | (3,5) | Ar·pur | 50 | −759.4257 | +311.7 | −20665.0365 | +8.4813 | 1e-7/173 | 12.0000 | — | 2e-14 | 2e-14 | 8.4e-15 | 5.9e-14 |  |
| CASSCF, frozen core | free | 3d / D | (3,5) | Ar·pur | 50 | −759.4014 | +336.0 | −20664.3742 | +9.1437 | 1e-3/109 | 10.8522 | — | 4e-01 | 9e-01 | 6.6e-01 | 5.6e+00 | not converged |
| SA-CASSCF | free | 3d / D | (3,5) | Ar·pur | 50 | −759.5715 | +165.9 | −20669.0027 | +4.5152 | 1e-5/180 | 11.9328 | 7e-04 | 3e-01 | 5e-01 | 2.1e-02 | 9.7e-02 |  |
| SA-CASSCF, frozen | free | 3d / D | (3,5) | Ar·pur | 50 | −759.4245 | +312.8 | −20665.0052 | +8.5127 | 1e-6/172 | 12.0000 | 6e-04 | 2e-14 | 2e-14 | 1.1e-14 | 6.5e-14 |  |
| SA-CASSCF, frozen core | free | 3d / D | (3,5) | Ar·pur | 50 | −759.4637 | +273.7 | −20666.0703 | +7.4476 | 1e-3/103 | 11.3322 | 6e-02 | 5e-01 | 5e-01 | 5.3e+01 | 9.5e+02 | not converged |
| CASSCF | free | 3d@r1 / D | (3,5) | Ar·pur | 50 | −759.5690 | +168.3 | −20668.9372 | +4.5807 | 1e-4/190 | 0.2093 | — | 4e+00 | 6e+00 | 9.6e+01 | 2.2e+02 | not converged; ; local minimum; bands 3d |
| CASSCF, frozen | free | 3d@r1 / D | (3,5) | Ar·pur | 50 | −759.4268 | +310.6 | −20665.0654 | +8.4524 | 1e-7/173 | 12.0000 | — | 2e-14 | 2e-14 | 8.3e-15 | 8.9e-14 | not converged; ; bands 3d |
| CASSCF, frozen core | free | 3d@r1 / D | (3,5) | Ar·pur | 50 | −759.3719 | +365.5 | −20663.5730 | +9.9449 | 1e-4/107 | 11.9131 | — | 2e-01 | 4e-01 | 4.8e-01 | 2.6e+00 | not converged; ; bands 3d |
| SA-CASSCF | free | 3d@r1 / D | (3,5) | Ar·pur | 50 | −759.5715 | +165.9 | −20669.0032 | +4.5147 | 1e-5/175 | 11.9423 | 7e-04 | 5e-01 | 4e-01 | 1.2e-02 | 6.0e-02 | not converged; ; bands 3d |
| SA-CASSCF, frozen | free | 3d@r1 / D | (3,5) | Ar·pur | 50 | −759.4256 | +311.8 | −20665.0338 | +8.4840 | 1e-6/169 | 12.0000 | 6e-04 | 2e-14 | 2e-14 | 7.5e-15 | 7.0e-14 | not converged; ; bands 3d |
| SA-CASSCF, frozen core | free | 3d@r1 / D | (3,5) | Ar·pur | 50 | −759.4363 | +301.1 | −20665.3252 | +8.1927 | 1e-5/133 | 11.8208 | 5e-02 | 4e-01 | 4e-01 | 2.6e+03 | 1.3e+04 | not converged; ; bands 3d |
| CASSCF | free | 3d@s1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d@s1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen core | free | 3d@s1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d@s1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d@s1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen core | free | 3d@s1 / D | (3,5) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d / MO | (3,5) | Ar·pur | 50 | −759.7674 | −30.1 | −20674.3359 | −0.8180 | 1e-10/14 | 0.0527 | — | 2e+00 | 3e+00 | 2.3e-01 | 5.9e-01 | m-incomplete |
| CASSCF, frozen | free | 3d / MO | (3,5) | Ar·pur | 50 | −759.7462 | −8.8 | −20673.7579 | −0.2400 | 0/6 | 0.0058 | — | 2e+00 | 3e+00 | 8.4e-02 | 2.2e-01 | m-incomplete |
| CASSCF, frozen core | free | 3d / MO | (3,5) | Ar·pur | 50 | −759.3798 | +357.6 | −20663.7880 | +9.7299 | 1e-11/103 | 1.8699 | — | 2e+00 | 3e+00 | 1.7e+00 | 4.2e+00 | not converged; m-incomplete |
| SA-CASSCF | free | 3d / MO | (3,5) | Ar·pur | 50 | −759.7628 | −25.4 | −20674.2089 | −0.6911 | 1e-8/33 | 0.0332 | 7e-02 | 2e+00 | 3e+00 | 1.6e-01 | 4.2e-01 | m-incomplete |
| SA-CASSCF, frozen | free | 3d / MO | (3,5) | Ar·pur | 50 | −759.7461 | −8.7 | −20673.7549 | −0.2370 | 1e-12/6 | 0.0058 | 8e-02 | 2e+00 | 3e+00 | 8.4e-02 | 2.2e-01 | m-incomplete |
| SA-CASSCF, frozen core | free | 3d / MO | (3,5) | Ar·pur | 50 | −759.7531 | −15.7 | −20673.9444 | −0.4265 | 0/103 | 0.0174 | 7e-02 | 2e+00 | 3e+00 | 1.2e-01 | 3.0e-01 | not converged; m-incomplete |
| CASSCF | free | 3d / MOC | (3,5) | Ar·pur | 50 | −759.7500 | −12.6 | −20673.8617 | −0.3438 | 1e-4/186 | 0.0264 | — | 4e+00 | 5e+00 | 1.2e+00 | 2.4e+00 | m-incomplete; local minimum |
| CASSCF, frozen | free | 3d / MOC | (3,5) | Ar·pur | 50 | −759.7387 | −1.3 | −20673.5536 | −0.0357 | 1e-8/111 | 0.7233 | — | 3e+00 | 4e+00 | 9.4e-01 | 2.8e+00 | m-incomplete |
| CASSCF, frozen core | free | 3d / MOC | (3,5) | Ar·pur | 50 | −759.7368 | +0.6 | −20673.5023 | +0.0155 | 1e-5/159 | 0.7452 | — | 3e+00 | 4e+00 | 8.3e+03 | 3.5e+04 | not converged; m-incomplete |
| SA-CASSCF | free | 3d / MOC | (3,5) | Ar·pur | 50 | −759.7305 | +6.8 | −20673.3317 | +0.1862 | 1e-4/163 | 3.8280 | 6e-02 | 3e+00 | 4e+00 | 1.7e+00 | 6.1e+00 | m-incomplete |
| SA-CASSCF, frozen | free | 3d / MOC | (3,5) | Ar·pur | 50 | −759.7387 | −1.3 | −20673.5525 | −0.0346 | 1e-7/115 | 0.7233 | 1e-01 | 3e+00 | 4e+00 | 9.4e-01 | 2.8e+00 | m-incomplete |
| SA-CASSCF, frozen core | free | 3d / MOC | (3,5) | Ar·pur | 50 | −759.7276 | +9.8 | −20673.2511 | +0.2668 | 1e-4/118 | 0.9089 | 1e-01 | 3e+00 | 4e+00 | 1.5e+02 | 5.0e+02 | not converged; m-incomplete |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d / D | (3,5) | Ar·lane | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASCI | free | 3d / D | (3,5) | Ar·lane | 50 | −759.3352 | +402.1 | −20662.5749 | +10.9430 | — | 12.0000 | 5e-02 | 2e-14 | 2e-14 | 6.3e-15 | 5.5e-14 |  |
| CASSCF | free | 3d / D | (3,5) | Ar·lane | 50 | −759.5714 | +166.0 | −20669.0010 | +4.5169 | 1e-5/159 | 10.6975 | — | 1e+00 | 2e+00 | 4.6e+00 | 2.7e+01 | not converged; ; local minimum |
| SA-CASSCF | free | 3d / D | (3,5) | Ar·lane | 50 | −759.5709 | +166.5 | −20668.9882 | +4.5296 | 1e-5/142 | 11.8667 | 8e-04 | 6e-01 | 6e-01 | 4.0e-02 | 2.2e-01 | not converged |
| SA-CASSCF | pinned | 3d / D | (3,5) | Ar·lane | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === | === |
| CASCI | free | 3d@s2 / D | (3,10) | Ar·raw | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d@s2 / D | (3,10) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d@s2 / D | (3,10) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d@s2 / D | (3,10) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen core | free | 3d@s2 / D | (3,10) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d@s2 / D | (3,10) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d@s2 / D | (3,10) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen core | free | 3d@s2 / D | (3,10) | Ar | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d@s2 / D | (3,10) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF | free | 3d@s2 / D | (3,10) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen | free | 3d@s2 / D | (3,10) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| CASSCF, frozen core | free | 3d@s2 / D | (3,10) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF | free | 3d@s2 / D | (3,10) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen | free | 3d@s2 / D | (3,10) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
| SA-CASSCF, frozen core | free | 3d@s2 / D | (3,10) | Ar·pur | — | — | — | — | — | — | — | — | — | — | — | — | pending |
<!-- methods:end -->
