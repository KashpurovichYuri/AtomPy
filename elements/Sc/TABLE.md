# Sc (4s^2 3d^1)

Review stage: the 3d4s window only -- CAS(3,6), plus the `@s2`
prior's CAS(3,12) (the budget stated per row); the other spans
rejoin as review proceeds. ΔE's absolute zero is the ROHF anchor
row (the UHF anchor follows). One block order in both tiers: the
D-family by (k, source), then the controls.

The starts: `D` = the reference's radial content in AO form -- on
this CONTRACTED basis the authors' 1:1 contraction, the equal-budget
baseline; `@r1` = the head configuration's density block (the head
is this table's own ROHF); `@s1` = the configuration-series mean;
`@s2` = the mean plus the configuration spread (the second radial)
at double CAS; `MO`/`MOC` = the SCF's canonical controls
(m-incomplete by construction); `NO` = the natural-orbital start
(its optimizer legs render; its CI row is energy-degenerate with D
and stays store-only).

The verdicts (2026-11-11 legs): FREE is start-blind EXACTLY -- `@r1`
and `@s1` both land the incumbent's own solution (−31.7, the same
basin to six decimals, 51–52 cycles beside its 54): a k=1 prior's
value is its START, never a new endpoint. SA: `@r1` reproduces the
incumbent's basin (−13.3) in 34 cycles against 129 -- the faster
road to the same state; `@s1` lands a cleaner-but-higher basin
(−3.5; its stored stack consistent, against the 0.8 Eh stack
anomaly the incumbent/`@r1` basin carries -- the comm columns state
it); `@s2` free is NOT CONVERGED (the note flags it; −37.8 is a cap
reading) and its SA (−11.9) sits between. MOC's free leg catches
the free tier (−31.5) but its SA stalls at its CI level (−1.5) --
the term average cannot repair the character-selected start. The
radial mask is undeclared on prior windows BY CONSTRUCTION (its
pair predicate needs member-labeled basis rows). The raw core group
is the anisotropy detector; the pur group the closure certificates;
the k evidence rides the payloads (`window_report`) and `logs/`.

<!-- methods:begin -->
| method | route | window | CAS | core | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −759.7363 | 0 | −20673.4871 | +0.0000 | — | — | — | — | — | — | — | forced filling |
| UHF (anchor) | — | — | — | — | — | −759.7379 | −1.6 | −20673.5317 | −0.0445 | — | — | — | — | — | — | — | ⟨S²⟩=0.7620, σ_min=0.9945 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (3,6) | Ar·raw | 90 | −759.7296 | +6.7 | −20673.3051 | +0.1820 | — | 6.0003 | 6e-02 | 2e-03 | 2e-02 | 2.8e-02 | 1.3e-01 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (3,6) | Ar | 90 | −759.7296 | +6.7 | −20673.3051 | +0.1820 | — | 5.9998 | 6e-02 | 2e-03 | 2e-02 | 8.9e-04 | 4.0e-02 |  |
| CASCI | free | 3d4s@r1 / D | (3,6) | Ar | 90 | −759.7366 | −0.4 | −20673.4978 | −0.0106 | — | 5.9997 | 9e-02 | 2e-03 | 2e-02 | 1.0e-03 | 4.7e-02 |  |
| CASCI | free | 3d4s@s1 / D | (3,6) | Ar | 90 | −759.7240 | +12.3 | −20673.1526 | +0.3346 | — | 5.9998 | 3e-02 | 2e-03 | 2e-02 | 6.8e-04 | 3.1e-02 |  |
| CASCI | free | 3d4s@s2 / D | (3,12) | Ar | 792 | −759.7460 | −9.8 | −20673.7527 | −0.2656 | — | 5.9998 | 5e-02 | 2e-03 | 3e-02 | 2.8e-03 | 1.4e-01 |  |
| CASCI | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7470 | −10.7 | −20673.7795 | −0.2924 | — | 0.0062 | 1e-01 | 3e+00 | 4e+00 | 1.5e-01 | 7.2e-01 | control: m-incomplete by construction |
| CASCI | free | 3d4s / MOC | (3,6) | Ar | 90 | −759.7385 | −2.3 | −20673.5487 | −0.0616 | — | 3.3532 | 2e-01 | 3e+00 | 4e+00 | 1.4e+00 | 1.6e+00 | control: m-incomplete by construction |
| CASSCF | free | 3d4s / D | (3,6) | Ar | 90 | −759.7680 | −31.7 | −20674.3496 | −0.8625 | 1e-8/54 | 0.0469 | — | 4e+00 | 4e+00 | 5.7e-01 | 1.2e+00 |  |
| CASSCF, virtual-mix | free | 3d4s / D | (3,6) | Ar | 90 | −759.7676 | −31.3 | −20674.3397 | −0.8525 | 1e-8/73 | 0.0465 | — | 4e+00 | 4e+00 | 5.9e-01 | 1.2e+00 |  |
| CASSCF, radial | free | 3d4s / D | (3,6) | Ar | 90 | −759.7431 | −6.9 | −20673.6741 | −0.1870 | 1e-12/15 | 3.5331 | — | 3e+00 | 2e+00 | 1.6e+00 | 4.0e+00 |  |
| CASSCF, frozen | free | 3d4s / D | (3,6) | Ar | 90 | −759.7297 | +6.6 | −20673.3075 | +0.1796 | 1e-11/5 | 5.9998 | — | 2e-03 | 2e-02 | 8.9e-04 | 4.0e-02 |  |
| SA-CASSCF | free | 3d4s / D | (3,6) | Ar | 90 | −759.7495 | −13.3 | −20673.8478 | −0.3607 | 1e-12/129 | 6.0020 | 8e-03 | 1e+00 | 2e+00 | 7.9e-01 | 3.8e+00 |  |
| SA-CASSCF, radial | free | 3d4s / D | (3,6) | Ar | 90 | −759.7370 | −0.8 | −20673.5085 | −0.0214 | 1e-12/9 | 5.9997 | 2e-03 | 2e-03 | 2e-02 | 1.2e-03 | 5.8e-02 |  |
| SA-CASSCF, frozen | free | 3d4s / D | (3,6) | Ar | 90 | −759.7291 | +7.2 | −20673.2924 | +0.1947 | 1e-11/6 | 5.9998 | 5e-04 | 2e-03 | 2e-02 | 8.9e-04 | 4.0e-02 |  |
| CASSCF | free | 3d4s@r1 / D | (3,6) | Ar | 90 | −759.7680 | −31.7 | −20674.3496 | −0.8625 | 1e-8/51 | 0.0469 | — | 4e+00 | 4e+00 | 5.6e-01 | 1.2e+00 |  |
| SA-CASSCF | free | 3d4s@r1 / D | (3,6) | Ar | 90 | −759.7495 | −13.3 | −20673.8479 | −0.3607 | 1e-8/34 | 6.0020 | 8e-03 | 1e+00 | 1e+00 | 8.4e-01 | 3.8e+00 |  |
| CASSCF | free | 3d4s@s1 / D | (3,6) | Ar | 90 | −759.7680 | −31.7 | −20674.3496 | −0.8625 | 1e-8/52 | 0.0470 | — | 4e+00 | 4e+00 | 5.6e-01 | 1.2e+00 |  |
| SA-CASSCF | free | 3d4s@s1 / D | (3,6) | Ar | 90 | −759.7397 | −3.5 | −20673.5819 | −0.0948 | 1e-8/65 | 6.0000 | 5e-07 | 3e-03 | 8e-03 | 6.6e-05 | 6.6e-04 |  |
| CASSCF | free | 3d4s@s2 / D | (3,12) | Ar | 792 | −759.7741 | −37.8 | −20674.5158 | −1.0286 | 1e-6/167 | 2.5991 | — | 5e+00 | 6e+00 | 7.6e+00 | 3.2e+01 | not converged;  |
| SA-CASSCF | free | 3d4s@s2 / D | (3,12) | Ar | 792 | −759.7482 | −11.9 | −20673.8119 | −0.3247 | 1e-9/16 | 6.0000 | 2e-07 | 3e-03 | 3e-03 | 3.3e-04 | 1.0e-02 |  |
| CASSCF | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7675 | −31.3 | −20674.3377 | −0.8506 | 1e-8/24 | 0.0533 | — | 3e+00 | 2e+00 | 3.8e-01 | 7.6e-01 | control: m-incomplete by construction |
| CASSCF, frozen | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7470 | −10.7 | −20673.7795 | −0.2924 | 1e-11/4 | 0.0062 | — | 3e+00 | 4e+00 | 1.5e-01 | 7.2e-01 | control: m-incomplete by construction |
| SA-CASSCF | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7626 | −26.3 | −20674.2032 | −0.7161 | 1e-11/14 | 0.0309 | 7e-02 | 3e+00 | 4e+00 | 3.1e-01 | 1.2e+00 | control: m-incomplete by construction |
| SA-CASSCF, frozen | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7469 | −10.7 | −20673.7774 | −0.2903 | 1e-10/4 | 0.0062 | 8e-02 | 3e+00 | 4e+00 | 1.5e-01 | 7.2e-01 | control: m-incomplete by construction |
| CASSCF | free | 3d4s / MOC | (3,6) | Ar | 90 | −759.7678 | −31.5 | −20674.3447 | −0.8575 | 1e-8/105 | 0.0526 | — | 2e+00 | 3e+00 | 3.5e-01 | 7.2e-01 | control: m-incomplete by construction |
| SA-CASSCF | free | 3d4s / MOC | (3,6) | Ar | 90 | −759.7378 | −1.5 | −20673.5282 | −0.0410 | 1e-7/135 | 5.0014 | 5e-02 | 2e+00 | 5e+00 | 1.1e+00 | 6.4e+00 | not converged; control: m-incomplete by construction |
| CASSCF | free | 3d4sno2 / NO | (3,6) | Ar | 2 | −759.7366 | −0.3 | −20673.4952 | −0.0081 | 0/152 | 0.0000 | — | 1e+00 | 4e+00 | 3.2e-06 | 3.6e-21 |  |
| CASSCF, frozen | free | 3d4sno2 / NO | (3,6) | Ar | 2 | −759.7266 | +9.6 | −20673.2251 | +0.2620 | 0/6 | 0.0000 | — | 3e-01 | 3e+00 | 5.9e-18 | 1.1e-32 |  |
| SA-CASSCF | free | 3d4sno2 / NO | (3,6) | Ar | 2 | −759.7366 | −0.3 | −20673.4952 | −0.0081 | 0/152 | 0.0000 | 4e+00 | 1e+00 | 4e+00 | 2.9e-06 | 9.1e-22 |  |
| SA-CASSCF, frozen | free | 3d4sno2 / NO | (3,6) | Ar | 2 | −759.7266 | +9.6 | −20673.2251 | +0.2620 | 0/6 | 0.0000 | 1e-01 | 3e-01 | 3e+00 | 5.9e-18 | 1.1e-32 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7289 | +7.4 | −20673.2867 | +0.2004 | — | 6.0000 | 6e-02 | 7e-15 | 9e-15 | 3.4e-14 | 1.7e-13 |  |
| CASCI | free | 3d4s / MOC | (3,6) | Ar·pur | 90 | −759.7377 | −1.5 | −20673.5276 | −0.0405 | — | 3.3532 | 2e-01 | 3e+00 | 4e+00 | 1.4e+00 | 1.6e+00 | control: m-incomplete by construction |
| CASSCF | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7678 | −31.5 | −20674.3447 | −0.8575 | 1e-8/95 | 0.0526 | — | 2e+00 | 2e+00 | 3.5e-01 | 7.2e-01 |  |
| CASSCF, radial | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7425 | −6.2 | −20673.6558 | −0.1686 | 1e-10/15 | 3.5384 | — | 3e+00 | 2e+00 | 1.6e+00 | 4.0e+00 |  |
| CASSCF, frozen | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7295 | +6.7 | −20673.3042 | +0.1830 | 0/6 | 6.0000 | — | 7e-15 | 9e-15 | 3.4e-14 | 1.7e-13 |  |
| SA-CASSCF | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7397 | −3.5 | −20673.5819 | −0.0948 | 1e-11/37 | 6.0000 | 3e-07 | 1e-04 | 2e-04 | 4.6e-05 | 3.0e-04 |  |
| SA-CASSCF, radial | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7362 | +0.0 | −20673.4869 | +0.0002 | 1e-12/9 | 6.0000 | 7e-09 | 3e-07 | 4e-06 | 4.5e-07 | 2.0e-05 |  |
| SA-CASSCF, frozen | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7290 | +7.3 | −20673.2891 | +0.1980 | 0/6 | 6.0000 | 5e-08 | 7e-15 | 9e-15 | 1.6e-13 | 5.9e-13 |  |
| CASSCF | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7675 | −31.3 | −20674.3377 | −0.8506 | 1e-8/24 | 0.0533 | — | 3e+00 | 2e+00 | 3.8e-01 | 7.6e-01 | control: m-incomplete by construction |
| CASSCF, frozen | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7470 | −10.7 | −20673.7786 | −0.2915 | 1e-10/5 | 0.0062 | — | 3e+00 | 4e+00 | 1.5e-01 | 7.2e-01 | control: m-incomplete by construction |
| SA-CASSCF | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7626 | −26.3 | −20674.2032 | −0.7161 | 1e-10/13 | 0.0309 | 7e-02 | 3e+00 | 4e+00 | 3.1e-01 | 1.2e+00 | control: m-incomplete by construction |
| SA-CASSCF, frozen | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7469 | −10.6 | −20673.7765 | −0.2894 | 1e-11/6 | 0.0063 | 8e-02 | 3e+00 | 4e+00 | 1.5e-01 | 7.2e-01 | control: m-incomplete by construction |
| CASSCF | free | 3d4sno2 / NO | (3,6) | Ar·pur | 90 | −759.7678 | −31.5 | −20674.3447 | −0.8575 | 1e-8/50 | 0.0526 | — | 2e+00 | 3e+00 | 3.4e-01 | 7.2e-01 |  |
| CASSCF, frozen | free | 3d4sno2 / NO | (3,6) | Ar·pur | 90 | −759.7295 | +6.7 | −20673.3042 | +0.1830 | 1e-9/6 | 6.0000 | — | 7e-15 | 9e-15 | 3.3e-14 | 2.3e-13 |  |
| SA-CASSCF | free | 3d4sno2 / NO | (3,6) | Ar·pur | 90 | −759.7397 | −3.5 | −20673.5819 | −0.0948 | 1e-8/34 | 6.0000 | 2e-07 | 3e-03 | 5e-03 | 2.2e-05 | 1.3e-04 |  |
| SA-CASSCF, frozen | free | 3d4sno2 / NO | (3,6) | Ar·pur | 90 | −759.7290 | +7.3 | −20673.2891 | +0.1980 | 1e-13/6 | 6.0000 | 6e-08 | 7e-15 | 9e-15 | 9.6e-14 | 7.2e-13 |  |
<!-- methods:end -->
