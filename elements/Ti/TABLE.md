# Ti (3d^2 4s^2)

The prior windows (the `@r1`/`@s1`/`@s2` labels in the window
column — the
basis-preparation layer): the tag names the RADIAL PRIOR the window
distills its columns from. `r1` = the head configuration's own density
block, one radial per shell; `s1`/`s2` = the configuration-series
average (one SCF per valence configuration, occupation-forced) at one
and two radials per shell — `s2`'s second radial IS the configuration
spread, so its window is DOUBLE BUDGET (the CAS column states it:
CAS(4,12) beside the incumbent CAS(4,6); its dominance reads as the
headroom the second radial buys, never same-budget superiority). The
`3d` prior family runs at the d² sector (2,0) — the configuration the
span can host without an incumbent twin at the case default; its own
block, compared within and against the anchor. Untagged rows are the
incumbent: the reference-source distillation, which on this CONTRACTED
basis is the basis authors' 1:1 contraction — the equal-budget
baseline the priors compare against; the ROHF anchor row is ΔE's
absolute zero throughout. The k evidence (per shell: the band list,
the radial spectrum, the discarded fraction, the per-configuration
captures) rides the cells' payloads (`window_report`) and the
provenance logs (`logs/`).

<!-- methods:begin -->
| method | route | window | CAS | core | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −848.4065 | 0 | −23086.3283 | +0.0000 | — | — | — | — | — | — | — | forced filling |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d / D | (2,5) | Ar | 10 | −848.4001 | +6.4 | −23086.1531 | +0.1752 | — | 12.0000 | 6e-02 | 2e-07 | 3e-07 | 1.6e-14 | 9.1e-14 |  |
| CASCI | free | 3d@r1 / D | (2,5) | Ar | 10 | −848.4063 | +0.2 | −23086.3236 | +0.0047 | — | 12.0000 | 6e-02 | 2e-07 | 3e-07 | 8.8e-15 | 7.7e-14 |  |
| CASCI | free | 3d@s1 / D | (2,5) | Ar | 10 | −848.3890 | +17.4 | −23085.8536 | +0.4748 | — | 12.0000 | 5e-02 | 2e-07 | 3e-07 | 7.2e-15 | 8.7e-14 |  |
| CASCI | free | 3d@s2 / D | (2,10) | Ar | 45 | −848.4096 | −3.1 | −23086.4138 | −0.0854 | — | 12.0000 | 6e-02 | 2e-07 | 3e-07 | 5.3e-09 | 2.7e-08 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4p / D | (2,8) | Ar | 28 | −848.4001 | +6.4 | −23086.1531 | +0.1752 | — | 12.0000 | 6e-02 | 1e-03 | 2e-03 | 3.6e-05 | 1.0e-04 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4p4d / D | (2,13) | Ar | 78 | −848.4096 | −3.1 | −23086.4138 | −0.0854 | — | 12.0000 | 6e-02 | 1e-03 | 2e-03 | 1.3e-04 | 5.0e-04 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (4,6) | Ar·raw | 120 | −848.3959 | +10.5 | −23086.0413 | +0.2871 | — | 12.0000 | 5e-02 | 6e-09 | 7e-09 | 4.9e-02 | 2.6e-01 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (4,6) | Ar | 120 | −848.3959 | +10.5 | −23086.0413 | +0.2871 | — | 12.0000 | 5e-02 | 6e-09 | 7e-09 | 4.2e-09 | 2.3e-08 |  |
| CASCI | free | 3d4s@r1 / D | (4,6) | Ar | 120 | −848.4066 | −0.1 | −23086.3306 | −0.0023 | — | 12.0000 | 6e-02 | 6e-09 | 7e-09 | 4.6e-09 | 2.5e-08 |  |
| CASCI | free | 3d4s@s1 / D | (4,6) | Ar | 120 | −848.3890 | +17.5 | −23085.8515 | +0.4768 | — | 12.0000 | 2e-02 | 6e-09 | 7e-09 | 3.8e-09 | 2.0e-08 |  |
| CASCI | free | 3d4s@s2 / D | (4,12) | Ar | 2,640 | −848.4174 | −10.9 | −23086.6251 | −0.2968 | — | 12.0000 | 2e-02 | 6e-09 | 8e-09 | 2.1e-08 | 1.5e-07 |  |
| CASCI | free | 3d4sno2 / NO | (4,6) | Ar | 120 | −848.3959 | +10.5 | −23086.0413 | +0.2871 | — | 12.0000 | 5e-02 | 6e-09 | 7e-09 | 4.1e-09 | 2.2e-08 |  |
| CASCI | free | 3d4s / MOC | (4,6) | Ar | 120 | −848.4092 | −2.7 | −23086.4021 | −0.0738 | — | 2.0925 | 2e-01 | 3e+00 | 4e+00 | 2.5e+00 | 4.3e+00 | control: m-incomplete by construction |
| CASSCF | free | 3d4s / D | (4,6) | Ar | 120 | −848.4377 | −31.3 | −23087.1787 | −0.8504 | 1e-8/38 | 0.0503 | — | 3e+00 | 4e+00 | 4.3e-01 | 1.1e+00 |  |
| CASSCF, virtual-mix | free | 3d4s / D | (4,6) | Ar | 120 | −848.4377 | −31.2 | −23087.1779 | −0.8496 | 1e-9/42 | 0.0496 | — | 3e+00 | 4e+00 | 4.3e-01 | 1.1e+00 |  |
| CASSCF, radial | free | 3d4s / D | (4,6) | Ar | 120 | −848.4112 | −4.7 | −23086.4576 | −0.1292 | 1e-11/17 | 6.9323 | — | 3e+00 | 4e+00 | 3.2e+00 | 1.3e+01 |  |
| CASSCF, frozen | free | 3d4s / D | (4,6) | Ar | 120 | −848.3961 | +10.4 | −23086.0459 | +0.2825 | 1e-11/5 | 12.0000 | — | 6e-09 | 7e-09 | 4.2e-09 | 2.3e-08 |  |
| SA-CASSCF | free | 3d4s / D | (4,6) | Ar | 120 | −848.4166 | −10.1 | −23086.6025 | −0.2741 | 1e-9/33 | 12.0029 | 3e-03 | 1e+00 | 1e+00 | 1.1e+00 | 6.0e+00 |  |
| SA-CASSCF, radial | free | 3d4s / D | (4,6) | Ar | 120 | −848.4069 | −0.4 | −23086.3399 | −0.0115 | 1e-10/8 | 12.0000 | 4e-03 | 2e-03 | 2e-03 | 1.3e-03 | 6.9e-03 |  |
| SA-CASSCF, frozen | free | 3d4s / D | (4,6) | Ar | 120 | −848.3948 | +11.7 | −23086.0090 | +0.3194 | 1e-10/6 | 12.0000 | 2e-07 | 6e-09 | 7e-09 | 4.2e-09 | 2.3e-08 |  |
| CASSCF | free | 3d4s / MO | (4,6) | Ar | 120 | −848.4377 | −31.3 | −23087.1787 | −0.8504 | 1e-8/12 | 0.0503 | — | 3e+00 | 4e+00 | 4.3e-01 | 1.1e+00 | control: m-incomplete by construction |
| CASSCF, frozen | free | 3d4s / MO | (4,6) | Ar | 120 | −848.4171 | −10.6 | −23086.6177 | −0.2893 | 0/4 | 0.0051 | — | 3e+00 | 4e+00 | 1.6e-01 | 3.9e-01 | control: m-incomplete by construction |
| SA-CASSCF | free | 3d4s / MO | (4,6) | Ar | 120 | −848.4316 | −25.1 | −23087.0118 | −0.6835 | 1e-11/9 | 0.0255 | 9e-02 | 3e+00 | 4e+00 | 3.0e-01 | 7.4e-01 | control: m-incomplete by construction |
| SA-CASSCF, frozen | free | 3d4s / MO | (4,6) | Ar | 120 | −848.4170 | −10.5 | −23086.6151 | −0.2868 | 1e-10/4 | 0.0051 | 9e-02 | 3e+00 | 4e+00 | 1.6e-01 | 3.9e-01 | control: m-incomplete by construction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (4,6) | Ar·pur | 120 | −848.3946 | +11.9 | −23086.0044 | +0.3239 | — | 12.0000 | 5e-02 | 6e-15 | 8e-15 | 3.1e-14 | 5.4e-13 |  |
| CASCI | free | 3d4sno2 / NO | (4,6) | Ar·pur | 120 | −848.3946 | +11.9 | −23086.0044 | +0.3239 | — | 12.0000 | 5e-02 | 6e-15 | 8e-15 | 5.5e-14 | 3.8e-13 |  |
| CASCI | free | 3d4s / MOC | (4,6) | Ar·pur | 120 | −848.4078 | −1.3 | −23086.3639 | −0.0355 | — | 2.0925 | 2e-01 | 3e+00 | 4e+00 | 2.5e+00 | 4.3e+00 | control: m-incomplete by construction |
| CASSCF | free | 3d4s / D | (4,6) | Ar·pur | 120 | −848.4377 | −31.3 | −23087.1787 | −0.8504 | 1e-11/41 | 0.0503 | — | 3e+00 | 4e+00 | 4.3e-01 | 1.1e+00 |  |
| CASSCF, radial | free | 3d4s / D | (4,6) | Ar·pur | 120 | −848.4100 | −3.5 | −23086.4234 | −0.0951 | 0/12 | 6.9638 | — | 3e+00 | 4e+00 | 3.2e+00 | 1.2e+01 |  |
| CASSCF, frozen | free | 3d4s / D | (4,6) | Ar·pur | 120 | −848.3961 | +10.4 | −23086.0459 | +0.2825 | 1e-12/6 | 12.0000 | — | 6e-15 | 8e-15 | 8.9e-14 | 4.7e-13 |  |
| SA-CASSCF | free | 3d4s / D | (4,6) | Ar·pur | 120 | −848.4166 | −10.1 | −23086.6025 | −0.2741 | 1e-8/94 | 12.0029 | 3e-03 | 1e+00 | 2e+00 | 1.0e+00 | 6.0e+00 |  |
| SA-CASSCF, radial | free | 3d4s / D | (4,6) | Ar·pur | 120 | −848.4055 | +1.0 | −23086.3017 | +0.0266 | 1e-10/8 | 12.0000 | 6e-08 | 6e-05 | 4e-05 | 1.4e-04 | 5.7e-04 |  |
| SA-CASSCF, frozen | free | 3d4s / D | (4,6) | Ar·pur | 120 | −848.3948 | +11.7 | −23086.0090 | +0.3194 | 1e-11/6 | 12.0000 | 1e-07 | 6e-15 | 8e-15 | 4.8e-14 | 4.0e-13 |  |
| CASSCF | free | 3d4s / MO | (4,6) | Ar·pur | 120 | −848.4377 | −31.3 | −23087.1787 | −0.8504 | 1e-12/14 | 0.0503 | — | 3e+00 | 4e+00 | 4.3e-01 | 1.1e+00 | control: m-incomplete by construction |
| CASSCF, frozen | free | 3d4s / MO | (4,6) | Ar·pur | 120 | −848.4171 | −10.6 | −23086.6177 | −0.2893 | 0/5 | 0.0051 | — | 3e+00 | 4e+00 | 1.6e-01 | 3.9e-01 | control: m-incomplete by construction |
| SA-CASSCF | free | 3d4s / MO | (4,6) | Ar·pur | 120 | −848.4316 | −25.1 | −23087.0119 | −0.6835 | 1e-11/10 | 0.0255 | 9e-02 | 3e+00 | 4e+00 | 3.0e-01 | 7.4e-01 | control: m-incomplete by construction |
| SA-CASSCF, frozen | free | 3d4s / MO | (4,6) | Ar·pur | 120 | −848.4170 | −10.5 | −23086.6151 | −0.2868 | 1e-13/6 | 0.0051 | 9e-02 | 3e+00 | 4e+00 | 1.6e-01 | 3.9e-01 | control: m-incomplete by construction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d4s / D | (4,6) | Ar·lane | 120 | −848.3931 | +13.4 | −23085.9634 | +0.3649 | — | 12.0000 | 5e-15 | — | — | 2.7e-14 | 2.4e-13 | SCF core; offset = bookkeeping form |
| CASCI | free | 3d4s / D | (4,6) | Ar·lane | 120 | −848.3931 | +13.4 | −23085.9634 | +0.3649 | — | 12.0000 | 4e-02 | 6e-15 | 8e-15 | 2.8e-14 | 2.5e-13 |  |
| CASSCF | free | 3d4s / D | (4,6) | Ar·lane | 120 | −848.4122 | −5.8 | −23086.4849 | −0.1566 | 1e-8/40 | 4.5730 | — | 3e+00 | 4e+00 | 3.6e+00 | 7.8e+00 |  |
| SA-CASSCF | free | 3d4s / D | (4,6) | Ar·lane | 120 | −848.4166 | −10.1 | −23086.6025 | −0.2741 | 1e-9/41 | 12.0029 | 3e-03 | 1e+00 | 1e+00 | 1.2e+00 | 6.0e+00 |  |
| SA-CASSCF | pinned | 3d4s / D | (4,6) | Ar·lane | 120 | −848.4094 | −2.9 | −23086.4067 | −0.0784 | 0/43 | 12.0000 | 3e-07 | 1e-06 | 2e-06 | 3.1e-06 | 2.4e-05 | SCF core; offset = bookkeeping form |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d4d / D | (4,10) | Ar·lane | 1,200 | −848.2156 | +190.9 | −23081.1339 | +5.1945 | — | 6.0–30.0 | 5e-02 | — | — | 3.4e-13 | 2.1e-12 | SCF core; offset = bookkeeping form |
| CASCI | free | 3d4d / D | (4,10) | Ar·lane | 1,200 | −848.2371 | +169.4 | −23081.7186 | +4.6097 | — | 6.0000 | 5e-02 | 7e-15 | 1e-14 | 3.3e-13 | 1.9e-12 |  |
| CASSCF | free | 3d4d / D | (4,10) | Ar·lane | 1,200 | −848.4472 | −40.7 | −23087.4364 | −1.1081 | 1e-8/89 | 0.0592 | — | 5e+00 | 7e+00 | 2.8e+00 | 7.5e+00 |  |
| SA-CASSCF | free | 3d4d / D | (4,10) | Ar·lane | 1,200 | −848.4382 | −31.7 | −23087.1918 | −0.8635 | 1e-8/66 | 12.0000 | 1e-06 | 6e-03 | 5e-03 | 1.0e-02 | 5.2e-02 |  |
| SA-CASSCF | pinned | 3d4d / D | (4,10) | Ar·lane | 1,200 | −848.2807 | +125.8 | −23082.9042 | +3.4241 | 1e-8/115 | 6.0–29.9 | 4e-02 | 2e-01 | 3e-01 | 1.2e+00 | 9.1e+00 | SCF core; offset = bookkeeping form |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d4s4p4d / D | (4,14) | Ar·lane | 5,096 | −848.4374 | −31.0 | −23087.1708 | −0.8424 | — | 12.0000 | 2e-13 | — | — | 6.2e-13 | 3.4e-12 | SCF core; offset = bookkeeping form |
| CASCI | free | 3d4s4p4d / D | (4,14) | Ar·lane | 5,096 | −848.4374 | −31.0 | −23087.1708 | −0.8424 | — | 12.0000 | 5e-02 | 2e-13 | 4e-13 | 8.6e-13 | 4.0e-12 |  |
| CASSCF | free | 3d4s4p4d / D | (4,14) | Ar·lane | 5,096 | −848.4512 | −44.7 | −23087.5451 | −1.2167 | 1e-8/47 | 5.8927 | — | 5e+00 | 8e+00 | 3.4e+01 | 1.6e+02 |  |
| SA-CASSCF | free | 3d4s4p4d / D | (4,14) | Ar·lane | 5,096 | −848.4424 | −35.9 | −23087.3048 | −0.9765 | 1e-11/10 | 12.0000 | 3e-07 | 2e-03 | 2e-03 | 1.9e-02 | 9.6e-02 |  |
| SA-CASSCF | pinned | 3d4s4p4d / D | (4,14) | Ar·lane | 5,096 | −848.4424 | −35.9 | −23087.3048 | −0.9765 | 1e-13/10 | 12.0000 | 1e-07 | 2e-04 | 1e-03 | 6.7e-04 | 1.9e-02 | SCF core; offset = bookkeeping form |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d3p / D | (8,8) | semi·lane | 3,136 | −848.4520 | −45.5 | −23087.5662 | −1.2379 | — | 12.0000 | 3e-14 | — | — | 5.2e-13 | 4.1e-12 | SCF core; offset = bookkeeping form |
| CASCI | free | 3d3p / D | (8,8) | semi·lane | 3,136 | −848.4520 | −45.5 | −23087.5662 | −1.2379 | — | 12.0000 | 5e-02 | 9e-14 | 1e-13 | 5.7e-13 | 4.2e-12 |  |
| CASSCF | free | 3d3p / D | (8,8) | semi·lane | 3,136 | −848.4941 | −87.6 | −23088.7115 | −2.3831 | 1e-8/14 | 10.4529 | — | 1e+00 | 2e+00 | 1.4e+01 | 9.3e+01 |  |
| SA-CASSCF | free | 3d3p / D | (8,8) | semi·lane | 3,136 | −848.4756 | −69.1 | −23088.2097 | −1.8814 | 1e-9/16 | 12.0000 | 2e-07 | 7e-06 | 1e-05 | 1.6e-05 | 1.9e-04 |  |
| SA-CASSCF | pinned | 3d3p / D | (8,8) | semi·lane | 3,136 | −848.4756 | −69.1 | −23088.2097 | −1.8814 | 1e-8/16 | 12.0000 | 4e-07 | 2e-05 | 4e-05 | 6.7e-05 | 6.2e-04 | SCF core; offset = bookkeeping form |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d3p4d / D | (8,13) | semi·lane | 368,082 | −848.4991 | −92.6 | −23088.8489 | −2.5206 | — | 12.0000 | 1e-13 | — | — | 2.8e-14ψ | 2.3e-13ψ | SCF core; offset = bookkeeping form |
| CASCI | free | 3d3p4d / D | (8,13) | semi·lane | 368,082 | −848.4991 | −92.6 | −23088.8489 | −2.5206 | — | 12.0000 | 5e-02 | 9e-14 | 1e-13 | 4.0e-15ψ | 7.4e-14ψ |  |
| CASSCF | free | 3d3p4d / D | (8,13) | semi·lane | 368,082 | −848.5300 | −123.5 | −23089.6888 | −3.3605 | 1e-8/35 | 10.1972 | — | 5e+00 | 7e+00 | 3.5e-01ψ | 2.9e+00ψ |  |
| SA-CASSCF | free | 3d3p4d / D | (8,13) | semi·lane | 368,082 | −848.5108 | −104.3 | −23089.1658 | −2.8375 | 1e-8/29 | 12.0000 | 2e-08 | 2e-04 | 3e-04 | 1.2e-05ψ | 8.7e-05ψ |  |
| SA-CASSCF | pinned | 3d3p4d / D | (8,13) | semi·lane | 368,082 | −848.5108 | −104.3 | −23089.1658 | −2.8375 | 1e-8/48 | 12.0000 | 1e-08 | 4e-06 | 2e-05 | 1.2e-07ψ | 1.4e-06ψ | SCF core; offset = bookkeeping form |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3s3p3d / D | (10,9) | {1s,2s,2p,4s}·lane | 10,584 | −848.4549 | −48.4 | −23087.6452 | −1.3169 | — | 12.0000 | 6e-14 | — | — | 7.2e-14ψ | 3.5e-13ψ | SCF core; offset = bookkeeping form |
| CASCI | free | 3s3p3d / D | (10,9) | {1s,2s,2p,4s}·lane | 10,584 | −848.4549 | −48.4 | −23087.6452 | −1.3169 | — | 12.0000 | 5e-02 | 9e-14 | 1e-13 | 2.9e-14ψ | 1.2e-13ψ |  |
| CASSCF | free | 3s3p3d / D | (10,9) | {1s,2s,2p,4s}·lane | 10,584 | −848.4993 | −92.8 | −23088.8535 | −2.5251 | 1e-10/23 | 10.4149 | — | 1e+00 | 2e+00 | 2.3e-01ψ | 3.8e-01ψ |  |
| SA-CASSCF | free | 3s3p3d / D | (10,9) | {1s,2s,2p,4s}·lane | 10,584 | −848.4793 | −72.8 | −23088.3104 | −1.9820 | 1e-8/34 | 12.0000 | 4e-06 | 2e-04 | 2e-04 | 8.9e-06ψ | 6.4e-05ψ |  |
| SA-CASSCF | pinned | 3s3p3d / D | (10,9) | {1s,2s,2p,4s}·lane | 10,584 | −848.4793 | −72.8 | −23088.3098 | −1.9815 | 1e-8/52 | 12.0000 | 5e-08 | 4e-06 | 9e-06 | 4.1e-07ψ | 3.3e-06ψ | SCF core; offset = bookkeeping form |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d3p4s / D | (10,9) | {1s,2s,2p,3s}·lane | 10,584 | −848.4529 | −46.4 | −23087.5918 | −1.2635 | — | 12.0000 | 3e-13 | — | — | 2.9e-14ψ | 1.8e-13ψ | SCF core; offset = bookkeeping form |
| CASCI | free | 3d3p4s / D | (10,9) | {1s,2s,2p,3s}·lane | 10,584 | −848.4529 | −46.4 | −23087.5918 | −1.2635 | — | 12.0000 | 5e-02 | 9e-14 | 1e-13 | 1.7e-15ψ | 3.6e-14ψ |  |
| CASSCF | free | 3d3p4s / D | (10,9) | {1s,2s,2p,3s}·lane | 10,584 | −848.4993 | −92.8 | −23088.8535 | −2.5251 | 1e-8/53 | 10.4149 | — | 1e+00 | 2e+00 | 5.7e-01ψ | 3.1e-01ψ |  |
| SA-CASSCF | free | 3d3p4s / D | (10,9) | {1s,2s,2p,3s}·lane | 10,584 | −848.4793 | −72.8 | −23088.3104 | −1.9821 | 1e-8/39 | 12.0000 | 1e-05 | 4e-04 | 6e-04 | 6.5e-05ψ | 2.5e-04ψ |  |
| SA-CASSCF | pinned | 3d3p4s / D | (10,9) | {1s,2s,2p,3s}·lane | 10,584 | −848.4793 | −72.8 | −23088.3103 | −1.9820 | 1e-8/40 | 12.0000 | 3e-06 | 2e-05 | 3e-05 | 4.1e-06ψ | 1.8e-05ψ | SCF core; offset = bookkeeping form |

Declared, not yet run -- the optimizer cells are void until their CASCI anchors are read (the anchor-first scenario). One regen away: `regen_cells.py Ti --window 3d4d`.

Declared, not yet run -- the optimizer cells are void until their CASCI anchors are read (the anchor-first scenario). One regen away: `regen_cells.py Ti --window 3d4s4p4d`.
<!-- methods:end -->
