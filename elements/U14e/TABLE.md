# U14e (f³d¹ outside the ECP core)

Stuttgart RLC ECP78 + cut-Dyall valence basis; CAMPAIGN key `U14e`,
symbol `U` (the ECP removes 78 core electrons → 14). Reference ROHF
−50.381852 Eh (gated Newton; rechecked from scratch 2026-11-06: same
basin, |dE| = 4e-09; the open-shell gauge rotates per draw — the
determinant ⟨L²⟩ = 40.85, the span certificates broken;
`logs/rohf_recheck.log`). Sector M_S = 2. The pre-campaign record:
`elements/U/` (frozen; `logs/u_run*.log`, `logs/u_sa*.log`,
`logs/u_six13.log`).

Window labels are exponent-grid coordinates, not principal numbers:
the physical 5f sits on bands 9f+10f+11f (10f dominant), the 6d on
12d+13d, the core 6s on 14s+15s, the 6p trio on 13p+14p+15p (14p),
the 7s on 18s+19s — "10f13d" IS the simplest 5f6d window this basis
expresses. The window column's slash names the START: AO = the raw
shell bands (exactly closed), MO-char = the character-selected
canonical columns (the SCF's own valence + the virtual complement
ranked by the named bands' S-projection), MO = the energy-ordered
canonical slice (the control: the reference echo). The optimizer rows
carry their start in the same column; the optimized spans stay
on-class in l-space (the SA exactly; the free CASSCF leaks ~0.1
columns; the 7s MO-start draws wandered 1.4 while cap-limited —
`logs/window_audit_ao.log`, `logs/window_audit_moc.log`).

<!-- ladder:begin -->
| method | route | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −50.3819 | 0 | −1370.9607 | +0.0000 | — | — | — | — | — | — | — | forced filling |
| CASCI | free | 9f10f11f12d13d / AO | (4,31) | Ne·raw | 31,465 | −50.2115 | +170.3 | −1366.3261 | +4.6347 | — | 33.9840 | 8e-03 | 7e-02 | 4e-01 | 3.1e-03ψ | 5.6e-02ψ |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 9f10f11f12d13d / AO | (4,31) | Ne | 31,465 | −50.2115 | +170.3 | −1366.3261 | +4.6347 | — | 33.9840 | 8e-03 | 7e-02 | 4e-01 | 3.1e-03ψ | 5.6e-02ψ |  |
| CASCI | free | 9f10f11f12d13d / MO-char | (4,31) | Ne | 31,465 | −50.3831 | −1.2 | −1370.9944 | −0.0337 | — | 0.0187 | 1e+00 | 7e+00 | 9e+00 | 1.7e-01ψ | 1.1e+00ψ |  |
| SA-CASSCF | free | 9f10f11f12d13d / AO | (4,31) | Ne | 31,465 | −50.3864 | −4.5 | −1371.0839 | −0.1231 | 1e-6/200 | 71.9999 | 6e-02 | 3e-01 | 4e-01 | 2.2e-03ψ | 2.9e-02ψ |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 9f10f11f12d13d / AO | (4,31) | Ne·pur | 31,465 | −50.2065 | +175.4 | −1366.1885 | +4.7722 | — | 29.5769 | 4e-03 | 1e-12 | 2e-12 | 3.0e-14ψ | 2.6e-13ψ |  |
| CASCI | free | 9f10f11f12d13d / MO-char | (4,31) | Ne·pur | 31,465 | −50.3793 | +2.6 | −1370.8912 | +0.0695 | — | 0.0189 | 1e+00 | 7e+00 | 9e+00 | 1.7e-01ψ | 1.1e+00ψ |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 10f13d / AO | (4,12) | Ne·raw | 495 | −49.4588 | +923.1 | −1345.8428 | +25.1179 | — | 38.9451 | 2e-03 | 7e-02 | 4e-01 | 8.1e-02 | 2.6e+00 |  |
| CASCI | free | 10f13d / AO | (4,12) | Ne | 495 | −49.4588 | +923.1 | −1345.8428 | +25.1179 | — | 41.9649 | 2e-03 | 7e-02 | 4e-01 | 3.0e-02 | 6.0e-01 |  |
| CASCI | free | 10f13d / AO | (4,12) | Ne·pur | 495 | −49.4587 | +923.1 | −1345.8412 | +25.1195 | — | 42.0000 | 1e-03 | 2e-13 | 3e-13 | 9.7e-14 | 1.8e-12 |  |
| CASCI | free | 10f13d / MO-char | (4,12) | Ne | 495 | −50.3821 | −0.3 | −1370.9680 | −0.0073 | — | 0.0159 | 1e+00 | 7e+00 | 9e+00 | 8.2e-01 | 1.1e+01 |  |
| CASCI | free | 10f13d / MO-char | (4,12) | Ne·pur | 495 | −50.3783 | +3.5 | −1370.8648 | +0.0959 | — | 0.0160 | 1e+00 | 7e+00 | 9e+00 | 8.3e-01 | 1.1e+01 |  |
| CASSCF | free | 10f13d / AO | (4,12) | Ne | 495 | −50.4297 | −47.8 | −1372.2627 | −1.3019 | 1e-5/200 | 64.5991 | — | 5e+00 | 7e+00 | 1.7e+01 | 1.2e+02 |  |
| CASSCF | free | 10f13d / MO-char | (4,12) | Ne | 495 | −50.4307 | −48.8 | −1372.2894 | −1.3287 | 1e-8/176 | 64.7359 | — | 6e+00 | 6e+00 | 9.2e+00 | 5.5e+01 |  |
| CASSCF | free | 10f13d / MO-char | (4,12) | Ne·pur | 495 | −50.4276 | −45.8 | −1372.2060 | −1.2452 | 1e-5/198 | 49.8403 | — | 7e+00 | 8e+00 | 8.3e+00 | 6.0e+01 |  |
| SA-CASSCF | free | 10f13d / AO | (4,12) | Ne | 495 | −50.4238 | −42.0 | −1372.1025 | −1.1418 | 1e-8/56 | 72.0000 | 3e-06 | 6e-04 | 7e-04 | 1.6e-04 | 1.2e-03 |  |
| SA-CASSCF | free | 10f13d / MO-char | (4,12) | Ne | 495 | −50.4238 | −42.0 | −1372.1026 | −1.1419 | 1e-8/160 | 72.0000 | 2e-05 | 2e-04 | 7e-04 | 1.1e-04 | 2.2e-03 |  |
| SA-CASSCF | free | 10f13d / MO-char | (4,12) | Ne·pur | 495 | −50.4226 | −40.8 | −1372.0707 | −1.1099 | 1e-4/200 | 59.6074 | 4e-02 | 5e+00 | 5e+00 | 1.6e+01 | 2.4e+02 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 10f13d7s / AO | (6,13) | Ne·raw | 16,731 | −49.3191 | +1062.7 | −1342.0428 | +28.9180 | — | 11.9979 | 2e-02 | 5e-02 | 1e-01 | 7.3e-04ψ | 5.9e-03ψ |  |
| CASCI | free | 10f13d7s / AO | (6,13) | Ne | 16,731 | −49.3191 | +1062.7 | −1342.0428 | +28.9180 | — | 11.9979 | 2e-02 | 5e-02 | 1e-01 | 7.3e-04ψ | 5.9e-03ψ |  |
| CASCI | free | 10f13d7s / AO | (6,13) | Ne·pur | 16,731 | −49.3164 | +1065.5 | −1341.9681 | +28.9926 | — | 12.0000 | 2e-02 | 2e-13 | 3e-13 | 4.6e-15ψ | 3.6e-14ψ |  |
| CASCI | free | 10f13d7s / MO-char | (6,13) | Ne | 16,731 | −50.3820 | −0.2 | −1370.9660 | −0.0052 | — | 0.0000 | 1e+00 | 5e+00 | 9e+00 | 1.4e-03ψ | 3.6e-03ψ |  |
| CASCI | free | 10f13d7s / MO-char | (6,13) | Ne·pur | 16,731 | −50.3802 | +1.6 | −1370.9169 | +0.0439 | — | 0.0000 | 1e+00 | 5e+00 | 9e+00 | 1.4e-03ψ | 3.6e-03ψ |  |
| CASSCF | free | 10f13d7s / AO | (6,13) | Ne | 16,731 | −50.4535 | −71.7 | −1372.9116 | −1.9509 | 1e-4/200 | 82.9982 | — | 3e+00 | 6e+00 | 4.1e+00ψ | 1.3e+01ψ |  |
| CASSCF | free | 10f13d7s / MO-char | (6,13) | Ne | 16,731 | −50.4561 | −74.3 | −1372.9814 | −2.0207 | 1e-4/200 | 50.8669 | — | 5e+00 | 6e+00 | 6.4e+00ψ | 1.3e+01ψ |  |
| CASSCF | free | 10f13d7s / MO-char | (6,13) | Ne·pur | 16,731 | −50.4231 | −41.2 | −1372.0825 | −1.1218 | 1e-4/200 | 23.8576 | — | 6e+00 | 8e+00 | 1.6e+01ψ | 1.1e+02ψ |  |
| SA-CASSCF | free | 10f13d7s / AO | (6,13) | Ne | 16,731 | −50.4529 | −71.0 | −1372.8927 | −1.9320 | 1e-9/190 | 89.9987 | 4e-04 | 2e-02 | 2e-02 | 3.9e-04ψ | 6.8e-03ψ |  |
| SA-CASSCF | free | 10f13d7s / MO-char | (6,13) | Ne | 16,731 | −50.4089 | −27.0 | −1371.6966 | −0.7358 | 1e-4/200 | 24.8171 | 8e-02 | 4e+00 | 9e+00 | 7.1e+00ψ | 6.9e+01ψ |  |
| SA-CASSCF | free | 10f13d7s / MO-char | (6,13) | Ne·pur | 16,731 | −50.3860 | −4.1 | −1371.0726 | −0.1119 | 1e-4/198 | 1.3487 | 1e-01 | 6e+00 | 8e+00 | 2.6e+01ψ | 5.6e+01ψ |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 10f13d14p / AO | (4,15) | Ne·raw | 1,365 | −49.4588 | +923.0 | −1345.8436 | +25.1172 | — | 38.9378 | 2e-03 | 1e-01 | 4e-01 | 1.5e+00 | 3.8e+01 |  |
| CASCI | free | 10f13d14p / AO | (4,15) | Ne | 1,365 | −49.4588 | +923.0 | −1345.8436 | +25.1172 | — | 41.9649 | 2e-03 | 1e-01 | 4e-01 | 1.5e+00 | 3.6e+01 |  |
| CASCI | free | 10f13d14p / AO | (4,15) | Ne·pur | 1,365 | −49.4587 | +923.1 | −1345.8416 | +25.1192 | — | 42.0000 | 1e-03 | 3e-13 | 5e-13 | 1.4e-12 | 1.8e-11 |  |
| CASSCF | free | 10f13d14p / AO | (4,15) | Ne | 1,365 | −50.4310 | −49.2 | −1372.2988 | −1.3380 | 1e-5/200 | 66.8672 | — | 6e+00 | 7e+00 | 2.8e+02 | 1.9e+03 |  |
| SA-CASSCF | free | 10f13d14p / AO | (4,15) | Ne | 1,365 | −50.4264 | −44.5 | −1372.1726 | −1.2118 | 1e-7/200 | 71.7790 | 2e-03 | 5e+00 | 7e+00 | 8.0e+00 | 6.0e+01 |  |
<!-- ladder:end -->
