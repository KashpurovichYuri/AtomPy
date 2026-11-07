# Tiunc (3d² 4s², the uncontracted cc-pVDZ split)

The primitive split of Ti/cc-pVDZ (data/ti_ccpvdz_unc.txt — every
contraction a single-primitive block); the shell layer declares each
chemical shell's band realization (the composition certificate:
logs/band_composition.log). The windows are the COMPRESSED-SHELL
form (k=1 reference radial × the m-set per shell, `shell_radials`):
the CAS labels mirror the certified contracted record's shapes —
CAS(2,5) [3d at the free-ion sector], CAS(4,6) [3d4s] — with the
capture certificates in logs/regen_compressed.log (the 3d family
holds its MOs to 1.000 — the lattice row's ΔE −0.0 next to the
contracted record's +6.4 is that fidelity read directly). The folds
stay FULL-family. The reference is the SEEDED deep-basin ROHF
(−848.40648875, logs/basin_probe.log — the plain aufbau lands 11 mEh
high on this basis). Validation atom of the shell-layer wave.
Sector M_S = 1. Certified sibling: `elements/Ti/` (contracted).

<!-- ladder:begin -->
| method | route | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −848.4065 | 0 | −23086.3283 | +0.0000 | — | — | — | — | — | — | — | forced filling |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d / AO | (2,5) | Ar | 10 | −848.4065 | −0.0 | −23086.3283 | −0.0000 | — | 12.0000 | 6e-02 | 2e-07 | 3e-07 | 1.1e-14 | 3.9e-14 |  |
| CASCI | free | 3d / MO | (2,5) | Ar | 10 | −848.4065 | −0.0 | −23086.3283 | −0.0000 | — | 0.0000 | 9e-01 | 3e+00 | 4e+00 | 5.9e-03 | 4.2e-09 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (4,6) | Ar·raw | 120 | −848.4030 | +3.5 | −23086.2325 | +0.0959 | — | 12.0000 | 6e-02 | 6e-09 | 7e-09 | 5.0e-02 | 2.6e-01 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (4,6) | Ar | 120 | −848.4030 | +3.5 | −23086.2325 | +0.0959 | — | 12.0000 | 6e-02 | 6e-09 | 7e-09 | 4.5e-09 | 2.4e-08 |  |
| CASCI | free | 3d4s / MO | (4,6) | Ar | 120 | −848.4171 | −10.6 | −23086.6176 | −0.2893 | — | 0.0051 | 1e-01 | 3e+00 | 4e+00 | 1.6e-01 | 3.9e-01 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (4,6) | Ar·pur | 120 | −848.4016 | +4.9 | −23086.1942 | +0.1341 | — | 12.0000 | 6e-02 | 1e-14 | 2e-14 | 1.9e-14 | 1.9e-13 |  |
| CASCI | free | 3d4s / MO | (4,6) | Ar·pur | 120 | −848.4157 | −9.2 | −23086.5795 | −0.2511 | — | 0.0051 | 1e-01 | 3e+00 | 4e+00 | 1.6e-01 | 3.9e-01 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d4s / AO | (4,6) | Ar·lane | 120 | −848.2832 | +123.3 | −23082.9732 | +3.3551 | — | 12.0000 | 6e-15 | 1e-14 | 2e-14 | 2.9e-14 | 2.0e-13 |  |
<!-- ladder:end -->
