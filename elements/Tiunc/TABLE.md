# Tiunc (3d² 4s², the uncontracted cc-pVDZ split)

The primitive split of Ti/cc-pVDZ (data/ti_ccpvdz_unc.txt — every
contraction a single-primitive block); the shell layer declares each
chemical shell's band realization (the composition certificate:
logs/band_composition.log). The windows are the DISTILLED form
(window_distilled: k=1 radial per shell, distilled from the
reference's occupied MOs by the member-summed radial Gram,
replicated on the m-set, `shell_radials`): the CAS labels mirror the
certified contracted record's shapes — CAS(2,5) [3d at the free-ion
sector], CAS(4,6) [3d4s] — and the whole ladder is construction-
uniform with `elements/Ti/` (the contracted sibling's windows are
the 1:1 identity of the same construction — its verification legs
reproduced the certified records at 0.00e+00 mEh). The d radial is
a three-band 6d/7d/8d mixture (logs/distilled_probe.log — no
label-picked band rotation approaches it; the 2-band scans sit
+753 mEh above). The folds stay FULL-family. The reference is the
SEEDED deep-basin ROHF (−848.40744343, logs/census.log — the plain
aufbau lands 11 mEh high on this basis). Validation atom of the
shell-layer wave. Sector M_S = 1. Certified sibling: `elements/Ti/`
(contracted).

**Field disclosure.** Every row anchors on the TRUE neutral ROHF (the
ROHF row — the tables are NOT re-anchored). The lane's pinned row
runs on the neutral field with the confined core on the FULL-family
span — its offset vs the free rows (+124.3 mEh) is STRUCTURAL (the
confined core ≠ the SCF core on a wide span; measured
basis-stable on the contracted control too, logs/folded_review.log)
— the lane's value is the independent machine-exact term
certification (³F: L² = 12 ×7, ptp 5e-15). The folded-reference
round (logs/matrix_folded.log, draw-documented) proved the field
fix: on the folded field the lane closes to the free route at
5.5e-13 Eh; folded-field cells carry the `folded` tag in their
addresses and sit BY ZERO above the anchors (+399.7 mEh reference,
+135.3 lane).

<!-- ladder:begin -->
| method | route | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −848.4074 | 0 | −23086.3543 | +0.0000 | — | — | — | — | — | — | — | forced filling |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d / AO | (2,5) | Ar | 10 | −848.4074 | −0.0 | −23086.3543 | −0.0000 | — | 12.0000 | 6e-02 | 1e-08 | 9e-09 | 4.5e-15 | 3.7e-14 |  |
| CASCI | free | 3d / MO | (2,5) | Ar | 10 | −848.4074 | −0.0 | −23086.3543 | −0.0000 | — | 0.0000 | 9e-01 | 3e+00 | 4e+00 | 5.7e-03 | 1.0e-09 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (4,6) | Ar·raw | 120 | −848.4037 | +3.8 | −23086.2515 | +0.1028 | — | 12.0000 | 6e-02 | 4e-09 | 4e-09 | 6.5e-02 | 3.4e-01 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (4,6) | Ar | 120 | −848.4037 | +3.8 | −23086.2515 | +0.1028 | — | 12.0000 | 6e-02 | 4e-09 | 4e-09 | 2.6e-09 | 1.3e-08 |  |
| CASCI | free | 3d4s / MO | (4,6) | Ar | 120 | −848.4104 | −2.9 | −23086.4335 | −0.0792 | — | 4.9041 | 1e-01 | 3e+00 | 4e+00 | 3.7e+00 | 1.4e+01 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (4,6) | Ar·pur | 120 | −848.4018 | +5.6 | −23086.2018 | +0.1525 | — | 12.0000 | 6e-02 | 1e-14 | 2e-14 | 1.1e-14 | 2.6e-13 |  |
| CASCI | free | 3d4s / MO | (4,6) | Ar·pur | 120 | −848.4085 | −1.1 | −23086.3839 | −0.0296 | — | 4.9041 | 1e-01 | 3e+00 | 4e+00 | 3.7e+00 | 1.4e+01 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d4s / AO | (4,6) | Ar·lane | 120 | −848.2832 | +124.3 | −23082.9732 | +3.3811 | — | 12.0000 | 5e-15 | 1e-14 | 2e-14 | 3.9e-14 | 2.1e-13 |  |
<!-- ladder:end -->
