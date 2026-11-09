# Scunc (3d¹ 4s², the uncontracted cc-pVDZ split)

The primitive split of Sc/cc-pVDZ (data/sc_ccpvdz_unc.txt — every
contraction a single-primitive block); the shell layer declares each
chemical shell's band realization (the composition certificate:
logs/band_composition.log). The windows are the DISTILLED form
(window_distilled: k=1 radial per shell, distilled from the
reference's occupied MOs by the member-summed radial Gram,
replicated on the m-set, `shell_radials`): the CAS labels mirror the
certified contracted record's shapes — CAS(3,5) [3d], CAS(3,6)
[3d4s] — and the whole ladder is construction-uniform with
`elements/Sc/` (the contracted sibling's windows are the 1:1
identity of the same construction — its verification legs
reproduced the certified records at 0.00e+00 mEh). The captures
(3d 0.999; 4s 0.81 — the node-tail honesty of a raw family vs a
baked-in contraction) are the window-shape readout. The folds stay
FULL-family (the mini-Roothaan grammar). Validation atom of the
shell-layer wave. Sector M_S = ½. Certified sibling:
`elements/Sc/` (contracted).

**Field disclosure.** Every row anchors on the TRUE neutral ROHF (the
ROHF row — the tables are NOT re-anchored). The lane's pinned row
runs on the neutral field with the confined core on the FULL-family
span — its offset vs the free rows (+38.7 mEh) is STRUCTURAL (the
confined core ≠ the SCF core on a wide span; measured basis-stable
on the contracted control too, elements/Ti/logs/folded_review.log)
— the lane's value is the independent machine-exact term
certification (²D: L² = 6 ×5, ptp 3e-15). The folded-reference
round (logs/matrix_folded.log, draw-documented) proved the field
fix: on the folded field the lane closes to the free route at
0.085 mEh (the dominance declaration 0.4 was that round's cut);
folded-field cells carry the `folded` tag in their addresses and
sit BY ZERO above the anchors (+150.9 mEh reference, +48.1 lane).

<!-- ladder:begin -->
| method | route | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −759.7374 | 0 | −20673.5179 | +0.0000 | — | — | — | — | — | — | — | forced filling |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d / AO | (3,5) | Ar | 50 | −759.3360 | +401.3 | −20662.5968 | +10.9210 | — | 11.9994 | 5e-02 | 3e-02 | 4e-02 | 1.2e-02 | 7.2e-02 |  |
| CASCI | free | 3d / MO | (3,5) | Ar | 50 | −759.7466 | −9.2 | −20673.7694 | −0.2515 | — | 0.0057 | 1e-01 | 2e+00 | 3e+00 | 8.4e-02 | 2.2e-01 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (3,6) | Ar·raw | 90 | −759.7342 | +3.2 | −20673.4319 | +0.0860 | — | 5.9990 | 9e-02 | 3e-02 | 4e-02 | 4.6e-02 | 2.5e-01 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (3,6) | Ar | 90 | −759.7342 | +3.2 | −20673.4319 | +0.0860 | — | 5.9989 | 9e-02 | 3e-02 | 4e-02 | 1.4e-02 | 8.4e-02 |  |
| CASCI | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7395 | −2.1 | −20673.5753 | −0.0574 | — | 1.4464 | 1e-01 | 2e+00 | 4e+00 | 1.5e+00 | 6.5e+00 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (3,6) | Ar·pur | 90 | −759.7327 | +4.7 | −20673.3903 | +0.1276 | — | 6.0000 | 9e-02 | 2e-14 | 2e-14 | 6.3e-15 | 9.3e-14 |  |
| CASCI | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7380 | −0.6 | −20673.5355 | −0.0176 | — | 1.4467 | 1e-01 | 2e+00 | 4e+00 | 1.5e+00 | 6.5e+00 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d4s / AO | (3,6) | Ar·lane | 90 | −759.6987 | +38.7 | −20672.4652 | +1.0526 | — | 6.0000 | 3e-15 | 1e-14 | 2e-14 | 2.8e-14 | 1.0e-13 |  |
<!-- ladder:end -->
