# Tiunc (3d² 4s², the uncontracted cc-pVDZ split)

The primitive split of Ti/cc-pVDZ (data/ti_ccpvdz_unc.txt — every
contraction a single-primitive block); the shell layer declares each
chemical shell's band realization (the composition certificate:
logs/band_composition.log), the 3d as the dominant-5 bands (the full
family is 530k determinants at the (3,1) sector — hours per cell).
Validation atom of the shell-layer wave: the multi-column fold's
machinery record, not a campaign table. Sector M_S = 1.
Certified sibling: `elements/Ti/` (contracted).

<!-- ladder:begin -->
| method | route | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −848.3964 | 0 | −23086.0535 | +0.0000 | — | — | — | — | — | — | — | forced filling |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (4,28) | Ar·raw | 91,728 | −848.3709 | +25.5 | −23085.3600 | +0.6935 | — | 11.9904 | 1e-01 | 2e-02 | 2e-02 | 1.2e-03ψ | 7.3e-03ψ |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (4,28) | Ar | 91,728 | −848.3709 | +25.5 | −23085.3600 | +0.6935 | — | 11.9904 | 1e-01 | 2e-02 | 2e-02 | 1.2e-03ψ | 7.3e-03ψ |  |
| CASCI | free | 3d4s / MO | (4,28) | Ar | 91,728 | −848.4073 | −10.9 | −23086.3512 | −0.2976 | — | 4.0898 | 2e-01 | 3e+00 | 4e+00 | 3.6e-01ψ | 1.1e+00ψ |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (4,28) | Ar·pur | 91,728 | −848.3721 | +24.3 | −23085.3933 | +0.6602 | — | 12.0169 | 2e-01 | 2e-13 | 2e-13 | 2.6e-14ψ | 2.0e-13ψ |  |
| CASCI | free | 3d4s / MO | (4,28) | Ar·pur | 91,728 | −848.4053 | −8.9 | −23086.2950 | −0.2415 | — | 4.0885 | 2e-01 | 3e+00 | 4e+00 | 3.6e-01ψ | 1.1e+00ψ |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d4s / AO | (4,28) | Ar·lane | 91,728 | −848.2786 | +117.7 | −23082.8496 | +3.2039 | — | 12.0000 | 1e-13 | 2e-13 | 2e-13 | 4.5e-14ψ | 1.9e-13ψ |  |
<!-- ladder:end -->
