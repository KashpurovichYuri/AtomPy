# Scunc (3d¹ 4s², the uncontracted cc-pVDZ split)

The primitive split of Sc/cc-pVDZ (data/sc_ccpvdz_unc.txt — every
contraction a single-primitive block); the shell layer declares each
chemical shell's band realization (the composition certificate:
logs/band_composition.log). Validation atom of the shell-layer wave:
the multi-column fold's machinery record, not a campaign table.
Sector M_S = ½. Certified sibling: `elements/Sc/` (contracted).

<!-- ladder:begin -->
| method | route | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −759.7374 | 0 | −20673.5179 | +0.0000 | — | — | — | — | — | — | — | forced filling |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (3,43) | Ar·raw | 38,829 | −759.7479 | −10.5 | −20673.8041 | −0.2862 | — | 5.9990 | 8e-02 | 3e-02 | 4e-02 | 1.2e-03ψ | 9.7e-03ψ |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (3,43) | Ar | 38,829 | −759.7479 | −10.5 | −20673.8041 | −0.2862 | — | 5.9990 | 8e-02 | 3e-02 | 4e-02 | 1.2e-03ψ | 9.7e-03ψ |  |
| CASCI | free | 3d4s / MO | (3,43) | Ar | 38,829 | −759.7494 | −12.0 | −20673.8457 | −0.3279 | — | 5.9989 | 5e-02 | 3e-02 | 5e-02 | 1.2e-03ψ | 1.0e-02ψ |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (3,43) | Ar·pur | 38,829 | −759.7467 | −9.3 | −20673.7703 | −0.2524 | — | 6.0000 | 1e-01 | 8e-14 | 2e-13 | 3.6e-13ψ | 1.1e-12ψ |  |
| CASCI | free | 3d4s / MO | (3,43) | Ar·pur | 38,829 | −759.7482 | −10.8 | −20673.8115 | −0.2936 | — | 5.9999 | 4e-02 | 2e-02 | 3e-02 | 4.4e-04ψ | 2.4e-03ψ |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d4s / AO | (3,43) | Ar·lane | 38,829 | −759.7177 | +19.7 | −20672.9811 | +0.5367 | — | 6.0000 | 4e-12 | 8e-14 | 2e-13 | 3.2e-13ψ | 1.3e-12ψ |  |
<!-- ladder:end -->
