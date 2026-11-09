# Sc-unc (3d¹ 4s², the uncontracted cc-pVDZ split)

Primitive split of Sc/cc-pVDZ (data/sc_ccpvdz_unc.txt). Windows:
DISTILLED, k=1 per shell (`window_distilled` — the reference's own
radial × the m-set; CAS(3,5)/CAS(3,6), construction-uniform with
`elements/Sc/`, whose verification legs reproduced the certified
records at 0.00e+00 mEh). The `mo`/`moc` rows are CONTROLS —
m-incomplete by construction (the SCF's hybrid-open gauge). Reference:
the trusted ROHF (−759.73738486, rungs/rohf.chk). The pinned rows
freeze the reference's OWN core (the span-restricted SCF columns);
their residual vs the anchor (+5.1 mEh) is the projected-density
bookkeeping form. Sector M_S = ½. Sibling: `elements/Sc/`
(contracted).

<!-- ladder:begin -->
| method | route | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −759.7374 | 0 | −20673.5179 | +0.0000 | — | — | — | — | — | — | — | forced filling |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d / AO | (3,5) | Ar | 50 | −759.3360 | +401.3 | −20662.5968 | +10.9210 | — | 11.9994 | 5e-02 | 3e-02 | 4e-02 | 1.2e-02 | 7.2e-02 |  |
| CASCI | free | 3d / MO | (3,5) | Ar | 50 | −759.7466 | −9.2 | −20673.7694 | −0.2515 | — | 0.0057 | 1e-01 | 2e+00 | 3e+00 | 8.4e-02 | 2.2e-01 | control: m-incomplete by construction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (3,6) | Ar·raw | 90 | −759.7342 | +3.2 | −20673.4319 | +0.0860 | — | 5.9990 | 9e-02 | 3e-02 | 4e-02 | 4.6e-02 | 2.5e-01 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (3,6) | Ar | 90 | −759.7342 | +3.2 | −20673.4319 | +0.0860 | — | 5.9989 | 9e-02 | 3e-02 | 4e-02 | 1.4e-02 | 8.4e-02 |  |
| CASCI | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7395 | −2.1 | −20673.5753 | −0.0574 | — | 1.4464 | 1e-01 | 2e+00 | 4e+00 | 1.5e+00 | 6.5e+00 | control: m-incomplete by construction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / AO | (3,6) | Ar·pur | 90 | −759.7327 | +4.7 | −20673.3903 | +0.1276 | — | 6.0000 | 9e-02 | 2e-14 | 2e-14 | 6.3e-15 | 9.3e-14 |  |
| CASCI | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7380 | −0.6 | −20673.5355 | −0.0176 | — | 1.4467 | 1e-01 | 2e+00 | 4e+00 | 1.5e+00 | 6.5e+00 | control: m-incomplete by construction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d4s / AO | (3,6) | Ar·lane | 90 | −759.7322 | +5.1 | −20673.3781 | +0.1398 | — | 6.0000 | 9e-16 | 1e-14 | 2e-14 | 5.6e-15 | 4.3e-14 | SCF core; offset = bookkeeping form |
<!-- ladder:end -->
