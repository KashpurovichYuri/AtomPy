# Ti-unc (3d² 4s², the uncontracted cc-pVDZ split)

Primitive split of Ti/cc-pVDZ (data/ti_ccpvdz_unc.txt). Windows:
DISTILLED, k=1 per shell (`window_distilled` — the reference's own
radial × the m-set; CAS(2,5)/CAS(4,6), construction-uniform with
`elements/Ti/`, whose verification legs reproduced the certified
records at 0.00e+00 mEh). The kind block compares the four window
sources at both ruler folds: `MO`/`MOC` are CONTROLS (m-incomplete
by construction — the SCF's hybrid-open gauge; MO is l-impure here,
win_l 0.60-class), `NO` is the density leg, degenerate with AO
(+3.3 µEh, both folds). The optimizer class and the frac fold (no
measured fractional-occupation basin on this basis) are out of scope
this wave (void at the key, one `--full` away). Reference: the trusted
deep-basin ROHF (−848.40744343, logs/census.log). The pinned row and
its free twin (route agreement at 1e-9 Eh) freeze the reference's
OWN core; their residual vs the anchor (+6.1 mEh) is the
projected-density bookkeeping form. Sector M_S = 1. Sibling:
`elements/Ti/` (contracted).

<!-- ladder:begin -->
| method | route | window | CAS | fold | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −848.4074 | 0 | −23086.3543 | +0.0000 | — | — | — | — | — | — | — | forced filling |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d / D | (2,5) | Ar | 10 | −848.4074 | −0.0 | −23086.3543 | −0.0000 | — | 12.0000 | 6e-02 | 1e-08 | 9e-09 | 4.5e-15 | 3.7e-14 |  |
| CASCI | free | 3d / MO | (2,5) | Ar | 10 | −848.4074 | −0.0 | −23086.3543 | −0.0000 | — | 0.0000 | 9e-01 | 3e+00 | 4e+00 | 5.7e-03 | 1.0e-09 | control: m-incomplete by construction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (4,6) | Ar·raw | 120 | −848.4037 | +3.8 | −23086.2515 | +0.1028 | — | 12.0000 | 6e-02 | 4e-09 | 4e-09 | 6.5e-02 | 3.4e-01 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (4,6) | Ar | 120 | −848.4037 | +3.8 | −23086.2515 | +0.1028 | — | 12.0000 | 6e-02 | 4e-09 | 4e-09 | 2.6e-09 | 1.3e-08 |  |
| CASCI | free | 3d4s / MO | (4,6) | Ar | 120 | −848.4172 | −9.8 | −23086.6210 | −0.2667 | — | 0.0048 | 1e-01 | 3e+00 | 4e+00 | 1.5e-01 | 3.8e-01 | control: m-incomplete by construction |
| CASCI | free | 3d4s / MOC | (4,6) | Ar | 120 | −848.4104 | −2.9 | −23086.4335 | −0.0792 | — | 4.9041 | 1e-01 | 3e+00 | 4e+00 | 3.7e+00 | 1.4e+01 | control: m-incomplete by construction |
| CASCI | free | 3d4sno2 / NO | (4,6) | Ar | 120 | −848.4037 | +3.8 | −23086.2515 | +0.1029 | — | 12.0000 | 6e-02 | 4e-09 | 4e-09 | 2.8e-09 | 1.4e-08 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (4,6) | Ar·pur | 120 | −848.4018 | +5.6 | −23086.2018 | +0.1525 | — | 12.0000 | 6e-02 | 1e-14 | 2e-14 | 1.1e-14 | 2.6e-13 |  |
| CASCI | free | 3d4s / MO | (4,6) | Ar·pur | 120 | −848.4154 | −8.0 | −23086.5712 | −0.2169 | — | 0.0048 | 1e-01 | 3e+00 | 4e+00 | 1.5e-01 | 3.8e-01 | control: m-incomplete by construction |
| CASCI | free | 3d4s / MOC | (4,6) | Ar·pur | 120 | −848.4085 | −1.1 | −23086.3839 | −0.0296 | — | 4.9041 | 1e-01 | 3e+00 | 4e+00 | 3.7e+00 | 1.4e+01 | control: m-incomplete by construction |
| CASCI | free | 3d4sno2 / NO | (4,6) | Ar·pur | 120 | −848.4018 | +5.6 | −23086.2017 | +0.1526 | — | 12.0000 | 6e-02 | 1e-14 | 2e-14 | 2.2e-14 | 2.7e-13 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d4s / D | (4,6) | Ar·lane | 120 | −848.4014 | +6.1 | −23086.1894 | +0.1650 | — | 12.0000 | 3e-15 | 1e-14 | 2e-14 | 2.4e-14 | 1.5e-13 | SCF core; offset = bookkeeping form |
| CASCI | free | 3d4s / D | (4,6) | Ar·lane | 120 | −848.4014 | +6.1 | −23086.1894 | +0.1650 | — | 12.0000 | 6e-02 | 1e-14 | 2e-14 | 2.2e-14 | 2.7e-13 |  |
<!-- ladder:end -->
