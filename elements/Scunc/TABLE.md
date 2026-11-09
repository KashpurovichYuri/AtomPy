# Scunc (3d¹ 4s², the uncontracted cc-pVDZ split)

The primitive split of Sc/cc-pVDZ (data/sc_ccpvdz_unc.txt — every
contraction a single-primitive block); the shell layer declares each
chemical shell's band realization (the composition certificate:
logs/band_composition.log). The windows are the COMPRESSED-SHELL
form (k=1 reference radial × the m-set per shell, `shell_radials`):
the CAS labels mirror the certified contracted record's shapes —
CAS(3,5) [3d], CAS(3,6) [3d4s] — with the capture certificates in
logs/regen_compressed.log (the 3d family holds its MO to 0.999; the
4s family's single radial to 0.81 — the node-tail honesty of a raw
family vs a baked-in contraction). The folds stay FULL-family (the
mini-Roothaan grammar). Validation atom of the shell-layer wave.
Sector M_S = ½. Certified sibling: `elements/Sc/` (contracted).

**Field disclosure.** Every row anchors on the TRUE neutral ROHF (the
ROHF row — the tables are NOT re-anchored). The folded-reference round
adds a second reference class: the confined-core field with the
valence re-relaxed (`run_folded`, logs/folded_probe.log —
−759.58645602, +150.9 mEh above the anchor BY ZERO: the frozen
confined core is the free ion on a wide span, plus the variational
cost). On that field the lane closes the lane-vs-free gap to
0.085 mEh — the residual is the C1/D2h gauge of the window cut; the
free-route cut declares its dominance gate at 0.4 (the C1 gauge's 4s
signal 0.4628, D2h's gauge passes 0.5; measured honest, ⟨L²⟩ = 6
exact). Full record: logs/matrix_folded.log (the ²D lane total
−759.689259, +48.1 mEh vs the anchor); folded-field cells carry the
`folded` tag in their addresses (the ninth address slot). The
ladder's pinned row is the HISTORICAL neutral-field lane record — the
+38.7 mEh gap that motivated the fold — kept for provenance.

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
