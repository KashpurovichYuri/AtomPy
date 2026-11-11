# Sc-unc (3d¹ 4s², the uncontracted cc-pVDZ split)

Primitive split of Sc/cc-pVDZ (data/sc_ccpvdz_unc.txt). Windows:
DISTILLED, k=1 per shell (`window_distilled` — the reference's own
radial × the m-set; CAS(3,5)/CAS(3,6), construction-uniform with
`elements/Sc/`, whose verification legs reproduced the certified
records at 0.00e+00 mEh). The kind block compares the four window
sources at both ruler cores: `MO`/`MOC` are CONTROLS (m-incomplete
by construction — the SCF's hybrid-open gauge), the `NO` density
leg is store-only, degenerate with D at the closed substrate
(+3.2 µEh at Ar·pur -- the rendering diet drops the repeat row).
The prior axis lands `@r1`-ONLY here (the head configuration's own
density block — the same trusted ROHF the incumbent distills): on
this atom `D` and `@r1` are the SAME WINDOW (both +6.0 mEh at Ar;
+4.7-class at Ar·pur) — the DEGENERACY CERTIFICATE: no authors'
contraction exists on the primitive split, the incumbent IS the
span-Gram distillation of the head, and for the d¹ head the span and
density signals carry the same radial information (two independent
signal definitions, one window — the same-source controlled pair at
its degenerate point). The `@s` windows (the series average) are
BLOCKED on this basis: the series' non-head configurations need
fresh forced-occupation SCFs, and every fresh run diverges here (the
twins' fresh-SCF disease, measured 2026-11-11 — map Issues); the
port's @s half waits on that repair. The `3d`-span priors are
internal-only (the incumbent's own 3d window misses the 4s
correlation by ~0.4 Eh; `@r1` follows it to 1.5 mEh).
This atom's canonical core is the ANISOTROPIC class: the Ar (pol)
rows carry the assembly defect (closures ~3e-2, ⟨L²⟩ 5.9989 —
identical across kinds) and their ground-root energies are
gauge-sensitive at the split scale (the ²D block splits 0.7 mEh
under ‖[H,L²]‖ 8e-2); the CLOSED comparisons are the Ar·pur rows and
the lanes. The optimizer class and the frac core (no measured
fractional-occupation basin on this basis) are out of scope this
wave (void at the key, one `--full` away). Reference: the trusted ROHF
(−759.73738486, rungs/rohf.chk). The pinned row and its free twin
freeze the reference's OWN core; their residual vs the anchor
(+5.1 mEh) is the projected-density bookkeeping form. Sector
M_S = ½. Sibling: `elements/Sc/` (contracted).

<!-- methods:begin -->
| method | route | window | CAS | core | dets | E, Eh | ΔE, mEh | E, eV | ΔE, eV | conv | ⟨L²⟩ | ptp (Eh) | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | ‖[H,Lz]‖ | ‖[H,L²]‖ | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | — | — | — | — | — | −759.7374 | 0 | −20673.5179 | +0.0000 | — | — | — | — | — | — | — | forced filling |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d / D | (3,5) | Ar | 50 | −759.3360 | +401.3 | −20662.5968 | +10.9210 | — | 11.9994 | 5e-02 | 3e-02 | 4e-02 | 1.2e-02 | 7.2e-02 |  |
| CASCI | free | 3d@r1 / D | (3,5) | Ar·raw | 50 | −759.3375 | +399.9 | −20662.6368 | +10.8810 | — | 11.9994 | 5e-02 | 3e-02 | 4e-02 | 3.6e-02 | 2.1e-01 |  |
| CASCI | free | 3d@r1 / D | (3,5) | Ar | 50 | −759.3375 | +399.9 | −20662.6368 | +10.8810 | — | 11.9994 | 5e-02 | 3e-02 | 4e-02 | 1.2e-02 | 7.2e-02 |  |
| CASCI | free | 3d@r1 / D | (3,5) | Ar·pur | 50 | −759.3374 | +400.0 | −20662.6332 | +10.8847 | — | 12.0000 | 5e-02 | 2e-14 | 2e-14 | 6.6e-15 | 6.5e-14 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (3,6) | Ar·raw | 90 | −759.7342 | +3.2 | −20673.4319 | +0.0860 | — | 5.9990 | 9e-02 | 3e-02 | 4e-02 | 4.6e-02 | 2.5e-01 |  |
| CASCI | free | 3d4s@r1 / D | (3,6) | Ar·raw | 90 | −759.7342 | +3.1 | −20673.4324 | +0.0855 | — | 5.9990 | 9e-02 | 3e-02 | 4e-02 | 4.6e-02 | 2.5e-01 |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (3,6) | Ar | 90 | −759.7342 | +3.2 | −20673.4319 | +0.0860 | — | 5.9989 | 9e-02 | 3e-02 | 4e-02 | 1.4e-02 | 8.4e-02 |  |
| CASCI | free | 3d4s@r1 / D | (3,6) | Ar | 90 | −759.7342 | +3.1 | −20673.4324 | +0.0855 | — | 5.9989 | 9e-02 | 3e-02 | 4e-02 | 1.4e-02 | 8.4e-02 |  |
| CASCI | free | 3d4s / MO | (3,6) | Ar | 90 | −759.7470 | −9.6 | −20673.7794 | −0.2616 | — | 0.0060 | 1e-01 | 2e+00 | 3e+00 | 1.0e-01 | 2.6e-01 | m-incomplete |
| CASCI | free | 3d4s / MOC | (3,6) | Ar | 90 | −759.7395 | −2.1 | −20673.5753 | −0.0574 | — | 1.4464 | 1e-01 | 2e+00 | 4e+00 | 1.5e+00 | 6.5e+00 | m-incomplete |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | free | 3d4s / D | (3,6) | Ar·pur | 90 | −759.7327 | +4.7 | −20673.3903 | +0.1276 | — | 6.0000 | 9e-02 | 2e-14 | 2e-14 | 6.3e-15 | 9.3e-14 |  |
| CASCI | free | 3d4s@r1 / D | (3,6) | Ar·pur | 90 | −759.7327 | +4.7 | −20673.3909 | +0.1270 | — | 6.0000 | 9e-02 | 2e-14 | 2e-14 | 1.0e-14 | 9.2e-14 |  |
| CASCI | free | 3d4s / MO | (3,6) | Ar·pur | 90 | −759.7455 | −8.1 | −20673.7395 | −0.2216 | — | 0.0060 | 1e-01 | 2e+00 | 3e+00 | 1.0e-01 | 2.6e-01 | m-incomplete |
| CASCI | free | 3d4s / MOC | (3,6) | Ar·pur | 90 | −759.7380 | −0.6 | −20673.5355 | −0.0176 | — | 1.4467 | 1e-01 | 2e+00 | 4e+00 | 1.5e+00 | 6.5e+00 | m-incomplete |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASCI | pinned | 3d4s / D | (3,6) | Ar·lane | 90 | −759.7322 | +5.1 | −20673.3781 | +0.1398 | — | 6.0000 | 9e-16 | 1e-14 | 2e-14 | 5.6e-15 | 4.3e-14 | bookkeeping form |
| CASCI | free | 3d4s / D | (3,6) | Ar·lane | 90 | −759.7323 | +5.0 | −20673.3805 | +0.1374 | — | 6.0000 | 9e-02 | 4e-11 | 5e-11 | 6.7e-15 | 4.4e-14 |  |
<!-- methods:end -->
