# U — uranium (f³d¹ outside the ECP core)

Stuttgart RLC ECP78 + cut-Dyall valence basis; reference **ROHF
−50.381852 Eh** (gated Newton; the reference determinant is a
numerically pure ⁵S° state — NOTES). Sector M_S = 2; pur substrate
mandatory. Provenance: `results/u_14e/u_run*.log` (T3a), `u_sa*.log`,
`u_six13.log` (T3b). NOTES.md: the ticket-23 narrative.

## T3a — Window lattice at matched cost (frozen CASCI; ROHF reference
−50.382 Eh; the first five rows are equal-determinant pairs)

| window | CAS | dets | AO, Eh | ⟨L²⟩ | MO, Eh | ⟨L²⟩ | NO, Eh | ⟨L²⟩ |
|---|---|---|---|---|---|---|---|---|
| 1f+1d ‖ open+8 virt. | (4,12) | 495 | −49.459 | 42.000 | −50.382 | 0.000 | −49.890 | 72.000 |
| 2f+1d ‖ open+15 | (4,19) | 3,876 | −49.751 | 72.000 | −50.383 | 0.132 | — | — |
| 1f+3d ‖ open+18 | (4,22) | 7,315 | −49.733 | 42.000 | −50.384 | 0.192 | — | — |
| 2f+2d ‖ open+20 | (4,24) | 10,626 | −49.894 | 72.000 | −50.384 | 0.254 | −49.894 | 72.000 |
| 3f+1d ‖ open+22 | (4,26) | 14,950 | −49.993 | 72.000 | −50.384 | 0.258 | — | — |
| 3f+2d | (4,31) | 31,465 | −50.081 | 72.000 | — | — | — | — |
| 4f+1d | (4,33) | 40,920 | −50.332 | 72.000 | — | — | — | — |
| 4f+2d | (4,38) | 73,815 | **−50.392** | 72.000 | — | — | — | — |

Descent race: Δf pays f³d¹ −339/−330 mEh but f²d² only −187; Δd pays
f²d² −144 but f³d¹ −69/−60. At 4f+2d the ⁵L° (f³d¹) manifold stands
below the best f²d² root by > 37 mEh — ground term **⁵L° (f³d¹,
⟨L²⟩ = 72)**.

## T3b — Treatment ladder on the NO substrate (frozen-CASCI rows carry
our fold: E_core validated −0.000 mEh against the campaign totals;
optimizer rows are as-run totals; labels from exact eigh — PySCF ci
vectors straggle in the 17-fold ⁵L cluster)

| treatment | window | E, Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | note |
|---|---|---|---|---|---|---|
| CASCI (anchor) | NO d5+f7, 495 dets | −49.889731 | 0 | 72.0000 | ~1e-13 | machine-exact |
| frozen growth | all 4 NO groups, 10,626 | −49.894287 | −4.6 | 72.0000 | ~1e-13 | |
| CASSCF, core-relax | NO-12 | −49.899729 | −10.0 | 72.0000 | 1.3e-13 | machine-exact closure |
| CASSCF, core-relax | NO-17 | −49.900366 | −10.6 | 72.0000 | 1.3e-12 | |
| SA-CASSCF (17×⁵L) | NO-12 | **−50.423810** | −534.1 | **72.0000** | 1.3e-04 | **below the ⁵S° mean field at exact ⁵L — radial relaxation** |
| SA-CASSCF (17×⁵L) | NO-17 | ≈ −50.428 | −537 | 71.88 | 1.0e+01 | degraded: span opens |
| CASSCF, free | any start | ≈ −50.43 | ≈ −540 | broken | 1e+01…1e+02 | never converged; no term |
| CASCI (7s² promoted) | NO-12 + 7s, CAS(6,13) | −49.894866 | −5.1 | 72.0000 | ~1e-13 | 7s² inert (⟨N7s⟩ = 1.97) |

The SA-12 mechanism (measured): the window's l-content is exactly
preserved (l2 = 5.0, l3 = 7.0) while 3.53 columns of weight move to
f/d primitives outside the selected exponent bands — ticket-24's
missing-radial axis recovered by relaxation.

## Open questions

⁵S° mechanism (AOC hypothesis); inner core-pair valence (the −551 mEh
accidental 6p probe, to be done m-complete); f³d¹-vs-f²d² ordering in
one NO window; NIST; 3rd-d point.
