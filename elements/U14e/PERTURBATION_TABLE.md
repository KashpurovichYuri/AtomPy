# U14e -- PT2

<!-- pt2 begin -->
PT2 analyzation: U14e -- nevpt2-sc (pyscf.mrpt.NEVPT, strongly contracted), pyscf 2.13.1; payload collected 2026-11-20T08:39:21.

Reading key: README · Analysis (the reference-class axis, the gauge-spread disclosure, the basis-thinness caveat). E_ref within its own partition (kernel totals only); E_PT2 = E_ref + E(2).  The row reads the GROUND root; the block spread runs over the roots degenerate with it (the gauge block) -- broken spans store excited roots beside it (different states, not gauges; full per-root records live in the payload).


**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `u14e_10f13d_4-0_casci_ao_-_pur_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 3.428e-08) and CI vectors both fail the stored solution | — | — | — | — |
| `u14e_10f13d_4-0_casci_mo_-_pur_avg` | 9 | 0.001862 | -50.37811374 | -358.2550 | -9.7486 | -50.73636877 | -1380.6069 | 0.0e+00 |
| `u14e_10f13d_4-0_casci_moc_-_pur_avg` | 9 | 0.016028 | -50.37832810 | -352.0219 | -9.5790 | -50.73034996 | -1380.4431 | 0.0e+00 |
| `u14e_10f13d7s_5-1_casci_ao_-_pur_avg` | 9 | 12.000000 | -49.31639422 | -1082.7940 | -29.4643 | -50.39918824 | -1371.4318 | 1.1e+00 |
| `u14e_10f13d7s_5-1_casci_mo_-_pur_avg` | 9 | 0.034793 | -50.39099589 | -344.7942 | -9.3823 | -50.73579005 | -1380.5912 | 0.0e+00 |
| `u14e_10f13d7s_5-1_casci_moc_-_pur_avg` | 9 | 0.000005 | -50.38024000 | -345.4660 | -9.4006 | -50.72570595 | -1380.3168 | 0.0e+00 |
| `u14e_10f13d14p_4-0_casci_ao_-_pur_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 3.063e-08) and CI vectors both fail the stored solution | — | — | — | — |
| `u14e_10f13d14p_4-0_casci_mo_-_pur_avg` | 9 | 3.466372 | -50.37822402 | -359.0757 | -9.7709 | -50.73729969 | -1380.6323 | 0.0e+00 |
| `u14e_10f13d14p_4-0_casci_moc_-_pur_avg` | 9 | 0.006650 | -50.37808140 | -348.8181 | -9.4918 | -50.72689950 | -1380.3493 | 0.0e+00 |

**the NO twins (CASCI)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `u14e_10f13dno2_4-0_casci_no_-_pur_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 8.703e-08) and CI vectors both fail the stored solution | — | — | — | — |
| `u14e_10f13d14pno3_4-0_casci_no_-_pur_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 9.383e-09) and CI vectors both fail the stored solution | — | — | — | — |

**the distilled prior family (CASCI ao)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `u14e_5f6d@r1_4-0_casci_ao_-_pur_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 7.214e-08) and CI vectors both fail the stored solution | — | — | — | — |
| `u14e_5f6d@s1_4-0_casci_ao_-_pur_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 4.925e-08) and CI vectors both fail the stored solution | — | — | — | — |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `u14e_10f13d_4-0_casscf_ao_-_pur_avg` | 1 | 65.799091 | -50.42946634 | -311.0458 | -8.4640 | -50.74051212 | -1380.7197 | 0.0e+00 |
| `u14e_10f13d_4-0_casscf_mo_-_pur_avg` | 1 | 64.443470 | -50.43064119 | -317.8351 | -8.6487 | -50.74847625 | -1380.9364 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `u14e_10f13d_4-0_sa-casscf_ao_-_pur_avg` | 17 | 71.999999 | -50.42381143 | -341.6651 | -9.2972 | -50.76547655 | -1381.3990 | 8.1e-02 |
| `u14e_10f13d_4-0_sa-casscf_mo_-_pur_avg` | 17 | 71.999995 | -50.42381662 | -341.6782 | -9.2975 | -50.76549478 | -1381.3995 | 2.8e-03 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `u14e_10f13d7s_5-1_casscf_ao_-_pur_avg` | 1 | 87.474451 | -50.45379277 | -278.2798 | -7.5724 | -50.73207253 | -1380.4900 | 0.0e+00 |
| `u14e_10f13d7s_5-1_casscf_mo_-_pur_avg` | 1 | 26.761879 | -50.43979295 | -324.6899 | -8.8353 | -50.76448289 | -1381.3720 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `u14e_10f13d7s_5-1_sa-casscf_ao_-_pur_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 2.325e-05) exceed the same-state noise band (1e-06) | — | — | — | — |
| `u14e_10f13d7s_5-1_sa-casscf_mo_-_pur_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 3.695e-03) exceed the same-state noise band (1e-06) | — | — | — | — |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `u14e_10f13d14p_4-0_casscf_ao_-_pur_avg` | 1 | 64.708978 | -50.43446231 | -309.1237 | -8.4117 | -50.74358598 | -1380.8033 | 0.0e+00 |
| `u14e_10f13d14p_4-0_casscf_mo_-_pur_avg` | 1 | 64.306769 | -50.43445358 | -310.6789 | -8.4540 | -50.74513243 | -1380.8454 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `u14e_10f13d14p_4-0_sa-casscf_ao_-_pur_avg` | 17 | 71.785121 | -50.42637850 | -336.0125 | -9.1434 | -50.76239101 | -1381.3150 | 2.1e-03 |
| `u14e_10f13d14p_4-0_sa-casscf_mo_-_pur_avg` | 17 | 71.782785 | -50.42638234 | -336.0442 | -9.1442 | -50.76242657 | -1381.3160 | 0.0e+00 |

**the frozen pair (the core-relaxation disclosure)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `u14e_10f13d_4-0_sa-casscf_ao_frozen_pur_avg` | 17 | 41.999987 | -49.55618155 | -921.3769 | -25.0719 | -50.47755841 | -1373.5643 | 7.7e-07 |

**the pinned lane (the core-partition axis; the D2h components)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `u14e_5f6d@6s6p7s_4-0_casci_ao_-_pur_avg` | — | — | — | REFUSED: lane re-solve refusal: u14e_5f6d@6s6p7s_4-0_casci_ao_-_pur_avg block Au energies 9.495e-01 off the stored record | — | — | — | — |
| `u14e_5f6d@6s6p7s_4-0_sa-casscf_ao_-_pur_avg` ·Au | 5 | 72.000000 | -49.15486025 | -341.7480 | -9.2994 | -49.49660826 | -1346.8713 | 9.5e-02 |
| `u14e_5f6d@6s6p7s_4-0_sa-casscf_ao_-_pur_avg` ·B1u | 4 | 72.000000 | -49.15486022 | -341.7614 | -9.2998 | -49.49662163 | -1346.8717 | 9.5e-02 |
| `u14e_5f6d@6s6p7s_4-0_sa-casscf_ao_-_pur_avg` ·B2u | 4 | 72.000000 | -49.15486024 | -341.7434 | -9.2993 | -49.49660363 | -1346.8712 | 9.5e-02 |
| `u14e_5f6d@6s6p7s_4-0_sa-casscf_ao_-_pur_avg` ·B3u | 4 | 72.000000 | -49.15486023 | -341.7576 | -9.2997 | -49.49661784 | -1346.8716 | 9.5e-02 |
| `u14e_5f6d@6s6p7s_4-0_casci_ao_-_pol_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 1.938e-02) and CI vectors both fail the stored solution | — | — | — | — |
| `u14e_5f6d7s@6s6p_5-1_casci_ao_-_pur_avg` ·Au | 5 | 72.000000 | -48.76450438 | -1835.7932 | -49.9545 | -50.60029763 | -1376.9042 | 1.5e-01 |
| `u14e_5f6d7s@6s6p_5-1_casci_ao_-_pur_avg` ·B1u | 4 | 72.000000 | -48.76450438 | -1835.7306 | -49.9528 | -50.60023499 | -1376.9025 | 1.5e-01 |
| `u14e_5f6d7s@6s6p_5-1_casci_ao_-_pur_avg` ·B2u | 4 | 72.000000 | -48.76450438 | -1835.7551 | -49.9534 | -50.60025951 | -1376.9032 | 1.5e-01 |
| `u14e_5f6d7s@6s6p_5-1_casci_ao_-_pur_avg` ·B3u | 4 | 72.000000 | -48.76450438 | -1835.8711 | -49.9566 | -50.60037551 | -1376.9064 | 1.5e-01 |
| `u14e_5f6d7s@6s6p_5-1_sa-casscf_ao_-_pur_avg` | — | — | — | REFUSED: lane re-solve refusal: u14e_5f6d7s@6s6p_5-1_sa-casscf_ao_-_pur_avg block B1u energies 8.113e-08 off the stored record | — | — | — | — |
| `u14e_5f6d7s@6s6p_5-1_casci_ao_-_pol_avg` | 9 | 55.999879 | -48.84977469 | -1709.1719 | -46.5089 | -50.55894659 | -1375.7790 | 0.0e+00 |

<!-- payload sha256 7788dd34ad4d72ba... bands: energies 1e-11 absolute; CI amplitudes 1e-8 per-root; the degenerate block's joint projector (p_tol 1e-4, the solver-tail scale); the re-solve angle (davidson noise) 1e-3; the same-state energy noise 1e-6 -->
<!-- pt2 end -->
