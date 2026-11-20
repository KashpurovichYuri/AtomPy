# Sc -- PT2

<!-- pt2 begin -->
PT2 analyzation: Sc -- nevpt2-sc (pyscf.mrpt.NEVPT, strongly contracted), pyscf 2.13.1; payload collected 2026-11-20T06:55:58.

Reading key: README · Analysis (the reference-class axis, the gauge-spread disclosure, the basis-thinness caveat). E_ref within its own partition (kernel totals only); E_PT2 = E_ref + E(2).  The row reads the GROUND root; the block spread runs over the roots degenerate with it (the gauge block) -- broken spans store excited roots beside it (different states, not gauges; full per-root records live in the payload).


**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc_3d4s_2-1_casci_ao_-_pur_avg` | 9 | 6.000000 | -759.72889088 | -197.9095 | -5.3854 | -759.92680041 | -20678.6617 | 8.3e-04 |
| `sc_3d4s_2-1_casci_mo_-_pur_avg` | 9 | 0.006211 | -759.74621598 | -181.8475 | -4.9483 | -759.92806344 | -20678.6961 | 0.0e+00 |
| `sc_3d4s_2-1_casci_moc_-_pur_avg` | 9 | 3.353210 | -759.73774256 | -183.0813 | -4.9819 | -759.92082383 | -20678.4991 | 0.0e+00 |
| `sc_3d4s@r1_2-1_casci_ao_-_pur_avg` | 9 | 6.000000 | -759.73584067 | -184.2700 | -5.0142 | -759.92011063 | -20678.4796 | 4.1e-03 |
| `sc_3d4s@r1_2-1_casci_mo_-_pur_avg` | 9 | 0.006211 | -759.74621598 | -181.8475 | -4.9483 | -759.92806344 | -20678.6961 | 0.0e+00 |
| `sc_3d4s@r1_2-1_casci_moc_-_pur_avg` | 9 | 3.353210 | -759.73774256 | -183.0813 | -4.9819 | -759.92082383 | -20678.4991 | 0.0e+00 |
| `sc_3d4s@s1_2-1_casci_ao_-_pur_avg` | 9 | 6.000000 | -759.72342039 | -200.1072 | -5.4452 | -759.92352763 | -20678.5726 | 2.9e-03 |
| `sc_3d4s@s1_2-1_casci_mo_-_pur_avg` | 9 | 0.006211 | -759.74621598 | -181.8475 | -4.9483 | -759.92806344 | -20678.6961 | 0.0e+00 |
| `sc_3d4s@s1_2-1_casci_moc_-_pur_avg` | 9 | 3.353210 | -759.73774256 | -183.0813 | -4.9819 | -759.92082383 | -20678.4991 | 0.0e+00 |
| `sc_3d4s@s2_2-1_casci_ao_-_pur_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 4.674e-09) and CI vectors both fail the stored solution | — | — | — | — |
| `sc_3d4s@s2_2-1_casci_mo_-_pur_avg` | 9 | 1.631128 | -759.74847438 | -180.2857 | -4.9058 | -759.92876012 | -20678.7150 | 0.0e+00 |
| `sc_3d4s@s2_2-1_casci_moc_-_pur_avg` | 9 | 3.976734 | -759.74915898 | -179.5226 | -4.8851 | -759.92868161 | -20678.7129 | 0.0e+00 |

**the NO twins (CASCI)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc_3d4sno2_2-1_casci_no_-_pur_avg` | 9 | 6.000000 | -759.72889088 | -197.9093 | -5.3854 | -759.92680017 | -20678.6617 | 6.7e-04 |
| `sc_3d4s@r1no2_2-1_casci_no_-_pur_avg` | 9 | 6.000000 | -759.73584067 | -184.2699 | -5.0142 | -759.92011061 | -20678.4796 | 4.1e-03 |
| `sc_3d4s@s1no2_2-1_casci_no_-_pur_avg` | 9 | 6.000000 | -759.72342039 | -200.1069 | -5.4452 | -759.92352728 | -20678.5726 | 2.6e-03 |
| `sc_3d4s@s2no2_2-1_casci_no_-_pur_avg` | 9 | 6.000000 | -759.73577496 | -187.3751 | -5.0987 | -759.92315008 | -20678.5624 | 5.9e-03 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc_3d4s_2-1_casscf_ao_-_pur_avg` | 1 | 0.052212 | -759.76776901 | -165.3314 | -4.4989 | -759.93310040 | -20678.8331 | 0.0e+00 |
| `sc_3d4s_2-1_casscf_mo_-_pur_avg` | 1 | 0.052779 | -759.76751435 | -165.5189 | -4.5040 | -759.93303329 | -20678.8313 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc_3d4s_2-1_sa-casscf_ao_-_pur_avg` | 5 | 6.000000 | -759.73973940 | -178.8731 | -4.8674 | -759.91861249 | -20678.4389 | 1.9e-01 |
| `sc_3d4s_2-1_sa-casscf_mo_-_pur_avg` | 5 | 0.030352 | -759.76257017 | -169.0268 | -4.5995 | -759.93159693 | -20678.7922 | 0.0e+00 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc_3d4s@r1_2-1_casscf_ao_-_pur_avg` | 1 | 0.052214 | -759.76776890 | -165.3308 | -4.4989 | -759.93309966 | -20678.8331 | 0.0e+00 |
| `sc_3d4s@r1_2-1_casscf_mo_-_pur_avg` | 1 | 0.052779 | -759.76751435 | -165.5189 | -4.5040 | -759.93303329 | -20678.8313 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc_3d4s@r1_2-1_sa-casscf_ao_-_pur_avg` | 5 | 6.000000 | -759.73973934 | -178.7477 | -4.8640 | -759.91848701 | -20678.4355 | 2.2e-01 |
| `sc_3d4s@r1_2-1_sa-casscf_mo_-_pur_avg` | 5 | 0.030352 | -759.76257017 | -169.0268 | -4.5995 | -759.93159693 | -20678.7922 | 0.0e+00 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc_3d4s@s1_2-1_casscf_ao_-_pur_avg` | 1 | 0.052791 | -759.76751436 | -165.5181 | -4.5040 | -759.93303245 | -20678.8313 | 0.0e+00 |
| `sc_3d4s@s1_2-1_casscf_mo_-_pur_avg` | 1 | 0.052779 | -759.76751435 | -165.5189 | -4.5040 | -759.93303329 | -20678.8313 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc_3d4s@s1_2-1_sa-casscf_ao_-_pur_avg` | 5 | 6.000000 | -759.73973929 | -178.7235 | -4.8633 | -759.91846277 | -20678.4348 | 1.3e-01 |
| `sc_3d4s@s1_2-1_sa-casscf_mo_-_pur_avg` | 5 | 0.030352 | -759.76257017 | -169.0268 | -4.5995 | -759.93159693 | -20678.7922 | 0.0e+00 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc_3d4s@s2_2-1_casscf_ao_-_pur_avg` | 1 | 2.663913 | -759.77400184 | -160.0508 | -4.3552 | -759.93405260 | -20678.8590 | 0.0e+00 |
| `sc_3d4s@s2_2-1_casscf_mo_-_pur_avg` | 1 | 1.723734 | -759.77384809 | -160.2229 | -4.3599 | -759.93407102 | -20678.8595 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc_3d4s@s2_2-1_sa-casscf_ao_-_pur_avg` | 5 | 6.000000 | -759.74818864 | -167.3962 | -4.5551 | -759.91558488 | -20678.3565 | 1.7e-03 |
| `sc_3d4s@s2_2-1_sa-casscf_mo_-_pur_avg` | 5 | 4.013531 | -759.76852485 | -164.5925 | -4.4788 | -759.93311734 | -20678.8336 | 0.0e+00 |

**the frozen pair (the core-relaxation disclosure)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc_3d4s_2-1_sa-casscf_ao_frozen_pur_avg` | 5 | 6.000000 | -759.72897770 | -197.6128 | -5.3773 | -759.92659054 | -20678.6560 | 6.8e-04 |
| `sc_3d4s@r1_2-1_sa-casscf_ao_frozen_pur_avg` | 5 | 6.000000 | -759.73584134 | -184.2576 | -5.0139 | -759.92009890 | -20678.4793 | 4.1e-03 |

**the pinned lane (the core-partition axis; the D2h components)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc_3d4s@ar_2-1_casci_ao_-_pur_avg` ·Ag | 2 | 6.000000 | -759.72743563 | -200.2098 | -5.4480 | -759.92764541 | -20678.6847 | 1.9e+01 |
| `sc_3d4s@ar_2-1_casci_ao_-_pur_avg` ·B1g | 1 | 6.000000 | -759.72743563 | -182.1603 | -4.9568 | -759.90959589 | -20678.1935 | 1.9e+01 |
| `sc_3d4s@ar_2-1_casci_ao_-_pur_avg` ·B2g | 1 | 6.000000 | -759.72743563 | -200.6938 | -5.4612 | -759.92812942 | -20678.6978 | 1.9e+01 |
| `sc_3d4s@ar_2-1_casci_ao_-_pur_avg` ·B3g | 1 | 6.000000 | -759.72743563 | -200.5651 | -5.4577 | -759.92800076 | -20678.6943 | 1.9e+01 |
| `sc_3d4s@ar_2-1_sa-casscf_ao_-_pur_avg` ·Ag | 2 | 6.000000 | -764.54287351 | -178.8304 | -4.8662 | -764.72170389 | -20809.1377 | 5.1e+01 |
| `sc_3d4s@ar_2-1_sa-casscf_ao_-_pur_avg` ·B1g | 1 | 6.000000 | -764.54287206 | -229.4582 | -6.2439 | -764.77233024 | -20810.5153 | 5.1e+01 |
| `sc_3d4s@ar_2-1_sa-casscf_ao_-_pur_avg` ·B2g | 1 | 6.000000 | -764.54287220 | -179.3670 | -4.8808 | -764.72223920 | -20809.1522 | 5.1e+01 |
| `sc_3d4s@ar_2-1_sa-casscf_ao_-_pur_avg` ·B3g | 1 | 6.000000 | -764.54287198 | -179.5957 | -4.8870 | -764.72246768 | -20809.1584 | 5.1e+01 |
| `sc_3d4s@ar_2-1_casci_ao_-_pol_avg` | 9 | 6.000000 | -759.72743563 | -200.2103 | -5.4480 | -759.92764595 | -20678.6847 | 1.3e-03 |
| `sc_3d@ar_2-1_casci_ao_-_pur_avg` ·Ag | 2 | 12.000000 | -759.44716252 | -342.0353 | -9.3073 | -759.78919786 | -20674.9173 | 8.1e+01 |
| `sc_3d@ar_2-1_casci_ao_-_pur_avg` ·B1g | 1 | 12.000000 | -759.44716252 | -356.8592 | -9.7106 | -759.80402167 | -20675.3207 | 8.1e+01 |
| `sc_3d@ar_2-1_casci_ao_-_pur_avg` ·B2g | 1 | 12.000000 | -759.44716252 | -345.0769 | -9.3900 | -759.79223942 | -20675.0001 | 8.1e+01 |
| `sc_3d@ar_2-1_casci_ao_-_pur_avg` ·B3g | 1 | 12.000000 | -759.44716252 | -346.6376 | -9.4325 | -759.79380010 | -20675.0426 | 8.1e+01 |
| `sc_3d@ar_2-1_sa-casscf_ao_-_pur_avg` ·Ag | 2 | 5.981897 | -759.73374444 | -189.9948 | -5.1700 | -759.92373926 | -20678.5784 | 2.0e+01 |
| `sc_3d@ar_2-1_sa-casscf_ao_-_pur_avg` ·B1g | 1 | 5.981897 | -759.73522160 | -178.8383 | -4.8664 | -759.91405992 | -20678.3150 | 2.0e+01 |
| `sc_3d@ar_2-1_sa-casscf_ao_-_pur_avg` ·B2g | 1 | 5.981897 | -759.72400658 | -185.8290 | -5.0567 | -759.90983559 | -20678.2000 | 2.0e+01 |
| `sc_3d@ar_2-1_sa-casscf_ao_-_pur_avg` ·B3g | 1 | 5.981897 | -759.73522261 | -187.8679 | -5.1121 | -759.92309049 | -20678.5607 | 2.0e+01 |
| `sc_3d@ar_2-1_casci_ao_-_pol_avg` | 9 | 12.000000 | -759.44716252 | -347.4108 | -9.4535 | -759.79457332 | -20675.0636 | 2.1e+00 |
| `sc_3d3p@semi_4-3_casci_ao_-_pur_avg` | — | — | — | REFUSED: Lowdin on a near-linearly-dependent block: min overlap eigenvalue 1.04e-17 <= tol 1e-10; pad/drop the needle directions or raise tol | — | — | — | — |
| `sc_3d3p@semi_4-3_sa-casscf_ao_-_pur_avg` | — | — | — | REFUSED: Lowdin on a near-linearly-dependent block: min overlap eigenvalue 1.04e-17 <= tol 1e-10; pad/drop the needle directions or raise tol | — | — | — | — |
| `sc_3d3p@semi_4-3_casci_ao_-_pol_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 2.999e-09) and CI vectors both fail the stored solution | — | — | — | — |
| `sc_3d3p4d@semi_4-3_casci_ao_-_pur_avg` | — | — | — | REFUSED: Lowdin on a near-linearly-dependent block: min overlap eigenvalue 1.04e-17 <= tol 1e-10; pad/drop the needle directions or raise tol | — | — | — | — |
| `sc_3d3p4d@semi_4-3_sa-casscf_ao_-_pur_avg` | — | — | — | REFUSED: Lowdin on a near-linearly-dependent block: min overlap eigenvalue 1.04e-17 <= tol 1e-10; pad/drop the needle directions or raise tol | — | — | — | — |
| `sc_3d3p4d@semi_4-3_casci_ao_-_pol_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 6.457e-11) and CI vectors both fail the stored solution | — | — | — | — |
| `sc_3d3p4s@1s2s2p3s_5-4_casci_ao_-_pur_avg` ·Ag | 2 | 6.000000 | -759.78987426 | -132.5423 | -3.6067 | -759.92241658 | -20678.5424 | 6.0e+00 |
| `sc_3d3p4s@1s2s2p3s_5-4_casci_ao_-_pur_avg` ·B1g | 1 | 6.000000 | -759.78987426 | -135.3522 | -3.6831 | -759.92522644 | -20678.6189 | 6.0e+00 |
| `sc_3d3p4s@1s2s2p3s_5-4_casci_ao_-_pur_avg` ·B2g | 1 | 6.000000 | -759.78987426 | -138.5178 | -3.7693 | -759.92839208 | -20678.7050 | 6.0e+00 |
| `sc_3d3p4s@1s2s2p3s_5-4_casci_ao_-_pur_avg` ·B3g | 1 | 6.000000 | -759.78987426 | -137.4777 | -3.7410 | -759.92735200 | -20678.6767 | 6.0e+00 |
| `sc_3d3p4s@1s2s2p3s_5-4_sa-casscf_ao_-_pur_avg` ·Ag | 2 | 6.000000 | -764.84913763 | -93.9772 | -2.5572 | -764.94311480 | -20815.1626 | 3.9e+01 |
| `sc_3d3p4s@1s2s2p3s_5-4_sa-casscf_ao_-_pur_avg` ·B1g | 1 | 6.000000 | -764.84913449 | -130.3343 | -3.5466 | -764.97946874 | -20816.1518 | 3.9e+01 |
| `sc_3d3p4s@1s2s2p3s_5-4_sa-casscf_ao_-_pur_avg` ·B2g | 1 | 6.000000 | -764.84913456 | -132.6923 | -3.6107 | -764.98182690 | -20816.2160 | 3.9e+01 |
| `sc_3d3p4s@1s2s2p3s_5-4_sa-casscf_ao_-_pur_avg` ·B3g | 1 | 6.000000 | -764.84913426 | -132.1850 | -3.5969 | -764.98131923 | -20816.2021 | 3.9e+01 |
| `sc_3d3p4s@1s2s2p3s_5-4_casci_ao_-_pol_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 2.381e-09) and CI vectors both fail the stored solution | — | — | — | — |
| `sc_3s3p3d@1s2s2p4s_5-4_casci_ao_-_pur_avg` ·Ag | 2 | 6.000000 | -759.78998491 | -138.5834 | -3.7710 | -759.92856830 | -20678.7098 | 3.9e+00 |
| `sc_3s3p3d@1s2s2p4s_5-4_casci_ao_-_pur_avg` ·B1g | 1 | 6.000000 | -759.78998491 | -134.7027 | -3.6654 | -759.92468759 | -20678.6042 | 3.9e+00 |
| `sc_3s3p3d@1s2s2p4s_5-4_casci_ao_-_pur_avg` ·B2g | 1 | 6.000000 | -759.78998491 | -135.3222 | -3.6823 | -759.92530714 | -20678.6211 | 3.9e+00 |
| `sc_3s3p3d@1s2s2p4s_5-4_casci_ao_-_pur_avg` ·B3g | 1 | 6.000000 | -759.78998491 | -135.3022 | -3.6818 | -759.92528710 | -20678.6205 | 3.9e+00 |
| `sc_3s3p3d@1s2s2p4s_5-4_sa-casscf_ao_-_pur_avg` ·Ag | 2 | 6.000000 | -760.10795498 | -94.0887 | -2.5603 | -760.20204363 | -20686.1514 | 4.0e+00 |
| `sc_3s3p3d@1s2s2p4s_5-4_sa-casscf_ao_-_pur_avg` ·B1g | 1 | 6.000000 | -760.10795396 | -96.1397 | -2.6161 | -760.20409365 | -20686.2072 | 4.0e+00 |
| `sc_3s3p3d@1s2s2p4s_5-4_sa-casscf_ao_-_pur_avg` ·B2g | 1 | 6.000000 | -760.10795395 | -95.3687 | -2.5951 | -760.20332265 | -20686.1862 | 4.0e+00 |
| `sc_3s3p3d@1s2s2p4s_5-4_sa-casscf_ao_-_pur_avg` ·B3g | 1 | 6.000000 | -760.10795427 | -98.0531 | -2.6682 | -760.20600735 | -20686.2593 | 4.0e+00 |
| `sc_3s3p3d@1s2s2p4s_5-4_casci_ao_-_pol_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 3.838e-09) and CI vectors both fail the stored solution | — | — | — | — |

<!-- payload sha256 1c8ef1cd68123f56... bands: energies 1e-11 absolute; CI amplitudes 1e-8 per-root; the degenerate block's joint projector (p_tol 1e-6, the solver-tail scale); the re-solve angle (davidson noise) 1e-3 -->
<!-- pt2 end -->

