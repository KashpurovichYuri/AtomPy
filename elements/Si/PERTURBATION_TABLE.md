# Si -- PT2

<!-- pt2 begin -->
PT2 analyzation: Si -- nevpt2-sc (pyscf.mrpt.NEVPT, strongly contracted), pyscf 2.13.1; payload collected 2026-11-20T00:23:53.

Reading key: README · Analysis (the reference-class axis, the gauge-spread disclosure, the basis-thinness caveat). E_ref within its own partition (kernel totals only); E_PT2 = E_ref + E(2).  The row reads the GROUND root; the block spread runs over the roots degenerate with it (the gauge block) -- broken spans store excited roots beside it (different states, not gauges; full per-root records live in the payload).


**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `si_3p_2-0_casci_ao_-_pur_avg` | 3 | 2.000000 | -288.62480814 | -236.0901 | -6.4243 | -288.86089824 | -7860.3055 | 2.5e-06 |
| `si_3p_2-0_casci_mo_-_pur_avg` | 3 | 1.959719 | -288.84639706 | -55.9263 | -1.5218 | -288.90232334 | -7861.4327 | 0.0e+00 |
| `si_3p_2-0_casci_moc_-_pur_avg` | 3 | 0.040281 | -288.84639707 | -55.9263 | -1.5218 | -288.90232334 | -7861.4327 | 0.0e+00 |
| `si_3p@r1_2-0_casci_ao_-_pur_avg` | 3 | 2.000000 | -288.84639706 | -55.9263 | -1.5218 | -288.90232335 | -7861.4327 | 1.5e-08 |
| `si_3p@r1_2-0_casci_mo_-_pur_avg` | 3 | 1.959719 | -288.84639706 | -55.9263 | -1.5218 | -288.90232334 | -7861.4327 | 0.0e+00 |
| `si_3p@r1_2-0_casci_moc_-_pur_avg` | 3 | 0.040281 | -288.84639707 | -55.9263 | -1.5218 | -288.90232334 | -7861.4327 | 0.0e+00 |
| `si_3p@s1_2-0_casci_ao_-_pur_avg` | 3 | 2.000000 | -288.84542446 | -57.1258 | -1.5545 | -288.90255027 | -7861.4389 | 1.4e-06 |
| `si_3p@s1_2-0_casci_mo_-_pur_avg` | 3 | 1.959719 | -288.84639706 | -55.9263 | -1.5218 | -288.90232334 | -7861.4327 | 0.0e+00 |
| `si_3p@s1_2-0_casci_moc_-_pur_avg` | 3 | 0.040281 | -288.84639707 | -55.9263 | -1.5218 | -288.90232334 | -7861.4327 | 0.0e+00 |
| `si_3p@s2_2-0_casci_ao_-_pur_avg` | 9 | 2.000000 | -288.84769214 | -54.5464 | -1.4843 | -288.90223856 | -7861.4304 | 1.1e-12 |
| `si_3p@s2_2-0_casci_mo_-_pur_avg` | 9 | 1.957408 | -288.84769214 | -54.5464 | -1.4843 | -288.90223856 | -7861.4304 | 0.0e+00 |
| `si_3p@s2_2-0_casci_moc_-_pur_avg` | 9 | 1.997293 | -288.84769214 | -54.5464 | -1.4843 | -288.90223856 | -7861.4304 | 0.0e+00 |

**the NO twins (CASCI)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `si_3pno1_2-0_casci_no_-_pur_avg` | 3 | 2.000000 | -288.62480814 | -236.0901 | -6.4243 | -288.86089824 | -7860.3055 | 1.5e-06 |
| `si_3p@r1no1_2-0_casci_no_-_pur_avg` | 3 | 2.000000 | -288.84639706 | -55.9263 | -1.5218 | -288.90232335 | -7861.4327 | 1.8e-07 |
| `si_3p@s1no1_2-0_casci_no_-_pur_avg` | 3 | 2.000000 | -288.84542446 | -57.1258 | -1.5545 | -288.90255027 | -7861.4389 | 9.9e-09 |
| `si_3p@s2no1_2-0_casci_no_-_pur_avg` | 3 | 2.000000 | -288.84639688 | -55.9178 | -1.5216 | -288.90231468 | -7861.4325 | 1.8e-07 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `si_3p_2-0_casscf_ao_-_pur_avg` | 1 | 1.577297 | -288.84643692 | -55.7949 | -1.5183 | -288.90223182 | -7861.4302 | 0.0e+00 |
| `si_3p_2-0_casscf_mo_-_pur_avg` | 1 | 1.961217 | -288.84643692 | -55.7949 | -1.5183 | -288.90223183 | -7861.4302 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `si_3p_2-0_sa-casscf_ao_-_pur_avg` | 3 | 2.000000 | -288.84639779 | -55.9284 | -1.5219 | -288.90232614 | -7861.4328 | 7.0e-07 |
| `si_3p_2-0_sa-casscf_mo_-_pur_avg` | 3 | 2.000000 | -288.84639779 | -55.9284 | -1.5219 | -288.90232618 | -7861.4328 | 2.0e-05 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `si_3p@r1_2-0_casscf_ao_-_pur_avg` | 1 | 1.999982 | -288.84643692 | -55.7949 | -1.5183 | -288.90223183 | -7861.4302 | 0.0e+00 |
| `si_3p@r1_2-0_casscf_mo_-_pur_avg` | 1 | 1.961217 | -288.84643692 | -55.7949 | -1.5183 | -288.90223183 | -7861.4302 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `si_3p@r1_2-0_sa-casscf_ao_-_pur_avg` | 3 | 2.000000 | -288.84639779 | -55.9284 | -1.5219 | -288.90232619 | -7861.4328 | 1.0e-07 |
| `si_3p@r1_2-0_sa-casscf_mo_-_pur_avg` | 3 | 2.000000 | -288.84639779 | -55.9284 | -1.5219 | -288.90232618 | -7861.4328 | 2.0e-05 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `si_3p@s1_2-0_casscf_ao_-_pur_avg` | 1 | 1.998665 | -288.84643692 | -55.7949 | -1.5183 | -288.90223185 | -7861.4302 | 0.0e+00 |
| `si_3p@s1_2-0_casscf_mo_-_pur_avg` | 1 | 1.961217 | -288.84643692 | -55.7949 | -1.5183 | -288.90223183 | -7861.4302 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `si_3p@s1_2-0_sa-casscf_ao_-_pur_avg` | 3 | 2.000000 | -288.84639779 | -55.9284 | -1.5219 | -288.90232617 | -7861.4328 | 2.8e-07 |
| `si_3p@s1_2-0_sa-casscf_mo_-_pur_avg` | 3 | 2.000000 | -288.84639779 | -55.9284 | -1.5219 | -288.90232618 | -7861.4328 | 2.0e-05 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `si_3p@s2_2-0_casscf_ao_-_pur_avg` | 1 | 1.999997 | -288.84773043 | -54.4230 | -1.4809 | -288.90215343 | -7861.4281 | 0.0e+00 |
| `si_3p@s2_2-0_casscf_mo_-_pur_avg` | 1 | 1.958720 | -288.84773043 | -54.4229 | -1.4809 | -288.90215337 | -7861.4281 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `si_3p@s2_2-0_sa-casscf_ao_-_pur_avg` | 3 | 2.000000 | -288.84769322 | -54.5525 | -1.4844 | -288.90224568 | -7861.4306 | 1.0e-12 |
| `si_3p@s2_2-0_sa-casscf_mo_-_pur_avg` | 3 | 1.997292 | -288.84769263 | -54.5514 | -1.4844 | -288.90224403 | -7861.4305 | 0.0e+00 |

**the frozen pair (the core-relaxation disclosure)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `si_3p_2-0_sa-casscf_ao_frozen_pur_avg` | 3 | 2.000000 | -288.62879188 | -226.0793 | -6.1519 | -288.85487120 | -7860.1415 | 4.2e-08 |
| `si_3p@r1_2-0_sa-casscf_ao_frozen_pur_avg` | 3 | 2.000000 | -288.84639707 | -55.9263 | -1.5218 | -288.90232339 | -7861.4327 | 7.5e-07 |

**the pinned lane (the core-partition axis; the D2h components)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `si_3p@1s2s2p3s_2-0_casci_ao_-_pur_avg` ·B1g | 1 | 2.000000 | -288.48505323 | -381.9356 | -10.3930 | -288.86698886 | -7860.4712 | 1.5e+02 |
| `si_3p@1s2s2p3s_2-0_casci_ao_-_pur_avg` ·B2g | 1 | 2.000000 | -288.48505323 | -229.4711 | -6.2442 | -288.71452431 | -7856.3224 | 1.5e+02 |
| `si_3p@1s2s2p3s_2-0_casci_ao_-_pur_avg` ·B3g | 1 | 2.000000 | -288.48505323 | -229.4730 | -6.2443 | -288.71452623 | -7856.3225 | 1.5e+02 |
| `si_3p@1s2s2p3s_2-0_sa-casscf_ao_-_pur_avg` ·B1g | 1 | 2.000000 | -288.73887754 | -55.9284 | -1.5219 | -288.79480592 | -7858.5070 | 5.5e+00 |
| `si_3p@1s2s2p3s_2-0_sa-casscf_ao_-_pur_avg` ·B2g | 1 | 2.000000 | -288.73887754 | -61.4637 | -1.6725 | -288.80034126 | -7858.6576 | 5.5e+00 |
| `si_3p@1s2s2p3s_2-0_sa-casscf_ao_-_pur_avg` ·B3g | 1 | 2.000000 | -288.73887754 | -61.4637 | -1.6725 | -288.80034121 | -7858.6576 | 5.5e+00 |
| `si_3p@1s2s2p3s_2-0_casci_ao_-_pol_avg` | 3 | 2.000000 | -288.48505323 | -381.9356 | -10.3930 | -288.86698887 | -7860.4712 | 4.1e-06 |

<!-- payload sha256 ae40607ff4d28ae1... bands: energies 1e-11 absolute; CI amplitudes 1e-8 per-root; the degenerate block's joint projector (p_tol 1e-6, the solver-tail scale); the re-solve angle (davidson noise) 1e-3 -->
<!-- pt2 end -->


