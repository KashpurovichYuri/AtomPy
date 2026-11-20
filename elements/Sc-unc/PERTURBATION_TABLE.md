# Sc-unc -- PT2

<!-- pt2 begin -->
PT2 analyzation: Sc-unc -- nevpt2-sc (pyscf.mrpt.NEVPT, strongly contracted), pyscf 2.13.1; payload collected 2026-11-20T07:36:11.

Reading key: README · Analysis (the reference-class axis, the gauge-spread disclosure, the basis-thinness caveat). E_ref within its own partition (kernel totals only); E_PT2 = E_ref + E(2).  The row reads the GROUND root; the block spread runs over the roots degenerate with it (the gauge block) -- broken spans store excited roots beside it (different states, not gauges; full per-root records live in the payload).


**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc-unc_3d4s_2-1_casci_ao_-_pur_avg` | 9 | 6.000000 | -759.73269656 | -659.3525 | -17.9419 | -760.39204909 | -20691.3217 | 5.4e-04 |
| `sc-unc_3d4s_2-1_casci_mo_-_pur_avg` | 9 | 0.006034 | -759.74552859 | -646.4618 | -17.5911 | -760.39199039 | -20691.3201 | 0.0e+00 |
| `sc-unc_3d4s_2-1_casci_moc_-_pur_avg` | 9 | 1.446716 | -759.73803250 | -648.2375 | -17.6394 | -760.38626996 | -20691.1645 | 0.0e+00 |
| `sc-unc_3d4s@r1_2-1_casci_ao_-_pur_avg` | 9 | 6.000000 | -759.73271880 | -659.3749 | -17.9425 | -760.39209367 | -20691.3230 | 5.1e-04 |
| `sc-unc_3d4s@r1_2-1_casci_mo_-_pur_avg` | 9 | 0.006034 | -759.74552859 | -646.4618 | -17.5911 | -760.39199039 | -20691.3201 | 0.0e+00 |
| `sc-unc_3d4s@r1_2-1_casci_moc_-_pur_avg` | 9 | 1.446716 | -759.73803250 | -648.2373 | -17.6394 | -760.38626984 | -20691.1645 | 0.0e+00 |
| `sc-unc_3d4s@s1_2-1_casci_ao_-_pur_avg` | 9 | 6.000000 | -759.71816993 | -669.0429 | -18.2056 | -760.38721281 | -20691.1901 | 1.3e-04 |
| `sc-unc_3d4s@s1_2-1_casci_mo_-_pur_avg` | 9 | 0.006034 | -759.74552859 | -646.4618 | -17.5911 | -760.39199039 | -20691.3201 | 0.0e+00 |
| `sc-unc_3d4s@s1_2-1_casci_moc_-_pur_avg` | 9 | 1.446716 | -759.73803250 | -648.2372 | -17.6394 | -760.38626975 | -20691.1645 | 0.0e+00 |
| `sc-unc_3d4s@s2_2-1_casci_ao_-_pur_avg` | 9 | 6.000000 | -759.74429945 | -644.5652 | -17.5395 | -760.38886469 | -20691.2351 | 1.1e-02 |
| `sc-unc_3d4s@s2_2-1_casci_mo_-_pur_avg` | 9 | 1.528684 | -759.75273130 | -640.7237 | -17.4350 | -760.39345495 | -20691.3600 | 0.0e+00 |
| `sc-unc_3d4s@s2_2-1_casci_moc_-_pur_avg` | 9 | 5.354448 | -759.74293922 | -644.6574 | -17.5420 | -760.38759666 | -20691.2006 | 0.0e+00 |

**the NO twins (CASCI)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc-unc_3d4sno2_2-1_casci_no_-_pur_avg` | 9 | 6.000000 | -759.73269340 | -659.3603 | -17.9421 | -760.39205374 | -20691.3219 | 6.7e-04 |
| `sc-unc_3d4s@r1no2_2-1_casci_no_-_pur_avg` | 9 | 6.000000 | -759.73271880 | -659.3753 | -17.9425 | -760.39209411 | -20691.3230 | 6.0e-04 |
| `sc-unc_3d4s@s1no2_2-1_casci_no_-_pur_avg` | 9 | 6.000000 | -759.71816993 | -669.0429 | -18.2056 | -760.38721282 | -20691.1901 | 2.2e-04 |
| `sc-unc_3d4s@s2no2_2-1_casci_no_-_pur_avg` | 9 | 6.000000 | -759.73297818 | -658.1338 | -17.9087 | -760.39111194 | -20691.2962 | 3.5e-04 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc-unc_3d4s_2-1_casscf_ao_-_pur_avg` | 1 | 1.928348 | -759.76060852 | -634.1635 | -17.2565 | -760.39477197 | -20691.3958 | 0.0e+00 |
| `sc-unc_3d4s_2-1_casscf_mo_-_pur_avg` | 1 | 0.053488 | -759.76869522 | -629.6934 | -17.1348 | -760.39838860 | -20691.4943 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc-unc_3d4s_2-1_sa-casscf_ao_-_pur_avg` | 5 | 6.004746 | -759.74971301 | -652.3723 | -17.7520 | -760.40208528 | -20691.5948 | 0.0e+00 |
| `sc-unc_3d4s_2-1_sa-casscf_mo_-_pur_avg` | 5 | 0.030234 | -759.76364065 | -632.5197 | -17.2117 | -760.39616032 | -20691.4336 | 0.0e+00 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc-unc_3d4s@r1_2-1_casscf_ao_-_pur_avg` | 1 | 1.618035 | -759.75995792 | -635.4685 | -17.2920 | -760.39542645 | -20691.4136 | 0.0e+00 |
| `sc-unc_3d4s@r1_2-1_casscf_mo_-_pur_avg` | 1 | 0.053488 | -759.76869522 | -629.6934 | -17.1348 | -760.39838860 | -20691.4943 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc-unc_3d4s@r1_2-1_sa-casscf_ao_-_pur_avg` | 5 | 6.003676 | -759.74904516 | -654.2138 | -17.8021 | -760.40325896 | -20691.6268 | 0.0e+00 |
| `sc-unc_3d4s@r1_2-1_sa-casscf_mo_-_pur_avg` | 5 | 0.030234 | -759.76364065 | -632.5197 | -17.2117 | -760.39616032 | -20691.4336 | 0.0e+00 |

**the frozen pair (the core-relaxation disclosure)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `sc-unc_3d4s_2-1_sa-casscf_ao_frozen_pur_avg` | 5 | 6.000000 | -759.73274220 | -659.1692 | -17.9369 | -760.39191145 | -20691.3180 | 1.1e-03 |
| `sc-unc_3d4s@r1_2-1_sa-casscf_ao_frozen_pur_avg` | 5 | 6.000000 | -759.73276435 | -659.1963 | -17.9376 | -760.39196066 | -20691.3193 | 9.4e-04 |

<!-- payload sha256 d3349ce54639c055... bands: energies 1e-11 absolute; CI amplitudes 1e-8 per-root; the degenerate block's joint projector (p_tol 1e-4, the solver-tail scale); the re-solve angle (davidson noise) 1e-3; the same-state energy noise 1e-6 -->
<!-- pt2 end -->

