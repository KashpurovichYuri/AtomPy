# C -- PT2

<!-- pt2 begin -->
PT2 analyzation: C -- nevpt2-sc (pyscf.mrpt.NEVPT, strongly contracted), pyscf 2.13.1; payload collected 2026-11-20T05:47:42.

Reading key: README · Analysis (the reference-class axis, the gauge-spread disclosure, the basis-thinness caveat). E_ref within its own partition (kernel totals only); E_PT2 = E_ref + E(2).  The row reads the GROUND root; the block spread runs over the roots degenerate with it (the gauge block) -- broken spans store excited roots beside it (different states, not gauges; full per-root records live in the payload).


**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2p_2-0_casci_ao_-_pur_avg` | 3 | 2.000000 | -37.43160772 | -265.7585 | -7.2317 | -37.69736624 | -1025.7976 | 1.1e-13 |
| `c_2p_2-0_casci_mo_-_pur_avg` | 3 | 1.917610 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 0.0e+00 |
| `c_2p_2-0_casci_moc_-_pur_avg` | 3 | 0.082390 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 0.0e+00 |
| `c_2p@r1_2-0_casci_ao_-_pur_avg` | 3 | 2.000000 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 1.3e-13 |
| `c_2p@r1_2-0_casci_mo_-_pur_avg` | 3 | 1.917610 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 0.0e+00 |
| `c_2p@r1_2-0_casci_moc_-_pur_avg` | 3 | 0.082390 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 0.0e+00 |
| `c_2p@s1_2-0_casci_ao_-_pur_avg` | 3 | 2.000000 | -37.68219775 | -58.3930 | -1.5890 | -37.74059074 | -1026.9738 | 1.2e-13 |
| `c_2p@s1_2-0_casci_mo_-_pur_avg` | 3 | 1.917610 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 0.0e+00 |
| `c_2p@s1_2-0_casci_moc_-_pur_avg` | 3 | 0.082390 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 0.0e+00 |
| `c_2p_2-0_sa-casscf_ao_-_pur_avg` | 3 | 2.000000 | -37.68239129 | -58.1340 | -1.5819 | -37.74052528 | -1026.9720 | 2.7e-11 |
| `c_2p_2-0_sa-casscf_mo_-_pur_avg` | 3 | 2.000000 | -37.68239130 | -58.1340 | -1.5819 | -37.74052529 | -1026.9720 | 4.4e-05 |
| `c_2p@r1_2-0_sa-casscf_ao_-_pur_avg` | 3 | 2.000000 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 4.9e-14 |
| `c_2p@r1_2-0_sa-casscf_mo_-_pur_avg` | 3 | 2.000000 | -37.68239130 | -58.1340 | -1.5819 | -37.74052529 | -1026.9720 | 4.4e-05 |
| `c_2p@s1_2-0_sa-casscf_ao_-_pur_avg` | 3 | 2.000000 | -37.68239129 | -58.1339 | -1.5819 | -37.74052517 | -1026.9720 | 3.5e-14 |
| `c_2p@s1_2-0_sa-casscf_mo_-_pur_avg` | 3 | 2.000000 | -37.68239130 | -58.1340 | -1.5819 | -37.74052529 | -1026.9720 | 4.4e-05 |

**the NO twins (CASCI)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2pno1_2-0_casci_no_-_pur_avg` | 3 | 2.000000 | -37.43160772 | -265.7585 | -7.2317 | -37.69736624 | -1025.7976 | 2.8e-13 |
| `c_2p@r1no1_2-0_casci_no_-_pur_avg` | 3 | 2.000000 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 5.6e-14 |
| `c_2p@s1no1_2-0_casci_no_-_pur_avg` | 3 | 2.000000 | -37.68219775 | -58.3930 | -1.5890 | -37.74059074 | -1026.9738 | 8.3e-14 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2p_2-0_casscf_ao_-_pur_avg` | 1 | 1.683618 | -37.68241788 | -58.0119 | -1.5786 | -37.74042975 | -1026.9694 | 0.0e+00 |
| `c_2p_2-0_casscf_mo_-_pur_avg` | 1 | 1.917610 | -37.68241788 | -58.0119 | -1.5786 | -37.74042974 | -1026.9694 | 0.0e+00 |
| `c_2p@r1_2-0_casscf_ao_-_pur_avg` | 1 | 2.000000 | -37.68241788 | -58.0119 | -1.5786 | -37.74042974 | -1026.9694 | 0.0e+00 |
| `c_2p@r1_2-0_casscf_mo_-_pur_avg` | 1 | 1.917610 | -37.68241788 | -58.0119 | -1.5786 | -37.74042974 | -1026.9694 | 0.0e+00 |
| `c_2p@s1_2-0_casscf_ao_-_pur_avg` | 1 | 1.998957 | -37.68241788 | -58.0119 | -1.5786 | -37.74042976 | -1026.9694 | 0.0e+00 |
| `c_2p@s1_2-0_casscf_mo_-_pur_avg` | 1 | 1.917610 | -37.68241788 | -58.0119 | -1.5786 | -37.74042974 | -1026.9694 | 0.0e+00 |

**the frozen pair (the core-relaxation disclosure)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2p_2-0_sa-casscf_ao_frozen_pur_avg` | 3 | 2.000000 | -37.44038205 | -243.3347 | -6.6215 | -37.68371673 | -1025.4262 | 1.1e-13 |
| `c_2p@r1_2-0_sa-casscf_ao_frozen_pur_avg` | 3 | 2.000000 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 9.7e-14 |

**the pinned lane (the core-partition axis; the D2h components)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2p@1s2s_2-0_casci_ao_-_pur_avg` ·B1g | 1 | 2.000000 | -37.06779261 | -682.3334 | -18.5672 | -37.75012597 | -1027.2333 | 2.6e+02 |
| `c_2p@1s2s_2-0_casci_ao_-_pur_avg` ·B2g | 1 | 2.000000 | -37.06779261 | -421.3158 | -11.4646 | -37.48910838 | -1020.1306 | 2.6e+02 |
| `c_2p@1s2s_2-0_casci_ao_-_pur_avg` ·B3g | 1 | 2.000000 | -37.06779261 | -421.3158 | -11.4646 | -37.48910838 | -1020.1306 | 2.6e+02 |
| `c_2p@1s2s_2-0_sa-casscf_ao_-_pur_avg` ·B1g | 1 | 2.000000 | -37.39724685 | -58.1340 | -1.5819 | -37.45538085 | -1019.2128 | 1.7e+01 |
| `c_2p@1s2s_2-0_sa-casscf_ao_-_pur_avg` ·B2g | 1 | 2.000000 | -37.39724685 | -75.4218 | -2.0523 | -37.47266864 | -1019.6833 | 1.7e+01 |
| `c_2p@1s2s_2-0_sa-casscf_ao_-_pur_avg` ·B3g | 1 | 2.000000 | -37.39724685 | -75.4218 | -2.0523 | -37.47266864 | -1019.6833 | 1.7e+01 |
| `c_2p@1s2s_2-0_casci_ao_-_pol_avg` | 3 | 2.000000 | -37.06779261 | -682.3334 | -18.5672 | -37.75012597 | -1027.2333 | 8.9e-13 |

<!-- payload sha256 937d0edb86e7b280... bands: energies 1e-11 absolute; CI amplitudes 1e-8 per-root; the degenerate block's joint projector (p_tol 1e-6, the solver-tail scale); the re-solve angle (davidson noise) 1e-3 -->
<!-- pt2 end -->





