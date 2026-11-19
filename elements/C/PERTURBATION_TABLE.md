# C -- PT2

<!-- pt2 begin -->
PT2 analyzation: C -- nevpt2-sc (pyscf.mrpt.NEVPT, strongly contracted), pyscf 2.13.1; payload collected 2026-11-19T22:52:57.

Reading key: README · Analysis (the reference-class axis, the gauge-spread disclosure, the basis-thinness caveat). E_ref within its own partition (kernel totals only); E_PT2 = E_ref + E(2).  The row reads the GROUND root; the block spread runs over the roots degenerate with it (the gauge block) -- broken spans store excited roots beside it (different states, not gauges; full per-root records live in the payload).


**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2p_2-0_casci_ao_-_pur_avg` | 3 | 2.000000 | -37.43160772 | -265.7585 | -7.2317 | -37.69736624 | -1025.7976 | 1.1e-13 |
| `c_2p_2-0_casci_mo_-_pur_avg` | 3 | 1.917610 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 0.0e+00 |
| `c_2p_2-0_casci_moc_-_pur_avg` | 3 | 0.082390 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 0.0e+00 |
| `c_2p@r1_2-0_casci_ao_-_pur_avg` | 3 | 2.000000 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 6.9e-14 |
| `c_2p@r1_2-0_casci_mo_-_pur_avg` | 3 | 1.917610 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 0.0e+00 |
| `c_2p@r1_2-0_casci_moc_-_pur_avg` | 3 | 0.082390 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 0.0e+00 |
| `c_2p@s1_2-0_casci_ao_-_pur_avg` | 3 | 2.000000 | -37.68219775 | -58.3930 | -1.5890 | -37.74059074 | -1026.9738 | 1.4e-13 |
| `c_2p@s1_2-0_casci_mo_-_pur_avg` | 3 | 1.917610 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 0.0e+00 |
| `c_2p@s1_2-0_casci_moc_-_pur_avg` | 3 | 0.082390 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 0.0e+00 |

**the NO twins (CASCI)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2p@r1no1_2-0_casci_no_-_pur_avg` | 3 | 2.000000 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 9.7e-14 |
| `c_2p@s1no1_2-0_casci_no_-_pur_avg` | 3 | 2.000000 | -37.68219775 | -58.3930 | -1.5890 | -37.74059074 | -1026.9738 | 1.4e-13 |
| `c_2pno1_2-0_casci_no_-_pur_avg` | 3 | 2.000000 | -37.43160772 | -265.7585 | -7.2317 | -37.69736624 | -1025.7976 | 5.0e-13 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2p_2-0_casscf_ao_-_pur_avg` | 1 | 1.683618 | -37.68241788 | -58.0119 | -1.5786 | -37.74042975 | -1026.9694 | 0.0e+00 |
| `c_2p_2-0_casscf_mo_-_pur_avg` | 1 | 1.917610 | -37.68241788 | -58.0119 | -1.5786 | -37.74042974 | -1026.9694 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2p_2-0_sa-casscf_ao_-_pur_avg` | 3 | 2.000000 | -37.68239129 | -58.1340 | -1.5819 | -37.74052528 | -1026.9720 | 2.7e-11 |
| `c_2p_2-0_sa-casscf_mo_-_pur_avg` | 3 | 2.000000 | -37.68239130 | -58.1340 | -1.5819 | -37.74052529 | -1026.9720 | 4.4e-05 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2p@r1_2-0_casscf_ao_-_pur_avg` | 1 | 2.000000 | -37.68241788 | -58.0119 | -1.5786 | -37.74042974 | -1026.9694 | 0.0e+00 |
| `c_2p@r1_2-0_casscf_mo_-_pur_avg` | 1 | 1.917610 | -37.68241788 | -58.0119 | -1.5786 | -37.74042974 | -1026.9694 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2p@r1_2-0_sa-casscf_ao_-_pur_avg` | 3 | 2.000000 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 4.9e-14 |
| `c_2p@r1_2-0_sa-casscf_mo_-_pur_avg` | 3 | 2.000000 | -37.68239130 | -58.1340 | -1.5819 | -37.74052529 | -1026.9720 | 4.4e-05 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2p@s1_2-0_casscf_ao_-_pur_avg` | 1 | 1.998957 | -37.68241788 | -58.0119 | -1.5786 | -37.74042976 | -1026.9694 | 0.0e+00 |
| `c_2p@s1_2-0_casscf_mo_-_pur_avg` | 1 | 1.917610 | -37.68241788 | -58.0119 | -1.5786 | -37.74042974 | -1026.9694 | 0.0e+00 |

**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2p@s1_2-0_sa-casscf_ao_-_pur_avg` | 3 | 2.000000 | -37.68239129 | -58.1339 | -1.5819 | -37.74052517 | -1026.9720 | 4.9e-14 |
| `c_2p@s1_2-0_sa-casscf_mo_-_pur_avg` | 3 | 2.000000 | -37.68239130 | -58.1340 | -1.5819 | -37.74052529 | -1026.9720 | 4.4e-05 |

**the frozen pair (the core-relaxation disclosure)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `c_2p_2-0_sa-casscf_ao_frozen_pur_avg` | 3 | 2.000000 | -37.44038205 | -243.3347 | -6.6215 | -37.68371673 | -1025.4262 | 1.4e-13 |
| `c_2p@r1_2-0_sa-casscf_ao_frozen_pur_avg` | 3 | 2.000000 | -37.68239129 | -58.1341 | -1.5819 | -37.74052539 | -1026.9720 | 1.5e-13 |

<!-- payload sha256 714d5c20e1f82843... bands: energies 1e-11 absolute; CI amplitudes 1e-8 per-root; the degenerate block's joint projector; the re-solve angle (davidson noise) 1e-3 -->
<!-- pt2 end -->

