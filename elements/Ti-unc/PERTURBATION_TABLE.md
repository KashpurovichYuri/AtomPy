# Ti-unc -- PT2

<!-- pt2 begin -->
PT2 analyzation: Ti-unc -- nevpt2-sc (pyscf.mrpt.NEVPT, strongly contracted), pyscf 2.13.1; payload collected 2026-11-20T07:38:16.

Reading key: README · Analysis (the reference-class axis, the gauge-spread disclosure, the basis-thinness caveat). E_ref within its own partition (kernel totals only); E_PT2 = E_ref + E(2).  The row reads the GROUND root; the block spread runs over the roots degenerate with it (the gauge block) -- broken spans store excited roots beside it (different states, not gauges; full per-root records live in the payload).


**the span ladder (CASCI: the unoptimized windows)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `ti-unc_3d4s_3-1_casci_ao_-_pur_avg` | 9 | 12.000000 | -848.40183803 | -710.4805 | -19.3332 | -849.11231850 | -23105.5233 | 1.8e-01 |
| `ti-unc_3d4s_3-1_casci_mo_-_pur_avg` | 9 | 0.004805 | -848.41541461 | -698.1801 | -18.9984 | -849.11359470 | -23105.5580 | 0.0e+00 |
| `ti-unc_3d4s_3-1_casci_moc_-_pur_avg` | 9 | 4.904147 | -848.40853049 | -698.4165 | -19.0049 | -849.10694700 | -23105.3771 | 0.0e+00 |
| `ti-unc_3d4s@r1_3-1_casci_ao_-_pur_avg` | 9 | 12.000000 | -848.40183803 | -710.6202 | -19.3370 | -849.11245822 | -23105.5271 | 2.3e-01 |
| `ti-unc_3d4s@r1_3-1_casci_mo_-_pur_avg` | 9 | 0.004805 | -848.41541461 | -698.1809 | -18.9985 | -849.11359551 | -23105.5580 | 0.0e+00 |
| `ti-unc_3d4s@r1_3-1_casci_moc_-_pur_avg` | 9 | 4.904147 | -848.40853049 | -698.4283 | -19.0052 | -849.10695878 | -23105.3774 | 0.0e+00 |
| `ti-unc_3d4s@s1_3-1_casci_ao_-_pur_avg` | 9 | 12.000000 | -848.37657158 | -718.0650 | -19.5395 | -849.09463662 | -23105.0421 | 2.3e-01 |
| `ti-unc_3d4s@s1_3-1_casci_mo_-_pur_avg` | 9 | 0.004805 | -848.41541461 | -698.1819 | -18.9985 | -849.11359647 | -23105.5580 | 0.0e+00 |
| `ti-unc_3d4s@s1_3-1_casci_moc_-_pur_avg` | 9 | 4.904147 | -848.40853049 | -698.4280 | -19.0052 | -849.10695845 | -23105.3774 | 0.0e+00 |
| `ti-unc_3d4s@s2_3-1_casci_ao_-_pur_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 1.827e-08) and CI vectors both fail the stored solution | — | — | — | — |
| `ti-unc_3d4s@s2_3-1_casci_mo_-_pur_avg` | — | — | — | REFUSED: stored-roots refusal: the readout energies (max diff 3.178e-08) and CI vectors both fail the stored solution | — | — | — | — |
| `ti-unc_3d4s@s2_3-1_casci_moc_-_pur_avg` | 9 | 11.149619 | -848.41174989 | -696.2112 | -18.9449 | -849.10796106 | -23105.4047 | 0.0e+00 |
| `ti-unc_3d4s_3-1_sa-casscf_ao_-_pur_avg` | 7 | 12.004661 | -848.41702462 | -704.2259 | -19.1630 | -849.12125048 | -23105.7663 | 2.2e-03 |
| `ti-unc_3d4s_3-1_sa-casscf_mo_-_pur_avg` | 7 | 0.025306 | -848.43303155 | -683.5772 | -18.6011 | -849.11660871 | -23105.6400 | 0.0e+00 |
| `ti-unc_3d4s@r1_3-1_sa-casscf_ao_-_pur_avg` | 7 | 12.004523 | -848.41673278 | -705.6192 | -19.2009 | -849.12235200 | -23105.7963 | 9.1e-04 |
| `ti-unc_3d4s@r1_3-1_sa-casscf_mo_-_pur_avg` | 7 | 0.025306 | -848.43303155 | -683.5753 | -18.6010 | -849.11660684 | -23105.6400 | 0.0e+00 |

**the NO twins (CASCI)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `ti-unc_3d4sno2_3-1_casci_no_-_pur_avg` | 9 | 12.000000 | -848.40183476 | -710.4803 | -19.3332 | -849.11231503 | -23105.5232 | 2.0e-01 |
| `ti-unc_3d4s@r1no2_3-1_casci_no_-_pur_avg` | 9 | 12.000000 | -848.40183803 | -710.5897 | -19.3361 | -849.11242773 | -23105.5262 | 1.8e-01 |
| `ti-unc_3d4s@s1no2_3-1_casci_no_-_pur_avg` | 9 | 12.000000 | -848.37657158 | -717.7928 | -19.5321 | -849.09436440 | -23105.0347 | 2.4e-01 |
| `ti-unc_3d4s@s2no2_3-1_casci_no_-_pur_avg` | 9 | 12.000000 | -848.40194038 | -708.9264 | -19.2909 | -849.11086680 | -23105.4838 | 2.1e-01 |

**the reference axis (single-root CASSCF vs SA)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `ti-unc_3d4s_3-1_casscf_ao_-_pur_avg` | 1 | 3.188134 | -848.42782864 | -669.2511 | -18.2112 | -849.09707973 | -23105.1086 | 0.0e+00 |
| `ti-unc_3d4s_3-1_casscf_mo_-_pur_avg` | 1 | 0.050328 | -848.43881645 | -679.1696 | -18.4811 | -849.11798601 | -23105.6775 | 0.0e+00 |
| `ti-unc_3d4s@r1_3-1_casscf_ao_-_pur_avg` | 1 | 4.167398 | -848.42752793 | -669.2658 | -18.2116 | -849.09679371 | -23105.1008 | 0.0e+00 |
| `ti-unc_3d4s@r1_3-1_casscf_mo_-_pur_avg` | 1 | 0.050328 | -848.43881645 | -679.1703 | -18.4812 | -849.11798680 | -23105.6775 | 0.0e+00 |

**the frozen pair (the core-relaxation disclosure)**

| cell | nroots | ⟨L²⟩ | E_ref (Eh) | E(2) (mEh) | E(2) (eV) | E_PT2 (Eh) | E_PT2 (eV) | block spread (mEh) |
|---|---|---|---|---|---|---|---|---|
| `ti-unc_3d4s_3-1_sa-casscf_ao_frozen_pur_avg` | 7 | 12.000000 | -848.40187119 | -710.2452 | -19.3268 | -849.11211636 | -23105.5178 | 2.3e-01 |
| `ti-unc_3d4s@r1_3-1_sa-casscf_ao_frozen_pur_avg` | 7 | 12.000000 | -848.40187441 | -710.2666 | -19.3273 | -849.11214097 | -23105.5184 | 2.2e-01 |

<!-- payload sha256 fe71d9f66c392b1e... bands: energies 1e-11 absolute; CI amplitudes 1e-8 per-root; the degenerate block's joint projector (p_tol 1e-4, the solver-tail scale); the re-solve angle (davidson noise) 1e-3; the same-state energy noise 1e-6 -->
<!-- pt2 end -->





