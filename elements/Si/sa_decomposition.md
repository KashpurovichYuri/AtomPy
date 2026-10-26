# SA-deviation decomposition: Si

input-window baselines (one-particle closure defects):
  shell-start      |[P,Lz]| =  3.39e-03   |[P,Lp]| =  4.45e-03
  radial-start     |[P,Lz]| =  3.39e-03   |[P,Lp]| =  4.45e-03
  canonical-start  |[P,Lz]| =  2.43e-01   |[P,Lp]| =  3.19e-01
  purified core substrate (ticket 13; the fact-17 fix):
  pur shell-start  |[P,Lz]| =  1.09e-15   |[P,Lp]| =  1.77e-15
  pur radial-start |[P,Lz]| =  1.09e-15   |[P,Lp]| =  1.77e-15
  core-span closure: |[P,Lz]| 6.5e-16  |[P,Lp]| 1.1e-15  (max foreign-l contamination 9.8e-06)

rung                  |[H,Lz]|  |[H,L2]|  |[P,Lz]|  |[P,Lp]|  sv_min refold_lz  attribution
casci_avg             1.26e-15  0.00e+00  1.15e-15  1.91e-15  1.0000   0.0e+00  shell span (closed)
casscf_relax          6.06e-05  4.06e-21  3.39e-03  4.45e-03  0.9943   0.0e+00  at input baseline
casscf_radial         2.35e-03  2.86e-08  1.18e-01  1.18e-01  0.8833   0.0e+00  identity change; member-tilted
sa                    7.38e-10  1.65e-22  3.98e-06  4.04e-06  0.8833   2.5e-15  healed below input baseline
sa_radial             1.40e-05  2.92e-13  2.33e-03  3.04e-03  0.8832   0.0e+00  at input baseline
sa_pyscf              1.33e-10  2.50e-26  6.94e-08  9.11e-08  0.8833   0.0e+00  healed below input baseline
casci_avg_pur         6.28e-16  0.00e+00  1.09e-15  1.77e-15  0.9943   0.0e+00  shell span (closed)
casscf_relax_pur      0.00e+00  0.00e+00  1.09e-15  1.77e-15  0.9943   0.0e+00  shell span (closed)
casscf_radial_pur     7.85e-16  0.00e+00  1.48e-14  1.30e+00  0.8833   2.5e-15  identity change; member-tilted
sa_relax_pur          0.00e+00  0.00e+00  1.09e-15  1.77e-15  0.9943   0.0e+00  shell span (closed)
sa_pur                3.24e-16  2.96e-32  3.67e-15  1.65e-14  0.8833   0.0e+00  shell span (closed)
sa_radial_pur         1.57e-16  0.00e+00  1.25e-15  1.21e-14  0.8833   1.3e-15  shell span (closed)

AVAS window (Si 3p):
  span overlap 0.8828  |[P,Lz]| 3.2e-03  |[P,Lp]| 4.2e-03
  CASCI E -288.84636276 (-126.7 mEh vs casci_avg)  |[H,Lz]| 1.3e-05  |[H,L2]| 1.7e-21  <L^2> 2.0000 (3P)
