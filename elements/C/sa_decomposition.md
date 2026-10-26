# SA-deviation decomposition: C

input-window baselines (one-particle closure defects):
  shell-start      |[P,Lz]| =  7.35e-16   |[P,Lp]| =  9.86e-16
  radial-start     |[P,Lz]| =  6.56e-16   |[P,Lp]| =  9.28e-16
  canonical-start  |[P,Lz]| =  6.44e-02   |[P,Lp]| =  5.70e-01
  purified core substrate (ticket 13; the fact-17 fix):
  pur shell-start  |[P,Lz]| =  6.55e-16   |[P,Lp]| =  9.26e-16
  pur radial-start |[P,Lz]| =  6.55e-16   |[P,Lp]| =  9.26e-16
  core-span closure: |[P,Lz]| 0.0e+00  |[P,Lp]| 9.0e-17  (max foreign-l contamination 6.3e-06)

rung                  |[H,Lz]|  |[H,L2]|  |[P,Lz]|  |[P,Lp]|  sv_min refold_lz  attribution
casci_avg             1.26e-15  0.00e+00  6.55e-16  9.26e-16  1.0000   6.3e-16  shell span (closed)
casscf_relax          1.26e-15  1.39e-31  7.35e-16  9.86e-16  1.0000   1.3e-15  shell span (closed)
casscf_radial         5.07e-07  4.19e-34  8.25e-06  1.05e+00  0.9177   0.0e+00  member-tilted
sa                    1.00e-10  4.61e-20  7.74e-05  6.85e-05  0.9177   1.3e-15  radial relaxation
sa_radial             8.27e-12  2.63e-43  1.23e-06  9.54e-05  0.9177   6.3e-16  radial relaxation
sa_pyscf              6.81e-13  1.57e-28  6.33e-07  5.60e-06  0.9177   1.3e-15  radial relaxation
casci_avg_pur         0.00e+00  0.00e+00  6.55e-16  9.26e-16  1.0000   0.0e+00  shell span (closed)
casscf_relax_pur      0.00e+00  0.00e+00  6.55e-16  9.26e-16  1.0000   0.0e+00  shell span (closed)
casscf_radial_pur     0.00e+00  0.00e+00  1.31e-15  1.12e+00  0.9177   0.0e+00  member-tilted
sa_relax_pur          0.00e+00  0.00e+00  6.55e-16  9.26e-16  1.0000   0.0e+00  shell span (closed)
sa_pur                6.28e-16  0.00e+00  1.74e-11  2.60e-11  0.9177   0.0e+00  shell span (closed)
sa_radial_pur         0.00e+00  0.00e+00  5.50e-16  1.81e-15  0.9177   0.0e+00  shell span (closed)

AVAS window (C 2p):
  span overlap 0.9173  |[P,Lz]| 2.8e-04  |[P,Lp]| 2.5e-03
  CASCI E -37.68236462 (-250.8 mEh vs casci_avg)  |[H,Lz]| 1.2e-07  |[H,L2]| 3.7e-23  <L^2> 2.0000 (3P)
