# SA-deviation decomposition: Sc

input-window baselines (one-particle closure defects):
  shell-start      |[P,Lz]| =  5.66e-03   |[P,Lp]| =  2.07e-02
  radial-start     |[P,Lz]| =  5.66e-03   |[P,Lp]| =  2.07e-02
  canonical-start  |[P,Lz]| =  2.80e+00   |[P,Lp]| =  3.73e+00
  purified core substrate (ticket 13; the fact-17 fix):
  pur shell-start  |[P,Lz]| =  6.64e-15   |[P,Lp]| =  8.68e-15
  pur radial-start |[P,Lz]| =  6.64e-15   |[P,Lp]| =  8.68e-15
  core-span closure: |[P,Lz]| 4.5e-14  |[P,Lp]| 7.1e-14  (max foreign-l contamination 5.4e-05)

rung                  |[H,Lz]|  |[H,L2]|  |[P,Lz]|  |[P,Lp]|  sv_min refold_lz  attribution
casci_avg             3.44e-14  7.39e-13  6.64e-15  8.67e-15  1.0000   5.6e-15  shell span (closed)
casscf_relax          2.84e-03  4.04e-02  5.66e-03  2.07e-02  0.9932   6.6e-14  at input baseline
casscf_radial         1.58e+00  3.98e+00  2.70e+00  1.68e+00  0.7849   3.4e-14  identity change; member-tilted
sa                    8.27e-01  3.81e+00  1.39e+00  1.45e+00  0.0000   5.9e-15  identity change; member-tilted
sa_radial             3.28e-03  5.43e-02  6.01e-03  2.18e-02  0.9860   2.9e-14  at input baseline
sa_pyscf              3.41e-01  1.20e+00  2.83e+00  3.99e+00  0.0000   9.2e-14  identity change; member-tilted
casci_avg_pur         6.65e-14  2.67e-13  6.64e-15  8.68e-15  0.9933   3.4e-14  shell span (closed)
casscf_relax_pur      6.65e-14  2.40e-13  6.64e-15  8.68e-15  0.9933   6.6e-14  shell span (closed)
casscf_radial_pur     2.44e-05  4.01e+00  2.14e-05  3.16e+00  0.8143   3.4e-14  identity change; member-tilted
sa_relax_pur          1.59e-13  6.33e-13  6.64e-15  8.68e-15  0.9933   1.6e-13  shell span (closed)
sa_pur                9.20e-02  3.81e+00  1.55e-01  2.00e+00  0.0000   9.6e-14  identity change; member-tilted
sa_radial_pur         2.26e-07  1.61e-05  3.15e-07  3.09e-06  0.9860   3.4e-14  healed below input baseline

AVAS window (Sc 3d, Sc 4s):
  span overlap 0.9855  |[P,Lz]| 4.9e-01  |[P,Lp]| 4.0e-01
  CASCI E -759.73677921 (+80.8 mEh vs casci_avg)  |[H,Lz]| 5.9e-01  |[H,L2]| 2.6e+00  <L^2> 5.9037 (2D)

# Sc tilt-follow curves (frozen, deterministic)

delta  subset          dE(uEh)   |[H,Lz]|  |[H,L2]|  <L^2>  L
 0.05  (0, 1, 2)         -298.3  2.44e-01 1.16e+00   5.989  2.00
 0.10  (0, 1, 2)         -662.3  4.27e-01 2.03e+00   5.953  1.99
 0.15  (0, 1, 2)        -1078.6  5.46e-01 2.59e+00   5.893  1.98
 0.20  (0, 1, 2)        -1521.4  6.05e-01 2.81e+00   5.805  1.96
 0.25  (0, 1, 2)        -1958.3  6.23e-01 2.73e+00   5.682  1.94
 0.30  (0, 1, 2)        -2380.2  6.43e-01 2.48e+00   5.534  1.91
 0.35  (0, 1, 2)        -2745.3  7.29e-01 2.35e+00   5.358  1.87
 0.05  (0, 1, 3, 4)      -409.6  8.17e-15 8.91e-01   5.986  2.00
 0.10  (0, 1, 3, 4)      -910.0  9.60e-15 1.51e+00   5.942  1.99
 0.15  (0, 1, 3, 4)     -1479.6  1.54e-14 1.84e+00   5.867  1.97
 0.20  (0, 1, 3, 4)     -2079.3  2.32e-14 1.86e+00   5.759  1.95
 0.25  (0, 1, 3, 4)     -2660.8  3.49e-14 1.63e+00   5.616  1.92
 0.30  (0, 1, 3, 4)     -3179.5  4.60e-14 1.36e+00   5.437  1.88
 0.35  (0, 1, 3, 4)     -3602.1  4.82e-14 1.58e+00   5.224  1.84

m-symmetric reference (all members together):
  delta 0.05: dE    +2.83 mEh  <L^2> 6.000  (L 2.0)
  delta 0.15: dE   +15.02 mEh  <L^2> 6.000  (L 2.0)
  delta 0.30: dE   +51.55 mEh  <L^2> 6.000  (L 2.0)
