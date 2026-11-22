# U14e -- the V&V family

| method | route | CAS | &lt;L²&gt; | &lt;S²&gt; | E_method (Eh) | dE_method (mEh/eV) | E(2) (mEh/eV) | E_tot (Eh) | dE_tot (mEh/eV) |
|---|---|---|---|---|---|---|---|---|---|
| rohf | - | - | - | - | -50.3819 | +0.0000 / +0.0000 | — | -50.3819 | +0.0000 / +0.0000 |
| uhf | - | - | - | 6.2289 (det) | -50.3888 | -6.9611 / -0.1894 | — | -50.3888 | -6.9611 / -0.1894 |
| casci | free | (4,12) | 42.0000 | 6.0000 | -49.4587 | +923.1246 / +25.1195 | -1053.8257 / -28.6761 | -50.5126 | -130.7012 / -3.5566 |
| casscf | free | (4,12) | 71.8167 | 6.0000 | -50.4297 | -47.8788 / -1.3028 | -315.6922 / -8.5904 | -50.7454 | -363.5710 / -9.8933 |
| sa-casscf | free | (4,12) | 72.0000 | 6.0000 | -50.4238 | -41.9595 / -1.1418 | -341.6680 / -9.2973 | -50.7655 | -383.6276 / -10.4390 |
| ci | pinned | (4,12) | 72.0000 | 6.0000 | -48.9351 | +1446.7552 / +39.3682 | -1780.0141 / -48.4367 | -50.7151 | -333.2588 / -9.0684 |
| sa-casscf | pinned | (4,12) | 72.0000 | 6.0000 | -50.4238 | -41.9577 / -1.1417 | -341.7330 / -9.2990 | -50.7655 | -383.6907 / -10.4408 |
| casci | free | (4,15) | 42.0000 | 6.0000 | -49.4587 | +923.1116 / +25.1191 | -1056.2479 / -28.7420 | -50.5150 | -133.1363 / -3.6228 |
| casscf | free | (4,15) | 71.5981 | 6.0000 | -50.4327 | -50.8041 / -1.3824 | -309.8471 / -8.4314 | -50.7425 | -360.6511 / -9.8138 |
| sa-casscf | free | (4,15) | 71.7849..71.9654 | 6.0000 | -50.4264 | -44.5272 / -1.2116 | -336.0093 / -9.1433 | -50.7624 | -380.5365 / -10.3549 |
| ci | pinned | (4,15) | 0.0000..42.0000 | 6.0000 | -48.1294 | +2252.4212 / +61.2915 | -2504.5576 / -68.1525 | -50.6340 | -252.1364 / -6.8610 |
| casci | free | (6,13) | 18.0085 | 12.0000 | -49.3168 | +1065.0798 / +28.9823 | -1081.0321 / -29.4164 | -50.3978 | -15.9523 / -0.4341 |
| casscf | free | (6,13) | 22.2994 | 6.0000 | -50.4494 | -67.5574 / -1.8383 | -307.3006 / -8.3621 | -50.7567 | -374.8580 / -10.2004 |
| sa-casscf | free | (6,13) | 28.3109..44.1466 | 12.0000 | -50.4528 | -70.9662 / -1.9311 | -282.0857 / -7.6759 | -50.7349 | -353.0519 / -9.6070 |
| ci | pinned | (6,13) | 72.0000..90.0000 | 12.0000 | -48.7864 | +1595.4834 / +43.4153 | -1823.7784 / -49.6275 | -50.6101 | -228.2950 / -6.2122 |
| sa-casscf | pinned | (6,13) | 71.9985..89.9988 | 12.0000 | -50.4515 | -69.6548 / -1.8954 | -285.4712 / -7.7681 | -50.7370 | -355.1260 / -9.6635 |

E_method -- the method's own energy (the CASSCF/SA total; the pinned CI value is the TERM energy, the components' mean plus the core bookkeeping constant -- its dE is that constant's projected-density form, not comparable to the free rows).  E_tot = E_method + E(2) (SC-NEVPT2, pyscf.mrpt, strongly contracted); dE against the ROHF row.

&lt;L²&gt;/&lt;S²&gt; -- the ground root for the single-root methods; min..max over the ensemble's 17 roots (free SA) or over the term's pinned components (pinned rows).  The spread is the disclosure (draw-dependent partitions, closure costs, block-ground classes); the per-root record lives in `logs/u14e_vv_*.log`, the committed campaign payload in `elements/U14e/perturbation.json` -- two independent sources, compared here.

The pinned (4,15) row is the RE-ATTEMPT under the corrected symmetrizer, and it MEASURED the refusal: the 14p column's survival norm against the frozen core span is 0.28 (92% core-lean) and the pinned-CI block grounds come out NOT the 5L (L classes 0..42, the 5S/5I classes -- at (4,12) the same route reads 17 components at 72.000000).  The pinned (4,15) SA is measured out (it would optimize that wrong ensemble).  The free rows carry the (4,15) experiment; the project-then-validate promotion into lane_core is declined, src unchanged.
