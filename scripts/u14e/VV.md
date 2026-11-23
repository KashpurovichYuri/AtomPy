# U14e -- the V&V family

| method | route | CAS | &lt;L²&gt; | &lt;S²&gt; | conv | E_method (Eh) | dE_method (mEh/eV) | E(2) (mEh/eV) | E_tot (Eh) | dE_tot (mEh/eV) | term |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rohf | - | - | - | - | — | -50.3819 | +0.0000 / +0.0000 | — | -50.3819 | +0.0000 / +0.0000 | — |
| uhf | - | - | - | 6.2289 (det) | — | -50.3888 | -6.9611 / -0.1894 | — | -50.3888 | -6.9611 / -0.1894 | — |
| casci | free | (4,12) | 42.0000 | 6.0000 | — | -49.4587 | +923.1246 / +25.1195 | -1053.8257 / -28.6761 | -50.5126 | -130.7012 / -3.5566 | 5I input ground |
| casci | free | (4,15) | 42.0000 | 6.0000 | — | -49.4587 | +923.1116 / +25.1191 | -1056.2457 / -28.7419 | -50.5150 | -133.1341 / -3.6228 | 5I input ground |
| casci | free | (6,13) | 12.0000 | 12.0000 | 1.0e-08/— | -49.3168 | +1065.0798 / +28.9823 | -1081.0323 / -29.4164 | -50.3978 | -15.9525 / -0.4341 | definite-L: the 7F, S=3 (the d5 high-spin ground; the 7s column unoccupied) |
| casscf | free | (4,12) | 71.8167 | 6.0000 | — | -50.4297 | -47.8788 / -1.3028 | -315.6922 / -8.5904 | -50.7454 | -363.5710 / -9.8933 | escape draw |
| casscf | free | (4,15) | 71.5977 | 6.0000 | — | -50.4327 | -50.7983 / -1.3823 | -309.8437 / -8.4313 | -50.7425 | -360.6419 / -9.8136 | escape draw |
| casscf | free | (6,13) | 71.3194 | 6.0000 | 1.3e-04/51* | -50.4496 | -67.7098 / -1.8425 | -307.1157 / -8.3570 | -50.7567 | -374.8255 / -10.1995 | escape draw |
| sa-casscf | free | (4,12) | 72.0000 | 6.0000 | — | -50.4238 | -41.9595 / -1.1418 | -341.6680 / -9.2973 | -50.7655 | -383.6276 / -10.4390 | 5L |
| sa-casscf | free | (4,15) | 71.7846..71.9654 | 6.0000 | — | -50.4264 | -44.5275 / -1.2117 | -336.0104 / -9.1433 | -50.7624 | -380.5379 / -10.3550 | 5L, closure-costed |
| sa-casscf | free | (6,13) | 89.9987..89.9999 | 12.0000 | 3.7e-08/13 | -50.4529 | -71.0038 / -1.9321 | -282.0713 / -7.6756 | -50.7349 | -353.0752 / -9.6077 | the L=9 class (90), S=3 |
| casci | pinned | (4,12) | 72.0000 | 6.0000 | — | -48.9351 | +1446.7552 / +39.3682 | -1780.0141 / -48.4367 | -50.7151 | -333.2588 / -9.0684 | 5L (term energy) |
| casci | pinned | (4,15) | 0.0000..42.0000 | 6.0000 | — | -48.1294 | +2252.4212 / +61.2915 | -2504.5576 / -68.1525 | -50.6340 | -252.1364 / -6.8610 | grounds not the 5L (5S/5I) |
| casci | pinned | (6,13) | 72.0000..90.0000 | 12.0000 | 1.0e-12/— | -48.7864 | +1595.4834 / +43.4153 | -1823.7784 / -49.6275 | -50.6101 | -228.2950 / -6.2122 | L=9 ground + the 5L present |
| sa-casscf | pinned | (4,12) | 72.0000 | 6.0000 | — | -50.4238 | -41.9577 / -1.1417 | -341.7330 / -9.2990 | -50.7655 | -383.6907 / -10.4408 | 5L |
| sa-casscf | pinned | (6,13) | 71.9987..89.9989 | 12.0000 | 1.6e-08/94 | -50.4515 | -69.6548 / -1.8954 | -285.4393 / -7.7672 | -50.7369 | -355.0941 / -9.6626 | the 72-component ok |

E_method -- the method's own energy (the CASSCF/SA total; the pinned casci value is the TERM energy, the components' mean plus the core bookkeeping constant -- its dE is that constant's projected-density form, not comparable to the free rows).  E_tot = E_method + E(2) (SC-NEVPT2, pyscf.mrpt, strongly contracted); dE against the ROHF row.

E(2) comparability -- within a row, E_tot = E_method + E(2); across rows, compare E_tot at matched partitions.  The frozen-window rows' large E(2) is the UNOPTIMIZED gauge: the optimizer converts PT2-priced (dynamic) correlation into variationally-absorbed (static) correlation -- read that trade down the E_method/E(2) columns (at (4,12): casci E(2) -1054 mEh -&gt; SA -342 mEh while E_var falls +923 -&gt; -42 mEh; the net is the E_tot column).

&lt;L²&gt;/&lt;S²&gt; -- the ground root for the single-root methods; min..max over the ensemble's 17 roots (free SA) or over the term's pinned components (pinned rows).  The spread is the disclosure (draw-dependent partitions, closure costs, block-ground classes); the per-root record lives in `logs/u14e_vv_*.log`, the committed campaign payload in `elements/U14e/perturbation.json` -- two independent sources, compared here.

conv -- "&lt;precision&gt;/&lt;n_iters&gt;": the optimizer's achieved |dE| (the callback's final macro-step delta; the pinned SA's conv_dE) over its macro-cycle count; the CI solves carry the requested CI tolerance with no iteration structure (—).  ``*`` marks a cycle-cap draw (not converged; the run's log discloses).  Cells read — where a run's record predates the column (the convergence flag lives in that run's log).

The pinned (4,15) row is the RE-ATTEMPT under the corrected symmetrizer, and it MEASURED the refusal: the declared core span of @6s6p7s realizes 6p as THREE bands (13p/14p/15p), and the pinned route's frozen core (the C1 reference's m-mixed occupied 6p projected onto that span) spans the 14p slots ENTIRELY -- the promoted columns are core-span content by construction, the pinned SA has no window at this partition (the typed VanishedColumnError, re-measured every run), and the pinned casci's block grounds come out NOT the 5L (L classes 0..42; at (4,12) the same route reads 17 components at 72.000000).  The pinned (4,15) is the pinned (4,12) here.  The FREE route's core (the pur substrate) does not contain the slots (survival 0.28) and hosts the promotion -- its measured cost +3.09 mEh (fresh, twice) is the experiment's answer; the lane_core promotion is declined, src unchanged.

The (6,13) rows carry the CORRECTED meters (the 2026-11-22 slice incident: the leg's meter slices ran one column off the key's core-pair count; the kernels were never affected).  The frozen ground is the 7F -- L = 3, S = 3 definite (the d5 high-spin class + one f; the 7s column UNOCCUPIED at 3e-03): the as-built window's variational ground abandons the intended 5f3 6d1 7s2 configuration, driven by the weak distillation captures (10f 0.82, 13d 0.92, 7s 0.27-0.71).  The window-rework candidate is the s distillation's k (2 or 3 of the 3 s bands) -- a declaration change, the user's call.  The window itself is healthy: assembled span = projected span exactly, closure meters ~1e-13 (the campaign cell agreed all along).
