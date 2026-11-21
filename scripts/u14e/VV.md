# U14e -- the V&V family

One row per method x route x CAS from `scripts/u14e/u14e_vv_*.py` -- standalone, PySCF-direct derivations of the committed campaign answers.  dE is against the ROHF reference (the first row -- the zero of every dE), before the perturbation (dE_var) and after (dE_tot); each in `mEh / eV`. E(2) is the SC-NEVPT2 correction (pyscf.mrpt, strongly contracted).  The fresh-vs-committed assertions and the full per-root record live in `logs/u14e_vv_*.log`; the committed payload is `elements/U14e/perturbation.json`.

| method | route | CAS | &lt;L²&gt; | &lt;S²&gt; | E_var (Eh) | dE_var (mEh/eV) | E(2) (mEh/eV) | E_tot (Eh) | dE_tot (mEh/eV) | date |
|---|---|---|---|---|---|---|---|---|---|---|
| rohf | - | - | - | - | -50.3819 | +0.0000 / +0.0000 | — | -50.3819 | +0.0000 / +0.0000 | 2026-11-21 |
| uhf | - | - | - | 6.2289 (det) | -50.3888 | -6.9611 / -0.1894 | — | -50.3888 | -6.9611 / -0.1894 | 2026-11-21 |
| casscf | free | (4,12) | 71.8171 | 6.0000 | -50.4298 | -47.9941 / -1.3060 | -301.9729 / -8.2171 | -50.4298 | -47.9912 / -1.3059 | 2026-11-21 |
| sa-casscf | free | (4,12) | 72.0000 | 6.0000 | -50.4238 | -41.9595 / -1.1418 | -341.6680 / -9.2973 | -50.7655 | -383.6276 / -10.4390 | 2026-11-21 |
| ci | pinned | (4,12) | 72.0000 | 6.0000 | -48.9351 | +1446.7552 / +39.3682 | -1779.9931 / -48.4361 | -50.7151 | -333.2378 / -9.0679 | 2026-11-21 |
| sa-casscf | pinned | (4,12) | 72.0000 | 6.0000 | -50.4238 | -41.9577 / -1.1417 | -341.7316 / -9.2990 | -50.7655 | -383.6893 / -10.4407 | 2026-11-21 |
| casscf | free | (4,15) | 71.5981 | 6.0000 | -50.4327 | -50.8003 / -1.3823 | -309.8353 / -8.4310 | -50.4327 | -50.8003 / -1.3823 | 2026-11-21 |
| sa-casscf | free | (4,15) | 71.7842..71.9655 | 6.0000 | -50.4264 | -44.5284 / -1.2117 | -336.0203 / -9.1436 | -50.7624 | -380.5487 / -10.3553 | 2026-11-21 |
| casscf | free | (6,13) | 22.3282 | 6.0000 | -50.4494 | -67.5599 / -1.8384 | -307.2710 / -8.3613 | -50.4494 | -67.5182 / -1.8373 | 2026-11-21 |
| sa-casscf | free | (6,13) | 23.3943..44.9792 | 12.0000 | -50.4529 | -71.0036 / -1.9321 | -282.0763 / -7.6757 | -50.7349 | -353.0799 / -9.6078 | 2026-11-21 |
| ci | pinned | (6,13) | 72.0000..90.0000 | 12.0000 | -48.7875 | +1594.3326 / +43.3840 | -1823.7784 / -49.6275 | -50.6113 | -229.4458 / -6.2435 | 2026-11-21 |
| sa-casscf | pinned | (6,13) | 71.9988..89.9986 | 12.0000 | -50.4515 | -69.6548 / -1.8954 | -285.4645 / -7.7679 | -50.7370 | -355.1193 / -9.6633 | 2026-11-21 |

The meters are atomterm's exact angular operators (double-certified); the pinned SA rows carry the run's own ensemble total; the pinned CI rows carry the TERM energy (the components' mean plus the core bookkeeping constant) -- their dE is that constant's projected-density form, not comparable to the free rows.
