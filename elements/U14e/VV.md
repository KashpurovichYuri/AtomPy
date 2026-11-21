# U14e -- the V&V family

One row per leg of `scripts/u14e/u14e_vv_*.py` (the standalone, PySCF-direct derivations), beside the committed payload anchor (`perturbation.json`) it asserts against. Re-rendered by every leg after it drops its row; the full record lives in `logs/u14e_vv_*.log`.

| leg | route | CAS | E_ref (Eh) | E(2) (mEh) | total (Eh) | anchor E(2) (mEh) | delta (mEh) | &lt;L²&gt; | &lt;S²&gt; | date |
|---|---|---|---|---|---|---|---|---|---|---|
| ref | uhf | CAS(4,12) | -50.38881343 | +0.0000 | -50.38881343 | — | -6.9611 | 0.032886..0.158732 | <S2> 6.2289 (UHF det) | 2026-11-21 |
| 4_12 | free | CAS(4,12) | -50.42381186 | -341.6680 | -50.76547989 | -341.6651 | -0.0029 | 71.999997..72.000000 | 5.3e-15 dev | 2026-11-21 |
| 4_12 | pinned | CAS(4,12) | -50.42381003 | -341.7481 | -50.76555816 | — | -0.0000 | 72.000000..72.000000 | 72-block exact | 2026-11-21 |
| 4_15 | free | CAS(4,15) | -50.42637972 | -336.0132 | -50.76239288 | -336.0125 | -0.0006 | 71.784763..71.965391 | 5.3e-15 dev | 2026-11-21 |
| 6_13 | free | CAS(6,13) | -50.45324468 | -308.4104 | -50.76165510 | — | +0.0000 | 38.336722..50.373537 | L=9 disclosure | 2026-11-21 |

Atomic units; E(2) = SC-NEVPT2 (pyscf.mrpt, strongly contracted); the meters are atomterm's exact angular operators.
