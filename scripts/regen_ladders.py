"""Regenerate the campaign ladders under the unified convention
(ticket 25/27): run_ladder for C/Si/Sc/Ti into results/, then print
the T3b cross-check values (max_de_vs_pyscf) for the element tables.

Usage: python scripts/regen_ladders.py
Threads: OPENBLAS_NUM_THREADS=64 OMP_NUM_THREADS=64 PYSCF_NUM_THREADS=64.
"""

import sys
from pathlib import Path

sys.path.insert(0, "src")
from atomterm.system import AtomConfig
from atomterm.ladder import run_ladder

CASES = {
    "C": (AtomConfig(), ["2p"], 2, [1 / 3] * 3, 2, 3),
    "Si": (AtomConfig("Si", "ccpvdz", spin=2), ["3p"], 2, [1 / 3] * 3, 2, 3),
    "Sc": (AtomConfig("Sc", "ccpvdz", spin=1), ["3d", "4s"], 3,
           [1 / 5] * 5, 3, 6),
    "Ti": (AtomConfig("Ti", "ccpvdz", spin=2), ["3d", "4s"], 4,
           [1 / 7] * 7, 4, 6),
}

import sys as _sys
only = _sys.argv[1:] or None
for name, (cfg, patterns, nelec, saw, fo, no) in CASES.items():
    if only and name not in only:
        continue
    print(f"\n===== ladder {name} =====", flush=True)
    out = Path("results") / f"{name}_ccpvdz_s{cfg.spin}_ladder"
    res = run_ladder(cfg, patterns, nelec, workdir=str(out),
                     sa_weights=saw, frac_nelec_open=fo, frac_n_open=no,
                     purified=True)
    rohf_e = None
    for tag, r in res.items():
        ours = r.get("e_ground")
        pyscf = (r.get("energies_pyscf") or [None])[0]
        if pyscf is not None and "rohf" in tag:
            rohf_e = pyscf
    for tag, r in res.items():
        ours = r.get("e_ground")
        pyscf = (r.get("energies_pyscf") or [None])[0]
        if ours is None or pyscf is None:
            continue
        xchk = r.get("max_de_vs_pyscf")
        dE = (pyscf - rohf_e) * 1e3 if rohf_e else float("nan")
        print(f"  {tag:18s} as-run {pyscf:.6f}  dE {dE:+8.1f} mEh"
              f"  xchk {xchk * 1e3:7.1f}" if xchk is not None else
              f"  {tag:18s} as-run {pyscf:.6f}  dE {dE:+8.1f}")

print("\ndone")
