"""The u14e V&V driver: runs the family legs in order, aggregates
the verdicts, and leaves the rendered table behind.  The derivations
live in the legs -- each one standalone and readable top-to-bottom
(PySCF-direct; the store-free form of the wave's V&V instrument).

This file began as the single-run orchestrator (the production-seam
wiring; its green record: logs/u14e_vv.log at 6e2ab39..30808bd,
2026-11-20) and is now the driver for the per-CAS family:
    u14e_vv_4_12.py   the CAS(4,12) leg, five routes
    u14e_vv_4_15.py   the CAS(4,15) leg (the pinned route measured:
                      the block grounds are not the 5L; the pinned
                      SA prices that ground ensemble, disclosed;
                      src unchanged)
    u14e_vv_6_13.py   the CAS(6,13) disclosure leg, five routes +
                      the PySCF spin_square cross-check
    u14e_vv_uhf.py    the reference comparison + the two table
                      reference rows

Run from the repo root:
    OMP_NUM_THREADS=16 python scripts/u14e/u14e_vv_orchestration.py \
        [--ground-only]  |  tee logs/u14e_vv_family.log
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import u14e_vv_4_12   # noqa: E402
import u14e_vv_4_15   # noqa: E402
import u14e_vv_6_13   # noqa: E402
import u14e_vv_uhf    # noqa: E402


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    bad = [a for a in argv if a not in ("--ground-only",)]
    if bad:
        sys.exit(f"unknown flags: {bad}")
    extra = ["--ground-only"] if "--ground-only" in argv else []
    legs = [("4_12", u14e_vv_4_12), ("4_15", u14e_vv_4_15),
            ("6_13", u14e_vv_6_13), ("uhf", u14e_vv_uhf)]
    bad_rc = []
    for name, leg in legs:
        print(f"\n===== leg {name} =====", flush=True)
        t0 = time.time()
        rc = leg.main(extra)
        print(f"===== leg {name}: {'GREEN' if rc == 0 else 'RED'} "
              f"({time.time() - t0:.0f} s) =====", flush=True)
        if rc != 0:
            bad_rc.append(name)
    if bad_rc:
        print(f"\nRED legs: {bad_rc}", flush=True)
        return 1
    print("\nall legs green", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
