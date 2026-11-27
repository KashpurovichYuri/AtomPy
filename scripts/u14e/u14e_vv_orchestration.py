"""The u14e V&V driver: the family's legs in order, each standalone
and readable top-to-bottom, the verdicts aggregated.

    u14e_vv_ladder.py   THE WRITER of the table's method rows: every
                        case of the CAS ladder through the free trio
                        and the pinned SA
    u14e_vv_uhf.py      the two reference rows + the UHF sibling's
                        spin honesty (notes only)

The library the legs read is ``scripts/u14e/vv``; the table is
``scripts/u14e/VV.md``.  This driver forwards its flags to the
ladder leg.

Run from the repo root:
    OMP_NUM_THREADS=64 python scripts/u14e/u14e_vv_orchestration.py \\
        [--cases TAG ...] [--routes ...] [--ground-only] [--base TAG] \\
        |  tee logs/u14e_vv_family.log
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import u14e_vv_ladder  # noqa: E402
import u14e_vv_uhf  # noqa: E402

LEGS = (("ladder", u14e_vv_ladder), ("uhf", u14e_vv_uhf))


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    red = []
    for name, leg in LEGS:
        print(f"\n===== leg {name} =====", flush=True)
        t0 = time.time()
        # only the ladder leg takes flags
        rc = leg.main(argv if name == "ladder" else [])
        print(f"===== leg {name}: {'GREEN' if rc == 0 else 'RED'} "
              f"({time.time() - t0:.0f} s) =====", flush=True)
        if rc != 0:
            red.append(name)
    if red:
        print(f"\nRED legs: {red}", flush=True)
        return 1
    print("\nall legs green", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
