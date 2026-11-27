"""The HF-virtual channel screen -- which of the reference's virtual
shells carries correlation?  Printed only (no store writes).

For each declared channel (``vv.channels.CHANNELS``): its radial set
beside the frozen base window, then the frozen CASCI (the free route,
no optimizer) with NEVPT2 on the ground root.  Read per channel
against the base built the same way without channels:

  the extended window's orthonormality and rotation closure; the
  frozen ground's term (<L²>, <S²> -- the base's);
  dE_CAS   the frozen variational energy's change;
  dE(2)    the perturbative correction's change;
  dE_tot   their sum.

A channel that converts PT2 into the CAS reads dE_CAS < 0 with dE(2)
> 0 and dE_tot ~ 0; a channel where PT2 was deficient reads
dE_tot < 0.  The existing (4,15)@r1 and (4,16)@r1 rows of the table
are the controls for the band-named 14p and 7s channels.

Run from the repo root:
    OMP_NUM_THREADS=64 python scripts/u14e/u14e_vv_screen.py \\
        [--base TAG] [--channels 7p 8s ...] [--spectrum] \\
        | tee logs/u14e_vv_screen.log
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import numpy as np  # noqa: E402

from vv import cases as ladder  # noqa: E402
from vv.channels import CHANNELS, channel_block, virtual_spectrum  # noqa: E402
from vv.reference import Substrate  # noqa: E402
from vv.routes import free_casci  # noqa: E402
from vv.windows import window_with_channels  # noqa: E402


def _parse(argv) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="the HF-virtual channel screen on a frozen base "
                    "window (frozen CASCI + NEVPT2 per channel)")
    p.add_argument("--base", default="(4,12)@r1", metavar="TAG",
                   help="the base case (default: (4,12)@r1)")
    p.add_argument("--channels", nargs="+", default=[],
                   metavar="SHELL", help="the channels to screen "
                   f"(default: {list(CHANNELS)})")
    p.add_argument("--spectrum", action="store_true",
                   help="print the reference's virtual spectrum and stop")
    return p.parse_args(argv)


def print_spectrum(sub) -> None:
    """The reference's virtuals by angular character (the inventory's
    measured basis)."""
    print(f"{'idx':>4} {'dE(Eh)':>8}  l(weight)        <L2>     <r2>")
    for r in virtual_spectrum(sub, below=0.9):
        print(f"{r['index']:4d} {r['de']:8.4f}  {r['l']} {r['weight']:.3f}      "
              f"{r['l2']:8.3f} {r['r2']:8.1f}")
    print("declared channels: " + ", ".join(
        f"{c.shell}{list(c.columns)}" for c in CHANNELS.values()))


def frozen(sub, space) -> dict:
    """The frozen reading of one window: the ground root priced."""
    t0 = time.time()
    priced = free_casci(sub, space).price().states[0]
    return {"ecas": priced.state.e_ref, "e2": priced.e2,
            "l2": priced.l2, "s2": priced.s2, "seconds": time.time() - t0}


def main(argv=None) -> int:
    args = _parse(sys.argv[1:] if argv is None else argv)
    sub = Substrate.load()
    if args.spectrum:
        print_spectrum(sub)
        return 0
    names = args.channels or list(CHANNELS)
    unknown = [n for n in names if n not in CHANNELS]
    if unknown:
        sys.exit(f"unknown channel(s) {unknown}; declared: {list(CHANNELS)}")
    base = ladder.CASES[args.base].window(sub)
    b = frozen(sub, base.space)
    print(f"\n===== the base {args.base}: {base.space.ncas} active columns; "
          f"frozen E_CAS {b['ecas']:.8f}, E(2) {b['e2'] * 1e3:+.4f} mEh, "
          f"<L2> {b['l2']:.4f}, <S2> {b['s2']:.4f} ({b['seconds']:.0f} s) =====",
          flush=True)
    rows = []
    for name in names:
        ch = CHANNELS[name]
        print(f"\n===== channel {name} =====", flush=True)
        try:
            block, lines = channel_block(sub, ch)
            built = window_with_channels(sub, base, [(block, lines)])
        except ValueError as e:
            print(f"  [refused] {e}", flush=True)
            rows.append((name, None))
            continue
        space = built.space
        for line in lines:
            print(f"  {line}", flush=True)
        ortho = space.certify(sub.mf)
        cl = space.closure(sub.mol)
        print(f"  window: {space.ncas} active columns; max|C^T S C - I| "
              f"{ortho:.1e}; ||[P,Lz]|| {cl['lz']:.2e}, ||[P,Lperp]|| "
              f"{cl['lperp']:.2e}", flush=True)
        r = frozen(sub, space)
        d = {"ecas": r["ecas"] - b["ecas"], "e2": r["e2"] - b["e2"]}
        d["tot"] = d["ecas"] + d["e2"]
        print(f"  frozen: E_CAS {r['ecas']:.8f}, E(2) {r['e2'] * 1e3:+.4f} mEh, "
              f"<L2> {r['l2']:.4f}, <S2> {r['s2']:.4f} ({r['seconds']:.0f} s) | "
              f"dE_CAS {d['ecas'] * 1e3:+.3f}, dE(2) {d['e2'] * 1e3:+.3f}, "
              f"dE_tot {d['tot'] * 1e3:+.3f} mEh", flush=True)
        rows.append((name, d))
    print("\n===== the screen (mEh against the frozen base) =====", flush=True)
    print(f"  {'channel':8s} {'ncas':>4} {'dE_CAS':>9} {'dE(2)':>9} {'dE_tot':>9}")
    for name, d in rows:
        if d is None:
            print(f"  {name:8s}   refused", flush=True)
        else:
            print(f"  {name:8s} {base.space.ncas + CHANNELS[name].l * 2 + 1:4d} "
                  f"{d['ecas'] * 1e3:+9.3f} {d['e2'] * 1e3:+9.3f} "
                  f"{d['tot'] * 1e3:+9.3f}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
