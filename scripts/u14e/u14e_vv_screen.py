"""The HF-virtual channel screen -- which of the reference's virtual
shells carries correlation?  Printed only (no store writes).

For each declared channel (``vv.channels.CHANNELS``): its radial set
beside the frozen base window, then the frozen CASCI (the free route,
no optimizer) with NEVPT2 on the ground root.  Read per channel
against the base built the same way without channels:

  the extended window's orthonormality and rotation closure; the
  frozen ground's term (<L²>, <S²> -- the base's) and the channel's
  natural occupations;
  dE_CAS   the frozen variational energy's change;
  dE(2)    the perturbative correction's change;
  dE_tot   their sum.

WHAT THE READING MEANS.  A channel is a SPECTATOR when dE_CAS ~ 0: it
does not couple at the frozen level and its row says nothing about
correlation.  A coupling channel's frozen gain is orbital relaxation
of the occupied active orbitals toward the channel (the wavefunction's
rank does not change), which an optimizer recovers on its own: the
screen ranks channels, it does not size the optimized rungs.  A
channel that converts PT2 into the CAS reads dE_CAS < 0, dE(2) > 0,
dE_tot ~ 0; one where PT2 was deficient reads dE_tot < 0.

THE PARTITION.  E(2) is read with the window priced AS SOLVED
(``vv.pricing.partition_as_solved``): PySCF's default prelude
re-partitions empty active orbitals with the virtuals, which would
make every spectator row the same re-partition number (the last
column shows that default reading's shift for comparison).

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
from vv.pricing import occupations, partition_as_solved, price_pyscf  # noqa: E402
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
    """The frozen reading of one window: the ground root's natural
    occupations and its E(2) in both partitions."""
    t0 = time.time()
    solved_casci = free_casci(sub, space)
    state = solved_casci.states[0]
    meters = solved_casci.meters()[0]
    default = price_pyscf(state)
    with partition_as_solved():
        solved = price_pyscf(state)
    return {"ecas": state.e_ref, "e2": solved, "e2_default": default,
            "l2": meters[0], "s2": meters[1], "occ": occupations(state),
            "seconds": time.time() - t0}


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
    n0 = base.space.ncas
    b = frozen(sub, base.space)
    print(f"\n===== the base {args.base}: {n0} active columns; frozen E_CAS "
          f"{b['ecas']:.8f}, E(2) as solved {b['e2'] * 1e3:+.4f} mEh "
          f"(PySCF's default partition: {b['e2_default'] * 1e3:+.4f}); "
          f"<L2> {b['l2']:.4f}, <S2> {b['s2']:.4f} ({b['seconds']:.0f} s) "
          f"=====", flush=True)
    print("  natural occupations: "
          + " ".join(f"{o:.4f}" for o in b["occ"]), flush=True)
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
        d = {"ecas": r["ecas"] - b["ecas"], "e2": r["e2"] - b["e2"],
             "default": r["e2_default"] - b["e2_default"]}
        d["tot"] = d["ecas"] + d["e2"]
        d["note"] = ""
        if d["ecas"] > -5e-7:
            d["note"] = "SPECTATOR: no coupling at the frozen level"
        d["occ"] = r["occ"][-(space.ncas - n0):]
        print(f"  frozen: E_CAS {r['ecas']:.8f}, <L2> {r['l2']:.4f}, "
              f"<S2> {r['s2']:.4f} ({r['seconds']:.0f} s) | dE_CAS "
              f"{d['ecas'] * 1e3:+.3f}, dE(2) {d['e2'] * 1e3:+.3f}, dE_tot "
              f"{d['tot'] * 1e3:+.3f} mEh (as solved) | E(2) shift under "
              f"PySCF's default partition {d['default'] * 1e3:+.3f} mEh | "
              f"the channel's natural occupations "
              f"{d['occ'].min():.1e}..{d['occ'].max():.1e}"
              + (f" | {d['note']}" if d["note"] else ""), flush=True)
        rows.append((name, d))
    print("\n===== the screen (mEh against the frozen base; E(2) with the "
          "window priced as solved) =====", flush=True)
    print(f"  {'channel':8s} {'ncas':>4} {'dE_CAS':>9} {'dE(2)':>9} "
          f"{'dE_tot':>9} {'chan. occ max':>13}")
    for name, d in rows:
        if d is None:
            print(f"  {name:8s}   refused", flush=True)
        else:
            print(f"  {name:8s} {n0 + CHANNELS[name].l * 2 + 1:4d} "
                  f"{d['ecas'] * 1e3:+9.3f} {d['e2'] * 1e3:+9.3f} "
                  f"{d['tot'] * 1e3:+9.3f} {d['occ'].max():13.1e}   "
                  f"{d['note']}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
