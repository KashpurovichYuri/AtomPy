"""The E(2) class decomposition of the composition-preserving state
average -- why does a NEAR-EMPTY second radial remove dynamical
correlation?  Printed only (no store writes).

THE QUESTION.  In the window-variant sweep (u14e_vv_variants.py) a
second 6d radial beside the (4,12)@r1 base lowers the state-averaged
energy by 1.9 mEh and raises E(2) by 25 mEh; a second 5f radial: 5.9
down, 12 up.  The added span holds 0.002..0.005 electrons.  Which
Koopmans class of the strongly contracted decomposition carries the
change -- and is it the window's (the classes' norms and energies on
the solved window) or the pricing partition's (the same state priced
as solved)?

Per variant: the constrained state average (``vv.routes.radial_sa``;
the base is the sweep's), then root 0 priced through the class
decomposition (``vv.pricing.price_classes``) in PySCF's partition and
in the solved one (``partition_as_solved``).

Run from the repo root:
    OMP_NUM_THREADS=32 python scripts/probes/probe_u14e_e2_radials.py \\
        | tee logs/probe_u14e_e2_radials.log
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "u14e"))

import warnings  # noqa: E402

warnings.simplefilter("ignore")

import u14e_vv_variants as sweep  # noqa: E402

from vv import cases as ladder  # noqa: E402
from vv.pricing import partition_as_solved, price_classes  # noqa: E402
from vv.reference import Substrate, rohf_energy  # noqa: E402
from vv.routes import radial_sa  # noqa: E402

VARIANTS = ((), ("7p",), ("7d",), ("5f",))
E_ROHF = rohf_energy()


def main() -> int:
    sub = Substrate.load()
    sweep._start_worker(sweep.BASE_TAG)
    rows = {}
    for combo in VARIANTS:
        name = sweep.name_of(combo)
        t0 = time.time()
        sol = radial_sa(sub, sweep._window(combo),
                        (sweep.L2_TERM, sweep.S2_TERM), ladder.NROOTS)
        st0 = sol.states[0]
        t_sa = time.time() - t0
        t0 = time.time()
        default = price_classes(st0)
        with partition_as_solved():
            solved = price_classes(st0)
        t_e2 = time.time() - t0
        rows[name] = (sol, default, solved)
        print(f"===== {name} (ncas {sol.states[0].ncas}): radial SA {sol.ensemble:.8f} "
              f"({sol.conv}, {t_sa:.0f} s); E(2) default {default.e_corr * 1e3:+.4f} mEh, "
              f"as solved {solved.e_corr * 1e3:+.4f} mEh ({t_e2:.0f} s) =====",
              flush=True)
    base = rows["base"]
    print("\n===== the Koopmans classes of root 0 (mEh; norm in brackets); "
          "the base beside each variant =====", flush=True)
    for name in rows:
        sol, default, solved = rows[name]
        static = (sol.states[0].e_ref - E_ROHF) * 1e3
        print(f"\n  --- {name}: static {static:+9.3f} | dynamic {default.e_corr * 1e3:+9.3f} "
              f"(as solved {solved.e_corr * 1e3:+9.3f}) | total {static + default.e_corr * 1e3:+9.3f} "
              f"(as solved {static + solved.e_corr * 1e3:+9.3f}) ---", flush=True)
        print("     class        note        norm      E(2)        d vs base     E(2) as solved   d vs base",
              flush=True)
        for c, c0, cs, cs0 in zip(default.classes, base[1].classes,
                                  solved.classes, base[2].classes):
            print(f"     {c.label:10s} {c.note:10s} {c.norm:9.4f} {c.energy * 1e3:+10.4f} "
                  f"{(c.energy - c0.energy) * 1e3:+11.4f}     {cs.energy * 1e3:+11.4f} "
                  f"{(cs.energy - cs0.energy) * 1e3:+11.4f}", flush=True)
        print(f"     {'sum':10s} {'':10s} {'':9s} {default.e_corr * 1e3:+10.4f} "
              f"{(default.e_corr - base[1].e_corr) * 1e3:+11.4f}     {solved.e_corr * 1e3:+11.4f} "
              f"{(solved.e_corr - base[2].e_corr) * 1e3:+11.4f}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
