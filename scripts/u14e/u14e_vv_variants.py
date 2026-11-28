"""The window-variant brute force -- for each active-space size, which
window of that size is the minimal one?  Printed only (no store
writes).

A size class (N active orbitals for the 4 valence electrons) has many
windows: the BASE plus any subset of the declared HF-virtual channels
(``vv.channels.CHANNELS``).  Every subset is a variant; this script
enumerates them and compares them where the comparison is exact: the
state average of the ground term CONSTRAINED to the window's
composition (``vv.routes.radial_sa``) -- its ensemble energy is
variational, so the lowest is a size's minimal window.  Each variant's
root 0 is priced (NEVPT2 in PySCF's partition; also with the window
priced as solved when the two differ).

THE BASE is the free state average's CONVERGED (4,12)@r1 window, not
the frozen key window.  The frozen window's lowest 17 roots are the
5L term of 5f^2 6d^2 (2.0 electrons in f, 2.0 in d -- observed); the
converged window holds the ground term, 5L of 5f^3 6d^1, as its
lowest level.  A channel is set beside that block, and the base
itself is the sweep's first variant (its reading must reproduce the
free state average's ensemble).

WHAT IS NOT A COMPARISON.
  - The FREE state average does not keep a composition: in C1 it
    rotates empty non-f orbitals into f (observed on the U14e base:
    an s and a p channel converge with f-weight 0.999..1.000, the
    term's roots splitting by 0.4..1.3 mEh).  It is recorded per
    variant as what the table's free route reads there (``--no-free``
    skips it) and never ranks.
  - A frozen CASCI is blind to a spectator and not variational between
    windows (``--stage frozen`` / ``both`` is a record).
  - A variant whose converged roots do not all carry the term
    (<L^2>, <S^2> within ``TERM_DEV``), or whose electron count in
    the base's f block leaves 3, is flagged and not ranked; a state
    average whose tracked roots are not the converged window's lowest
    level is refused by the route itself.
E(2) is not variational: a size whose lowest E_SA is not its lowest
total is flagged.

Each result is ONE record line (``[stage] key=value ...``); ``--from``
skips the variants a log already holds (a sweep resumes), and
``--summary`` ranks the records of any logs, with the table's own
free-SA rows of the same sizes beside them.

Run from the repo root -- one detached process whose ``--jobs``
workers take ``--threads`` each:
    OMP_NUM_THREADS=64 python scripts/u14e/u14e_vv_variants.py \\
        [--min-size 12 --max-size 19] [--jobs 4 --threads 32] \\
        | tee logs/u14e_vv_variants.log
    python scripts/u14e/u14e_vv_variants.py --summary \\
        --from logs/u14e_vv_variants.log
"""

import os
import sys


def _argv_value(flag):
    """``--flag value`` or ``--flag=value`` from the command line."""
    for i, a in enumerate(sys.argv):
        if a == flag and i + 1 < len(sys.argv):
            return sys.argv[i + 1]
        if a.startswith(flag + "="):
            return a.split("=", 1)[1]
    return None


# the BLAS/OMP pools form when numpy loads: a worker's count must be
# in the environment before the first import (a spawned worker re-runs
# this module with the parent's command line)
_THREADS = _argv_value("--threads")
if _THREADS:
    for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                 "MKL_NUM_THREADS"):
        os.environ[_var] = _THREADS

import argparse  # noqa: E402
import itertools  # noqa: E402
import multiprocessing  # noqa: E402
import re  # noqa: E402
import signal  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import numpy as np  # noqa: E402
from pyscf import fci  # noqa: E402

from atomterm.windows import (build_full_mo, occupied_core_idx,  # noqa: E402
                              shell_weights)

from vv import cases as ladder  # noqa: E402
from vv.channels import CHANNELS, channel_block  # noqa: E402
from vv.pricing import (ensemble_occupations, occupations,  # noqa: E402
                        partition_as_solved, price_pyscf, resorted)
from vv.reference import Substrate, rohf_energy  # noqa: E402
from vv.routes import (BUDGET, Route, free_casci, free_sa,  # noqa: E402
                       radial_sa)
from vv.table import load_rows  # noqa: E402
from vv.windows import (ActiveSpace, Built, window_with_channels)  # noqa: E402

BASE_TAG = "(4,12)@r1"
L2_TERM, S2_TERM = ladder.L2_TERM, ladder.S2_QUINTET
TERM_DEV = 1e-4     # the radial route's own band for "every root carries
                    # the term" (routes.SEED_DEV; its closed windows read
                    # <= 3e-06)
CFG_DEV = 0.1       # electrons: the base's f block must hold 3 and its
                    # d block 1 (the ground configuration 5f^3 6d^1)
CLOSURE = 1e-2      # a start window's rotation-closure meters: a construction
                    # check, not a certificate -- the base block carries the
                    # free optimizer's residual asymmetry (observed 1.3e-03);
                    # the lm-pure basis projects it out (observed 1e-15)
TIE = 0.02          # mEh of E_SA: variants this close share the minimum
STAGES = ("frozen", "sa", "both")
DISPLAY = ("5f", "6f", "7d", "8d", "7p", "8p", "8s", "9s")   # a name's channel order
F_BLOCK, D_BLOCK = 7, 5        # the base's columns: one 5f radial, one 6d radial
# the record's keys in print order; a float prints to 12 significant
# digits, a string with blanks quoted
_KEYS = ("name", "ncas", "stage", "done", "error", "conv", "t",
         "e_cas", "fpops", "e_sa", "e0", "spread", "l2_lo", "l2_hi",
         "s2_max", "dev", "pops", "added", "resorted", "e2", "e2s", "gap",
         "closure0", "closure", "ortho", "residual", "core_leak",
         "fconv", "fe", "fdev", "fkept")


def _parse(argv) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="the window-variant brute force on the (4,12)@r1 "
                    "base: the composition-preserving state average of "
                    "the term per variant, ranked per size")
    p.add_argument("--base", default=BASE_TAG, metavar="TAG",
                   help=f"the base case (default: {BASE_TAG})")
    p.add_argument("--channels", nargs="+", default=list(CHANNELS),
                   metavar="SHELL", help="the channels to combine "
                   f"(default: {list(CHANNELS)})")
    p.add_argument("--min-size", type=int, default=12, metavar="N",
                   help="the smallest active-space size (default: 12: the base)")
    p.add_argument("--max-size", type=int, default=19, metavar="N",
                   help="the largest active-space size (default: 19)")
    p.add_argument("--sizes", nargs="+", type=int, default=[],
                   metavar="N", help="only these sizes")
    p.add_argument("--variants", nargs="+", default=[], metavar="NAME",
                   help="only these variants (e.g. 5f+8s; the base is 'base')")
    p.add_argument("--stage", choices=STAGES, default="sa",
                   help="frozen CASCI, the state average, or both "
                        "(default: sa)")
    p.add_argument("--no-e2", action="store_true",
                   help="no NEVPT2 (the state averages and the meters only)")
    p.add_argument("--no-free", action="store_true",
                   help="skip the free state average beside the radial one")
    p.add_argument("--max-cycles", type=int, default=0, metavar="K",
                   help="cap the optimizers' macro cycles (default: the "
                        "library's budget, 200)")
    p.add_argument("--jobs", type=int, default=1, metavar="J",
                   help="concurrent worker processes (default: 1)")
    p.add_argument("--threads", default=None, metavar="T",
                   help="OMP/BLAS threads per process (default: the "
                        "environment's)")
    p.add_argument("--from", dest="resume", nargs="*", default=[],
                   metavar="LOG", help="logs whose recorded variants "
                        "are skipped")
    p.add_argument("--summary", nargs="*", default=None, metavar="LOG",
                   help="rank the records of these logs (plus --from) "
                        "and exit")
    p.add_argument("--list", action="store_true",
                   help="print the variants per size and exit")
    a = p.parse_args(argv)
    unknown = [c for c in a.channels if c not in CHANNELS]
    if unknown:
        p.error(f"unknown channel(s) {unknown}; declared: {list(CHANNELS)}")
    a.resume = list(a.resume)
    return a


# ----------------------------------------------------------- the enumeration

def name_of(combo) -> str:
    """A variant's name: its channels in ``DISPLAY`` order, joined by
    ``+`` (the base: ``base``)."""
    return "+".join(sorted(combo, key=lambda c: DISPLAY.index(c))) or "base"


def width_of(combo) -> int:
    """The orbitals a subset adds to the base."""
    return sum(CHANNELS[c].l * 2 + 1 for c in combo)


def variants(names=None) -> dict[int, list[tuple[str, ...]]]:
    """Every subset of the channels by the width it adds to the base
    (the empty subset included)."""
    names = list(names if names is not None else CHANNELS)
    out = {}
    for k in range(len(names) + 1):
        for combo in itertools.combinations(names, k):
            out.setdefault(width_of(combo), []).append(combo)
    return out


def _selected(args, base_width, done=()) -> list[tuple]:
    """The (name, combo, ncas) to run, by size then declaration order;
    the ones in ``done`` are skipped."""
    out = []
    table = variants(args.channels)
    for width in sorted(table):
        n = base_width + width
        if not args.min_size <= n <= args.max_size:
            continue
        if args.sizes and n not in args.sizes:
            continue
        for combo in table[width]:
            name = name_of(combo)
            if args.variants and name not in args.variants:
                continue
            if name not in done:
                out.append((name, combo, n))
    return out


# ------------------------------------------------------------- the readings

_SUB = _BASE = _BLOCKS = None      # the worker's substrate, base window, channel blocks


def _start_worker(base_tag):
    global _SUB, _BASE, _BLOCKS
    _SUB = Substrate.load()
    _BASE = _converged_base(_SUB, base_tag)
    _BLOCKS = {n: channel_block(_SUB, CHANNELS[n]) for n in CHANNELS}


def _converged_base(sub, base_tag) -> Built:
    """The free state average's converged ACTIVE block, assembled on
    the key window's purified core -- the ground term's basin (the
    module docstring's THE BASE) on an exactly rotation-closed core:
    the free optimizer also relaxes the core off that span, and the
    relaxed block splits the term's components by 3.4e-06 Eh where the
    purified one splits them by 6e-08 (the ensemble agrees to 1e-08).
    """
    frozen = ladder.CASES[base_tag].window(sub)
    sol = free_sa(sub, frozen.space, ladder.NROOTS)
    st = sol.states[0]
    act = np.array(st.active)
    mo = build_full_mo(sub.mf, occupied_core_idx(sub.mf, st.ncore), act,
                       c_core=frozen.space.core)
    return Built(ActiveSpace(mo, st.ncore, st.ncas, st.nelecas),
                 (f"base: the free state average's converged active block "
                  f"of {base_tag} on its purified core (ensemble "
                  f"{sol.ensemble:.8f}, {sol.conv})",))


def _window(combo):
    return window_with_channels(_SUB, _BASE, [_BLOCKS[c] for c in combo]).space


def _edges(combo, base_width) -> list[int]:
    """The block boundaries of a variant's active columns: the base's f
    block, its d block, then each channel."""
    edges = [0, F_BLOCK, base_width]
    for c in combo:
        edges.append(edges[-1] + CHANNELS[c].l * 2 + 1)
    return edges


def _pops(sol, edges) -> list[float]:
    """The ensemble's electrons per block (the mean 1-RDM's diagonal,
    summed block by block)."""
    dm = sum(fci.direct_spin1.make_rdm1(st.ci, st.ncas, st.nelecas)
             for st in sol.states) / len(sol.states)
    d = np.real(np.diag(dm))
    return [float(d[a:b].sum()) for a, b in zip(edges[:-1], edges[1:])]


def _added(occ, width) -> float:
    """The largest of the ``width`` lowest natural occupations -- how
    much the added span is used (NaN for the base)."""
    return float(occ[-width:].max()) if width else np.nan


def _dev(l2, s2) -> float:
    """The roots' largest distance from the term's (<L^2>, <S^2>)."""
    return max(max(abs(x - L2_TERM) for x in l2),
               max(abs(x - S2_TERM) for x in s2))


def _sa_reading(sol, edges, width, t_sa) -> dict:
    """What both state averages report: energies, the term's meters
    over the roots, the ensemble's block populations."""
    meters = sol.meters()
    l2, s2 = [a for a, _b in meters], [b for _a, b in meters]
    roots = [st.e_ref for st in sol.states]
    pops = _pops(sol, edges)
    return {"e_sa": float(sol.ensemble), "e0": roots[0],
            "spread": (max(roots) - min(roots)) * 1e3,
            "l2_lo": min(l2), "l2_hi": max(l2), "s2_max": max(s2),
            "dev": _dev(l2, s2), "conv": sol.conv, "t": t_sa,
            "added": _added(ensemble_occupations(sol.states), width),
            "pops": "/".join(f"{x:.3f}" for x in pops),
            "nf": pops[0], "nd": pops[1], "nadd": sum(pops[2:])}


def _frozen(space, combo, base_width) -> dict:
    """The frozen CASCI in the window as built: the ground's energy,
    term and block populations."""
    sol = free_casci(_SUB, space)
    (l2, s2), = sol.meters()
    pops = _pops(sol, _edges(combo, base_width))
    return {"e_cas": sol.states[0].e_ref, "l2_lo": l2, "l2_hi": l2,
            "s2_max": s2, "dev": _dev([l2], [s2]),
            "fpops": "/".join(f"{x:.3f}" for x in pops),
            "added": _added(occupations(sol.states[0]),
                            space.ncas - base_width),
            "resorted": resorted(sol.states[0])}


def _radial(space, combo, base_width, budget, price: bool) -> dict:
    """The state average constrained to the window's composition, with
    root 0's E(2)."""
    t0 = time.time()
    sol = radial_sa(_SUB, space, (L2_TERM, S2_TERM), ladder.NROOTS, budget)
    t_sa = time.time() - t0
    rec = _sa_reading(sol, _edges(combo, base_width),
                      space.ncas - base_width, t_sa)
    rec.update({"resorted": 0, "e2": np.nan, "e2s": np.nan,
                "gap": float(sol.certificate["level_gap"]),
                "closure": float(sol.certificate["closure"]),
                "ortho": sol.certificate["ortho"],
                "residual": sol.certificate["residual"],
                "core_leak": sol.certificate["core_leak"]})
    if price:
        st0 = sol.states[0]
        rec["e2"] = price_pyscf(st0) * 1e3
        rec["resorted"] = resorted(st0)
        if rec["resorted"]:            # the two partitions coincide otherwise
            with partition_as_solved():
                rec["e2s"] = price_pyscf(st0) * 1e3
        else:
            rec["e2s"] = rec["e2"]
    return rec


def _free(space, combo, base_width, budget) -> dict:
    """The table's own route on the same window: what the free state
    average makes of the composition (``fkept``: the share of the
    added columns still dominated by their channel's l)."""
    t0 = time.time()
    sol = free_sa(_SUB, space, ladder.NROOTS, budget)
    rec = _sa_reading(sol, _edges(combo, base_width),
                      space.ncas - base_width, time.time() - t0)
    st = sol.states[0]
    want = [CHANNELS[c].l for c in combo for _ in range(CHANNELS[c].l * 2 + 1)]
    if want:
        _, w = shell_weights(st.mol, st.active[:, base_width:], st.mf.get_ovlp())
        kept = sum(int(np.argmax(w[k])) == l for k, l in enumerate(want))
    else:
        kept = 0
    return {"fconv": rec["conv"], "fe": rec["e_sa"], "fdev": rec["dev"],
            "fkept": (kept / len(want)) if want else np.nan}


def _evaluate(task) -> dict:
    """One variant, one stage: its record (a refused or failed variant
    is recorded with its error; the sweep goes on)."""
    name, combo, ncas, stage, price, free, max_cycles, base_width = task
    t0 = time.time()
    try:
        space = _window(combo)
        closure0 = max(space.closure(_SUB.mol).values())
        if closure0 > CLOSURE:
            raise ValueError(f"the start window is not rotation-closed "
                             f"(max meter {closure0:.1e})")
        rec = {"ncas": space.ncas, "closure0": closure0}
        if stage in ("frozen", "both"):
            rec.update(_frozen(space, combo, base_width))
        if stage in ("sa", "both"):
            budget = ({**BUDGET, "max_cycle_macro": max_cycles}
                      if max_cycles else BUDGET)
            rec.update(_radial(space, combo, base_width, budget, price))
            if free:
                rec.update(_free(space, combo, base_width, budget))
        rec["done"] = 1
    except Exception as e:                    # noqa: BLE001 -- recorded, not raised
        text = " ".join(str(e).split())
        rec = {"ncas": ncas, "done": 0,
               "error": f"{type(e).__name__}: {text}"[:300]}
    rec["name"] = name
    rec["stage"] = stage
    rec.setdefault("t", time.time() - t0)
    return rec


# ------------------------------------------------------------- the records

def fmt_record(rec) -> str:
    """One record as one ``[stage] key=value`` line."""
    parts = []
    for k in _KEYS:
        if k not in rec or rec[k] is None:
            continue
        v = rec[k]
        if isinstance(v, float):
            if np.isnan(v):
                continue
            parts.append(f"{k}={v:.12g}")
        elif isinstance(v, str) and (" " in v or k in ("pops", "fpops")):
            parts.append(f'{k}="{v}"')
        else:
            parts.append(f"{k}={v}")
    return f"  [{rec['stage']}] " + " ".join(parts)


def _unquote(text):
    if text.startswith('"') and text.endswith('"'):
        return text[1:-1]
    if re.fullmatch(r"[-+]?\d+", text):
        return int(text)
    try:
        return float(text)
    except ValueError:
        return text


def read_records(logs) -> dict:
    """The LATEST record of every (stage, variant) in the logs, each
    carrying its position ``seq`` in the logs' order."""
    out, seq = {}, 0
    for log in logs:
        path = Path(log)
        if not path.exists():
            continue
        for line in path.read_text(errors="replace").splitlines():
            m = re.match(r"\s*\[(frozen|sa|both)\] (.*)", line)
            if not m:
                continue
            rec = {k: _unquote(v) for k, v in
                   re.findall(r'(\w+)=("[^"]*"|\S+)', m.group(2))}
            if "name" not in rec:
                continue
            rec["stage"], rec["seq"] = m.group(1), seq
            seq += 1
            out[(m.group(1), rec["name"])] = rec
    return out


def _l2_cell(r) -> str:
    if r["l2_hi"] - r["l2_lo"] < 5e-7:
        return f"{r['l2_lo']:.6f}"
    return f"{r['l2_lo']:.6f}..{r['l2_hi']:.6f}"


def _from_row(r) -> dict:
    """A table row (the free SA's, root 0) as a listing entry."""
    n = int(re.match(r"\(4,(\d+)\)", r.cas).group(1))
    lo, _, hi = r.l2.partition("..")
    l2_lo, l2_hi = float(lo), float(hi or lo)
    return {"name": f"table {r.cas}", "ncas": n, "seq": -1, "conv": r.conv,
            "e_sa": r.e_method, "e0": r.e_method, "spread": np.nan,
            "l2_lo": l2_lo, "l2_hi": l2_hi, "s2_max": float(r.s2),
            "dev": _dev([l2_lo, l2_hi], [float(r.s2)]),
            "added": np.nan, "resorted": 0, "e2": r.e2 * 1e3, "e2s": r.e2 * 1e3,
            "nf": np.nan, "nd": np.nan, "nadd": np.nan, "gap": np.nan,
            "fe": np.nan, "fdev": np.nan, "fkept": np.nan}


def _rankable(r) -> bool:
    """A variant ranks when its roots carry the term and its ensemble
    is the ground configuration's."""
    return (r["dev"] <= TERM_DEV and abs(r["nf"] - 3.0) <= CFG_DEV
            and abs(r["nd"] - 1.0) <= CFG_DEV)


def _note(r, best) -> tuple[str, bool]:
    """(the note, whether the variant is the size's minimum)."""
    if r["name"].startswith("table "):
        return "the table's free row", False
    if r["dev"] > TERM_DEV:
        return f"NOT the term (max dev {r['dev']:.1e})", False
    if abs(r["nf"] - 3.0) > CFG_DEV or abs(r["nd"] - 1.0) > CFG_DEV:
        return (f"NOT the ground configuration (f {r['nf']:.2f}, "
                f"d {r['nd']:.2f})"), False
    if r["conv"].endswith("*"):
        return "cycle cap: not converged", False
    unpriced = "" if not np.isnan(r.get("e2", np.nan)) else " (not priced)"
    if r["e_sa"] - best <= TIE * 1e-3:
        return "MINIMAL" + unpriced, True
    return f"+{(r['e_sa'] - best) * 1e3:.3f} mEh" + unpriced, False


def summary(logs, sizes=()) -> None:
    """Rank every size's variants by the constrained state average's
    E_SA and print the minimal ladder, with the table's own free-SA
    rows of the same sizes beside them."""
    e_rohf = rohf_energy()
    latest = {}
    for (stage, name), r in read_records(logs).items():
        if stage in ("sa", "both") and r.get("done") == 1:
            if name not in latest or r["seq"] > latest[name]["seq"]:
                blocks = [float(x) for x in r["pops"].split("/")]
                r.update(nf=blocks[0], nd=blocks[1], nadd=sum(blocks[2:]))
                latest[name] = r
    by_size = {}
    for r in latest.values():
        by_size.setdefault(r["ncas"], []).append(r)
    for key, r in load_rows().items():
        if (re.match(r"\(4,\d+\)", r.cas) and r.method == "sa-casscf"
                and r.route is Route.FREE and not r.withdrawn):
            by_size.setdefault(_from_row(r)["ncas"], []).append(_from_row(r))
    print("\n===== the variants by size =====", flush=True)
    print("  E_SA = the composition-preserving state average's ensemble energy (the\n"
          "  variational comparison); static / dynamic / total = root 0 against ROHF,\n"
          "  mEh (dynamic = E(2) in PySCF's partition; dyn_s / tot_s price the window\n"
          "  as solved); e(f/d/add) = the ensemble's electrons in the base's f block,\n"
          "  its d block, the added span; gap = the converged window's next level;\n"
          "  free = the table's free route on the same window: E_SA, its largest term\n"
          f"  deviation, the share of added columns keeping their l.  MINIMAL = the\n  lowest E_SA among the term-preserving, ground-configuration variants (ties within {TIE:g} mEh).",
          flush=True)
    ladder_rows = []
    for n in sorted(by_size):
        if sizes and n not in sizes:
            continue
        rows = sorted(by_size[n], key=lambda r: r["e_sa"])
        clean = [r for r in rows if _rankable(r)]
        best = min((r["e_sa"] for r in clean), default=None)
        state = "" if clean else "; NONE preserves the term and the ground configuration"
        print(f"\n  --- size {n} ({len(rows)} window(s){state}) ---", flush=True)
        print(f"  {'variant':16s} {'conv':>12s} {'E_SA':>14s} {'static':>9s} "
              f"{'dynamic':>9s} {'total':>9s} {'dyn_s':>9s} {'tot_s':>9s} "
              f"{'e(f/d/add)':>19s} {'gap':>7s} | {'free E_SA':>14s} {'dev':>8s} "
              f"{'kept':>5s}  note", flush=True)
        priced = []
        for r in rows:
            e2 = r.get("e2", np.nan)
            static = (r["e0"] - e_rohf) * 1e3
            tot = static + e2 if not np.isnan(e2) else np.nan
            e2s = r.get("e2s", e2)
            note, minimal = _note(r, best)
            if minimal and n not in [x[0] for x in ladder_rows]:
                ladder_rows.append((n, r, static, e2, tot))
            if r in clean and not np.isnan(e2):
                priced.append((r, tot))
            solved = "" if np.isnan(e2) else f"{e2s:+9.3f} {static + e2s:+9.3f}"
            cells = "  ".join("-" if np.isnan(x) else f"{x:.2f}"
                              for x in (r.get("nf", np.nan), r.get("nd", np.nan),
                                        r.get("nadd", np.nan)))
            gap = "-" if np.isnan(r.get("gap", np.nan)) else f"{r['gap']:7.2f}"
            free = ("-" * 14 + " " + "-" * 8 + " " + "-" * 5
                    if np.isnan(r.get("fe", np.nan))
                    else f"{r['fe']:14.8f} {r['fdev']:8.1e} "
                         f"{r.get('fkept', np.nan):5.2f}")
            print(f"  {r['name']:16s} {r['conv']:>12s} {r['e_sa']:14.8f} "
                  f"{static:+9.3f} {e2:+9.3f} {tot:+9.3f} {solved} "
                  f"{cells:>19s} {gap} | {free}  {note}", flush=True)
        if priced:
            r_tot, e_tot = min(priced, key=lambda x: x[1])
            r_sa = min(priced, key=lambda x: x[0]["e_sa"])[0]
            if r_tot is not r_sa and r_tot["e_sa"] - r_sa["e_sa"] > TIE * 1e-3:
                print(f"  NOTE: the lowest total is {r_tot['name']} "
                      f"({e_tot:+.3f} mEh), not the lowest E_SA", flush=True)
    print("\n===== the minimal ladder (increments against the line above) =====",
          flush=True)
    print(f"  {'size':>4s} {'variant':16s} {'conv':>12s} {'static':>9s} "
          f"{'dynamic':>9s} {'total':>9s} | {'d_static':>9s} {'d_dynamic':>9s} "
          f"{'d_total':>9s}", flush=True)
    prev = None
    for n, r, static, e2, tot in ladder_rows:
        inc = "" if prev is None else (f"| {static - prev[0]:+9.3f} "
                                       f"{cell(e2 - prev[1], '+9.3f')} "
                                       f"{cell(tot - prev[2], '+9.3f')}")
        def cell(x, fmt):
            return "-" * 9 if np.isnan(x) else format(x, fmt)
        print(f"  {n:4d} {r['name']:16s} {r['conv']:>12s} {static:+9.3f} "
              f"{cell(e2, '+9.3f')} {cell(tot, '+9.3f')} {inc}", flush=True)
        prev = (static, e2, tot)


# -------------------------------------------------------------------- main

def main(argv=None) -> int:
    args = _parse(sys.argv[1:] if argv is None else argv)
    if args.summary is not None:            # present (maybe with no log)
        summary(list(args.summary) + list(args.resume), args.sizes)
        return 0
    base_width = ladder.CASES[args.base].window(Substrate.load()).space.ncas
    if args.list:
        table = variants(args.channels)
        for width in sorted(table):
            n = base_width + width
            inside = args.min_size <= n <= args.max_size
            print(f"  size {n:2d}: {len(table[width]):2d} variant(s): "
                  + "  ".join(name_of(c) for c in table[width])
                  + ("" if inside else "   (outside the size range)"))
        return 0
    done = {name for stage, name in read_records(args.resume)
            if stage == args.stage or stage == "both"}
    tasks = [(name, combo, n, args.stage, not args.no_e2, not args.no_free,
              args.max_cycles, base_width)
             for name, combo, n in _selected(args, base_width, done)]
    if not tasks:
        print("nothing to run: every selected variant is recorded", flush=True)
        return 0
    by_n = {}
    for _name, _combo, n, *_rest in tasks:
        by_n[n] = by_n.get(n, 0) + 1
    print(f"===== the variant sweep: {len(tasks)} variant(s) on the free state "
          f"average's converged {args.base} base ({base_width} active "
          f"columns); stage {args.stage}; {args.jobs} job(s) x "
          f"{os.environ.get('OMP_NUM_THREADS', '?')} threads =====", flush=True)
    print("  per size: " + "  ".join(f"{n}: {k}" for n, k in sorted(by_n.items())),
          flush=True)
    t0 = time.time()
    if args.jobs == 1:
        _start_worker(args.base)
        for task in tasks:
            print(fmt_record(_evaluate(task)), flush=True)
    else:
        signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))
        pool = multiprocessing.get_context("spawn").Pool(
            args.jobs, initializer=_start_worker, initargs=(args.base,))
        try:
            for rec in pool.imap_unordered(_evaluate, tasks):
                print(fmt_record(rec), flush=True)
            pool.close()
        except BaseException:                  # a killed sweep takes its workers
            pool.terminate()
            raise
        finally:
            pool.join()
    print(f"===== the sweep ended: {len(tasks)} variant(s) in "
          f"{time.time() - t0:.0f} s =====", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
