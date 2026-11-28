"""The CAS-ladder leg -- THE WRITER of the table's method rows.

Every case of ``vv.cases`` goes through the routes the table keeps:
the free trio (PySCF-direct, on the purified core), the pinned state
average (the production wiring, on the same window and the same core)
and -- on the minimal family -- the constrained state average
(``vv.routes.radial_sa``: the window's composition held, the core
frozen; priced with the window AS SOLVED, where the free rows carry
PySCF's partition).  Each solution is priced (strongly contracted
NEVPT2), metered (<L²>, <S²>), checked and written.  The pinned route
first runs its frozen CI (print-only: the gate's reading and the
ledger's frozen point -- the table keeps only the pinned SA).

WHAT GATES (constructions and declared records, never a physics
outcome): each CONSTRUCTED window is S-orthonormal and
rotation-closed; every solution passes its birth certificate
(orbitals orthonormal; CI stationary in its own orbitals --
``vv.routes.assert_consistent``); the frozen ground carries a
definite L; and each case's declared expectations (``vv.cases``): the
frozen term, the free SA's term identity over ALL its roots, the
ground configuration where declared (the electrons of the base's f
and d blocks -- <L²> and <S²> cannot tell one configuration of a term
from another), the committed E(2) anchors, the pin-on-pur gate at
(4,12)@r1, the certified component class at (6,13), the s channel's
pair at (6,13)@r1).  The constrained route refuses at birth (its
seed, its trap check).

WHAT IS DISCLOSED: the per-root tables, the E(2) mean over the roots
that are priced, the free-vs-pinned energy difference, each
solution's core leak, what the free state average made of a
minimal-family window's composition (the share of added columns that
kept their l), and the static / E(2) / total deltas of every row --
in mEh and in eV.

A window its certificate refuses is logged ``[refused]`` and skipped
(the declared refusal class); nothing is stored.

The single-root ``casscf`` leaves the term (its optimum is not a term
component); on a minimal-family window it also runs to the cycle cap
(observed 200 cycles, <L²> 64.7), so those cases declare a 50-cycle
budget (``Case.casscf_cycles``; the conv cell carries the ``*``).

Run from the repo root -- one process, which alone writes the store
(``--jobs J``: its workers compute the cases; the cases go largest
first):
    OMP_NUM_THREADS=64 python scripts/u14e/u14e_vv_ladder.py \\
        [--cases "(4,12)" "(6,13)"] [--routes casci,casscf,sa,pinned,radial] \\
        [--ground-only] [--base "(4,12)m"] [--jobs 4 --threads 32] \\
        |  tee logs/u14e_vv_ladder.log
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
import contextlib  # noqa: E402
import io  # noqa: E402
import multiprocessing  # noqa: E402
import re  # noqa: E402
import signal  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from pyscf import ao2mo, fci  # noqa: E402
from pyscf.mcscf import addons as mcaddons  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from atomterm.irreps import orb_irreps  # noqa: E402
from atomterm.windows import shell_weights  # noqa: E402

from vv import cases as ladder  # noqa: E402
from vv.channels import CHANNELS  # noqa: E402
from vv.pricing import partition_as_solved  # noqa: E402
from vv.reference import Substrate, anchor, rohf_energy  # noqa: E402
from vv.routes import (BUDGET, free_casci, free_casscf,  # noqa: E402
                       free_sa, pinned_casci, pinned_sa, radial_sa)
from vv.symmetry import symmetrize_columns  # noqa: E402
from vv.table import EH_EV, Route, drop_row, load_rows  # noqa: E402
from vv.verdict import Verdict  # noqa: E402
from vv.windows import s_column  # noqa: E402

ROUTES = ("casci", "casscf", "sa", "pinned", "radial")

# the expectations' bands (declared constants; each cites the observed
# value -- the family's committed logs)
BAND = dict(
    l2=1e-5,        # the Casimir meters; a seeded readout's observed
                    # deviation 3.0e-06 (free) and 1.9e-07 (pinned)
    definite=1e-4,  # a frozen ground's distance to the nearest L(L+1)
                    # (observed <= 1e-06 on every closed window)
    e2_meh=2.0,     # a fresh SA draw's root-0 E(2) against its
                    # committed anchor; observed delta 0.0006-0.0029
    routes=1e-6,    # Eh: the free and the pinned SA energies on one
                    # window where the gate declares them equal
                    # (observed 3.3e-09 at (4,12))
    cls=0.01,       # a certified component class's distance to L(L+1)
    closure=1e-10,  # a constructed window's rotation-closure meters
                    # (observed ~1e-13)
    pair=1.5,       # the s channel's occupation in the frozen ground
                    # (the pair forms at ~1.99; the as-built read 3e-03)
    spin=1e-6,      # Eh-free: PySCF's spin_square against the meters
                    # (observed <= 1e-8)
    wrapper=1e-9,   # Eh: the window's ERIs on the D2h wrapper against
                    # the C1 molecule's (observed 0.0)
    pure=1e-6,      # the active columns' irrep impurity on the wrapper
                    # (the labeling seam refuses past it)
    config=0.1,     # electrons: the ensemble's count in the base's f
                    # block against 3 and in its d block against 1
                    # (observed 2.99..3.00 / 1.00)
)

_SUB = None         # the process's substrate (the parent's, a worker's)


def _parse(argv) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="the U14e V&V ladder: every case through the "
                    "table's routes")
    p.add_argument("--cases", nargs="+", default=[], metavar="TAG",
                   help="the case tags to run (default: the whole ladder)")
    p.add_argument("--routes", default=",".join(ROUTES),
                   help="comma-separated subset of "
                        f"{','.join(ROUTES)} (default: all)")
    p.add_argument("--ground-only", action="store_true",
                   help="price root 0 only on the state averages")
    p.add_argument("--base", default=None,
                   help="a case tag: the ledger prints every SA row's "
                        "increments against it, route by route")
    p.add_argument("--jobs", type=int, default=1, metavar="J",
                   help="worker processes computing the cases (default: "
                        "1); this process alone writes the store")
    p.add_argument("--threads", default=None, metavar="T",
                   help="OMP/BLAS threads per process (default: the "
                        "environment's)")
    a = p.parse_args(argv)
    bad = [r for r in a.routes.split(",") if r not in ROUTES]
    if bad:
        p.error(f"unknown route(s) {bad}; choose from {ROUTES}")
    a.route_list = [r for r in a.routes.split(",") if r]
    a.case_list = list(a.cases)
    return a


def _delta(de_eh: float) -> str:
    """One delta in both units: ``+x.xx mEh / +y.yyyy eV``."""
    meh = de_eh * 1e3
    return f"{meh:+.2f} mEh / {meh * EH_EV / 1e3:+.4f} eV"


def _print_solution(name, sol, priced, seconds, as_solved=False) -> None:
    c = sol.certificate
    if priced.resorted:
        if as_solved:
            print(f"  [{name}] NOTE: {priced.resorted} empty active "
                  f"orbital(s); this row is priced with the window AS SOLVED "
                  f"(PySCF's default prelude would re-partition them with "
                  f"the virtuals)", flush=True)
        else:
            print(f"  [{name}] NOTE: {priced.resorted} empty active orbital(s) "
                  f"re-partitioned with the virtuals before pricing (the "
                  f"lowest virtuals take their slots) -- the priced partition "
                  f"is not the solved window", flush=True)
    print(f"  [{name}] {seconds:.1f} s | conv {sol.conv} | certificate: "
          f"orbitals {c['ortho']:.1e}, eigenpairs {c['residual']:.1e}, "
          f"core leak {c['core_leak']:.1e}"
          + (f", readout shift {c['readout_shift']:.1e}"
             if "readout_shift" in c else "")
          + (f", orbitals as returned {c['ortho_in']:.1e}"
             if "ortho_in" in c else "")
          + (f", component mean vs ensemble {c['ensemble_dev']:.1e}"
             if "ensemble_dev" in c else "")
          + (f", next level {c['level_gap']:.2f} mEh"
             if "level_gap" in c else "")
          + (f", window closure {c['closure']:.1e}"
             if "closure" in c else ""), flush=True)
    for s in priced.states:
        block = f" {s.state.irrep}" if s.state.irrep else ""
        print(f"    root {s.state.root}{block}: e_ref {s.state.e_ref:.8f}, "
              f"E(2) {s.e2 * 1e3:+.4f} mEh, <L2> {s.l2:.6f}, "
              f"<S2> {s.s2:.6f}", flush=True)
    if len(priced.states) > 1:
        e2 = np.array(priced.e2) * 1e3
        print(f"    E(2) over {len(e2)} roots: mean {e2.mean():+.4f} mEh, "
              f"spread {np.ptp(e2):.4f} mEh", flush=True)
    if sol.ensemble is not None:
        print(f"    ensemble {sol.ensemble:.8f}", flush=True)


def _wrapper(v, sub, space) -> None:
    """The pinned route's substrate IS the C1 molecule: the two-electron
    integrals over its core + window block are the C1 molecule's (the
    overlap and the one-electron integrals are certified where the
    wrapper is built).  No free route can test this -- on the
    symmetry-tagged substrate PySCF pins the FCI to one irrep."""
    mol, _mf = sub.d2h
    c = space.mo[:, :space.ncore + space.ncas]
    a = ao2mo.restore(1, ao2mo.kernel(sub.mol, c, verbose=0), c.shape[1])
    b = ao2mo.restore(1, ao2mo.kernel(mol, c, verbose=0), c.shape[1])
    dev = float(np.abs(a - b).max())
    v.check("the D2h wrapper is the C1 molecule: the window's "
            "two-electron integrals coincide", dev <= BAND["wrapper"],
            f"{c.shape[1]} orbitals, max|<pq|rs>_C1 - <pq|rs>_wrapper| "
            f"{dev:.1e} (band {BAND['wrapper']:g})")


def _purity(v, sub, space, tag) -> None:
    """The active columns' irrep purity on the wrapper -- the pinned
    solvers' labeling reads it (``orb_irreps`` refuses an impure
    block)."""
    mol, mf = sub.d2h
    s = mf.get_ovlp()
    _, imp = symmetrize_columns(mol, s, space.active, tol=1.0)
    orb_irreps(mol, space.active)
    v.note(f"{tag}: the window's irrep purity on the wrapper",
           f"max impurity {max(imp):.2e} (the labeling seam's band "
           f"{BAND['pure']:g}); labels {orb_irreps(mol, space.active).tolist()}")


def _run(name, go, roots=None, as_solved=False):
    """One route: its solution priced (``as_solved``: with the window
    priced as solved -- the constrained route's rows)."""
    t0 = time.time()
    sol = go()
    if as_solved:
        with partition_as_solved():
            priced = sol.price(roots=roots)
    else:
        priced = sol.price(roots=roots)
    _print_solution(name, sol, priced, time.time() - t0, as_solved)
    return priced


def _window(v, case, built, sub) -> None:
    for line in built.report:
        print(f"  {line}", flush=True)
    space = built.space
    dev = space.certify(sub.mf)
    cl = space.closure(sub.mol)
    if not case.start_closed:
        # an optimizer's converged block carries its residual asymmetry;
        # the constrained route's lm-pure basis projects it out
        v.note(f"{case.tag}: the start window's rotation closure",
               f"{space.ncas} active columns on {space.ncore} core pairs, "
               f"sector {space.sector}; max|C^T S C - I| {dev:.1e}; "
               f"||[P,Lz]|| {cl['lz']:.2e}, ||[P,Lperp]|| {cl['lperp']:.2e} "
               f"(the free optimizer's residual; not a constructed window)")
        return
    v.check(f"{case.tag}: the window is S-orthonormal and rotation-closed "
            "(the definite-L certificate)",
            cl["lz"] <= BAND["closure"] and cl["lperp"] <= BAND["closure"],
            f"{space.ncas} active columns on {space.ncore} core pairs, "
            f"sector {space.sector}; max|C^T S C - I| {dev:.1e}; "
            f"||[P,Lz]|| {cl['lz']:.2e}, ||[P,Lperp]|| {cl['lperp']:.2e} "
            f"(band {BAND['closure']:g})")


def _electrons(states, edges) -> list[float]:
    """The ensemble's electrons per block of active columns (the mean
    1-RDM's diagonal, summed block by block)."""
    dm = np.zeros((states[0].ncas, states[0].ncas))
    for st in states:
        dm += fci.direct_spin1.make_rdm1(st.ci, st.ncas, st.nelecas)
    d = np.real(np.diag(dm)) / len(states)
    return [float(d[a:b].sum()) for a, b in zip(edges[:-1], edges[1:])]


def _term_dev(sol, term) -> float:
    """The solution's largest root deviation from the term's (<L^2>,
    <S^2>) over ALL its roots."""
    return max(max(abs(a - term[0]), abs(b - term[1]))
               for a, b in sol.meters())


def _expectations(v, case, sub, space, sols) -> None:
    """The case's declared expectations on the routes that ran."""
    if "casci" in sols:
        s = sols["casci"].states[0]
        near = min(abs(s.l2 - l * (l + 1)) for l in range(16))
        v.check(f"{case.tag}: the frozen ground carries a definite L",
                near <= BAND["definite"],
                f"<L2> {s.l2:.6f} -- nearest L(L+1) within {near:.2e} "
                f"(band {BAND['definite']:g}); <S2> {s.s2:.6f}")
        if case.frozen is not None:
            dev = max(abs(s.l2 - case.frozen[0]), abs(s.s2 - case.frozen[1]))
            v.check(f"{case.tag}: the frozen ground is the record's term",
                    dev <= BAND["l2"],
                    f"<L2> {s.l2:.6f}, <S2> {s.s2:.6f} vs "
                    f"{case.frozen} (dev {dev:.2e})")
        if case.frozen_below is not None:
            de = (s.state.e_ref - rohf_energy()) * 1e3
            v.check(f"{case.tag}: the frozen total sits near the reference",
                    de < case.frozen_below,
                    f"dE {de:+.1f} mEh (band < {case.frozen_below:g})")
        if case.s_pair:
            occ = float(np.diag(s.state.mc.make_rdm1())[
                s_column(space.active, sub.mol)])
            v.check(f"{case.tag}: the frozen ground doubly occupies the s "
                    "channel (the fix's acceptance)",
                    occ >= BAND["pair"], f"occupation {occ:.4f} (band >= {BAND['pair']:g})")
        if case.ground_config is not None:
            # <L^2> and <S^2> cannot tell one configuration of a term from
            # another: the electrons of the base's f and d blocks can
            pops = _electrons([s.state], _BASE_EDGES)
            dev = max(abs(pops[k] - case.ground_config[k]) for k in range(2))
            v.check(f"{case.tag}: the frozen ground is the ground "
                    "configuration (5f^3 6d^1)",
                    dev <= BAND["config"],
                    f"{pops[0]:.3f} / {pops[1]:.3f} electrons in the f / d "
                    f"block (declared {case.ground_config[0]:g} / "
                    f"{case.ground_config[1]:g}; band +-{BAND['config']:g})")
    if "sa" in sols:
        sa = sols["sa"]
        meters = sa.solution.meters()
        l2, s2 = [a for a, _b in meters], [b for _a, b in meters]
        if case.sa_term is not None:
            dev = _term_dev(sa.solution, case.sa_term)
            v.check(f"{case.tag}: every free-SA root carries the term",
                    dev <= BAND["l2"],
                    f"{len(l2)} roots, max dev {dev:.2e} (band "
                    f"{BAND['l2']:g}); <L2> in "
                    f"[{min(l2):.4f}, {max(l2):.4f}]")
        if case.radial:
            # the free state average on a minimal-family window: what it
            # made of the composition (it does not keep one -- the C1
            # optimizer rotates empty non-f orbitals into f)
            channels = ladder.channels_of(case.tag)
            want = [CHANNELS[c].l for c in channels
                    for _ in range(CHANNELS[c].l * 2 + 1)]
            st = sa.states[0].state
            kept = ""                       # the base: no added columns
            if want:
                _, w = shell_weights(st.mol, st.active[:, _BASE_WIDTH:],
                                     st.mf.get_ovlp())
                kept = (f"{sum(int(np.argmax(w[k])) == l for k, l in enumerate(want)) / len(want):.2f}"
                        " of the added columns still dominated by their "
                        "channel's l; ")
            v.note(f"{case.tag}: what the free state average made of the "
                   "composition",
                   f"<L2> in [{min(l2):.4f}, {max(l2):.4f}]; "
                   f"{kept}the constrained route holds the composition")
        k = int(np.argmax(sa.s2))
        state = sa.states[k].state
        ss, mult = mcaddons.spin_square(state.mc, ci=state.ci)
        v.check(f"{case.tag}: PySCF's spin_square and the meters agree on "
                "the free SA's max-<S2> root (two independent spin "
                "constructions)", abs(ss - sa.s2[k]) <= BAND["spin"],
                f"root {k}: spin_square (S^2, 2S+1) = ({ss:.8f}, "
                f"{mult:.4f}); the meters read <S2> {sa.s2[k]:.8f}")
        if case.anchor is not None:
            ref_e2 = anchor(case.anchor)["roots"][0]["e_corr"]
            d = (sa.states[0].e2 - ref_e2) * 1e3
            v.check(f"{case.tag}: the free SA's root-0 E(2) matches its "
                    "committed anchor", abs(d) <= BAND["e2_meh"],
                    f"fresh {sa.states[0].e2 * 1e3:+.4f} vs anchor "
                    f"{ref_e2 * 1e3:+.4f} mEh (delta {d:+.4f}, band "
                    f"+-{BAND['e2_meh']})")
    if "pinned-casci" in sols and case.pinned_frozen is not None:
        l2 = sols["pinned-casci"].l2
        dev = max(abs(x - case.pinned_frozen) for x in l2)
        v.check(f"{case.tag}: the pinned frozen CI reads the certified "
                "term (the gate)", dev <= BAND["l2"],
                f"{len(l2)} components, max dev {dev:.2e} (band "
                f"{BAND['l2']:g})")
    if "pinned" in sols:
        p = sols["pinned"]
        if case.pinned_sa is not None:
            v.check(f"{case.tag}: the pinned SA reads the term's ensemble "
                    "(the gate)",
                    max(abs(x - case.pinned_sa) for x in p.l2) <= BAND["l2"],
                    f"{len(p.states)} components, max dev "
                    f"{max(abs(x - case.pinned_sa) for x in p.l2):.2e}")
        if case.pinned_class is not None:
            near9 = sum(abs(x - ladder.L2_L9) <= BAND["cls"] for x in p.l2)
            near8 = sum(abs(x - ladder.L2_TERM) <= BAND["cls"] for x in p.l2)
            v.check(f"{case.tag}: the pinned components read the certified "
                    "class (the L=9 block ground, the 5L present once)",
                    (near9, near8) == case.pinned_class,
                    f"{len(p.l2)} components: {near9} at "
                    f"{ladder.L2_L9:g} +- {BAND['cls']:g}, {near8} at "
                    f"{ladder.L2_TERM:g} +- {BAND['cls']:g}, other "
                    f"{len(p.l2) - near9 - near8}")
    if "sa" in sols and "pinned" in sols:
        f, p = sols["sa"].solution, sols["pinned"].solution
        d = p.ensemble - f.states[0].state.e_ref
        detail = (f"free root 0 {f.states[0].state.e_ref:.8f} / free "
                  f"ensemble {f.ensemble:.8f} vs pinned ensemble "
                  f"{p.ensemble:.8f} (pinned - free root 0: {d:+.2e} Eh); "
                  f"E(2) means free {np.mean(sols['sa'].e2) * 1e3:+.4f} / "
                  f"pinned {np.mean(sols['pinned'].e2) * 1e3:+.4f} mEh")
        if case.pinned_sa is not None:
            v.check(f"{case.tag}: the pinned and the free SA land the same "
                    "energy (the gate)", abs(d) <= BAND["routes"],
                    detail + f" (band {BAND['routes']:g})")
        else:
            v.note(f"{case.tag}: the routes' energies", detail)
    if "radial" in sols:
        r = sols["radial"]
        c = r.solution.certificate
        if case.radial_term is not None:
            dev = _term_dev(r.solution, case.radial_term)
            v.check(f"{case.tag}: every constrained-SA root carries the term",
                    dev <= BAND["l2"],
                    f"{len(r.solution.states)} roots, max dev {dev:.2e} (band "
                    f"{BAND['l2']:g}); the route's trap check: the tracked "
                    f"roots are the converged window's lowest level")
        if case.ground_config is not None:
            pops = _electrons(list(r.solution.states), _BASE_EDGES)
            dev = max(abs(pops[k] - case.ground_config[k]) for k in range(2))
            v.check(f"{case.tag}: the constrained state average's ensemble "
                    "is the ground configuration (5f^3 6d^1)",
                    dev <= BAND["config"],
                    f"{pops[0]:.3f} / {pops[1]:.3f} electrons in the f / d "
                    f"block (declared {case.ground_config[0]:g} / "
                    f"{case.ground_config[1]:g}; band +-{BAND['config']:g})")
        v.note(f"{case.tag}: the constrained state average's converged "
               "window",
               f"the next level {c['level_gap']:.2f} mEh above the term; "
               f"largest closure meter {c['closure']:.1e}; core leak "
               f"{c['core_leak']:.1e}")


_BASE_EDGES = [0, ladder.F_BLOCK, ladder.F_BLOCK + ladder.D_BLOCK]
_BASE_WIDTH = ladder.F_BLOCK + ladder.D_BLOCK


def _ledger(tag) -> str:
    """The static / E(2) / total deltas of every row written for a tag
    (against the ROHF row, in mEh / eV)."""
    e0 = rohf_energy()
    parts = []
    for key, r in load_rows().items():
        if r.cas == tag and not r.reference and not r.withdrawn:
            parts.append(
                f"{r.method}{' pinned' if r.route is Route.PINNED else ''}"
                f"{' radial' if r.route is Route.RADIAL else ''}: "
                f"{_delta(r.e_method - e0)} | {_delta(r.e2)} | "
                f"{_delta(r.e_total - e0)}")
    return "\n    ".join(parts) if parts else "no rows"


def _increments(base) -> None:
    """Every SA row's static / E(2) / total increments against the base's
    row of its own route (a constrained row against the constrained
    base, a free one against the free base), in mEh / eV."""
    rows = load_rows()
    for route in (Route.FREE, Route.PINNED, Route.RADIAL):
        b = rows.get(f"sa-casscf|{route.value}|{base}")
        if b is None or b.withdrawn:
            continue
        print(f"  the {route.value} SA rows against {base} "
              "(static / E(2) / total increments):", flush=True)
        for key, r in rows.items():
            if (r.method == "sa-casscf" and r.route is route
                    and not r.withdrawn and r.cas != base):
                print(f"    {r.cas:16s} {_delta(r.e_method - b.e_method)} | "
                      f"{_delta(r.e2 - b.e2)} | "
                      f"{_delta(r.e_total - b.e_total)}", flush=True)


def _ladder_deltas() -> None:
    """The minimal family's ladders, route by route: each size's static /
    E(2) / total against ROHF and its increments against the line above
    (the size before), in mEh / eV."""
    rows = load_rows()
    e0 = rohf_energy()
    for route in (Route.RADIAL, Route.FREE):
        have = [rows[k] for k in (f"sa-casscf|{route.value}|{tag}"
                                  for tag in ladder.LADDER_MIN)
                if k in rows and not rows[k].withdrawn]
        if not have:
            continue
        print(f"\n  --- the minimal family on {route.value}: static / E(2) / "
              "total against ROHF, then the increments against the size "
              "before ---", flush=True)
        prev = None
        for r in have:
            s, d, t = ((r.e_method - e0) * 1e3, r.e2 * 1e3,
                       (r.e_total - e0) * 1e3)
            line = (f"    {r.cas:16s} {s:+9.3f} {d:+9.3f} {t:+9.3f} mEh   "
                    f"{s * EH_EV / 1e3:+8.4f} {d * EH_EV / 1e3:+8.4f} "
                    f"{t * EH_EV / 1e3:+8.4f} eV")
            if prev is not None:
                ds = (s - prev[0], d - prev[1], t - prev[2])
                line += (f"   | {ds[0]:+8.3f} {ds[1]:+8.3f} {ds[2]:+8.3f}"
                         f"   | {ds[0] * EH_EV / 1e3:+7.4f} "
                         f"{ds[1] * EH_EV / 1e3:+7.4f} "
                         f"{ds[2] * EH_EV / 1e3:+7.4f}")
            print(line, flush=True)
            prev = (s, d, t)


# -------------------------------------------------------------------- main

def _start() -> None:
    global _SUB
    if _SUB is None:
        _SUB = Substrate.load()


def _one_case(task) -> dict:
    """One case through its routes -- the whole procedure, its output
    captured as one text; the rows to write, in order; the verdict's
    rows.  (Runs in this process or in a worker: no store is touched
    here -- the parent alone writes.)"""
    case, want, ground_only, wrapper = task
    sub = _SUB
    v = Verdict()
    out = io.StringIO()
    sols = {}
    with contextlib.redirect_stdout(out):
        print(f"===== {case.tag} =====", flush=True)
        try:
            built = case.window(sub)
        except ValueError as e:
            print(f"  [refused] the window: {e}", flush=True)
            v.note(f"{case.tag}: the window refused by its certificate",
                   str(e)[:240])
        else:
            _window(v, case, built, sub)
            space = built.space
            if "pinned" in want and case.pinned:
                if wrapper:
                    _wrapper(v, sub, space)
                _purity(v, sub, space, case.tag)
            t = ladder.TERM
            if "casci" in want:
                sols["casci"] = _run("casci", lambda: free_casci(sub, space))
            if "casscf" in want:
                sols["casscf"] = _run("casscf", lambda: free_casscf(
                    sub, space, {**BUDGET, "max_cycle_macro": case.casscf_cycles}
                    if case.casscf_cycles else BUDGET))
            if "sa" in want:
                sols["sa"] = _run("sa", lambda: free_sa(sub, space,
                                                        ladder.NROOTS),
                                  roots=[0] if ground_only else None)
            if "pinned" in want and case.pinned:
                sols["pinned-casci"] = _run(
                    "pinned-casci (print-only)",
                    lambda: pinned_casci(sub, space, t))
                sols["pinned"] = _run("pinned-sa",
                                      lambda: pinned_sa(sub, space, t))
            if "radial" in want and case.radial:
                sols["radial"] = _run(
                    "radial-sa",
                    lambda: radial_sa(sub, space, case.radial_term,
                                      ladder.NROOTS),
                    roots=[0] if ground_only else None, as_solved=True)
            _expectations(v, case, sub, space, sols)
    rows = [(name, sols[name].row(case.tag)) for name in ROUTES
            if name in sols]
    return {"tag": case.tag, "text": out.getvalue(), "rows": rows,
            "verdict": v.rows}


def main(argv=None) -> int:
    args = _parse(sys.argv[1:] if argv is None else argv)
    v = Verdict()
    cases = ladder.select_cases(args.case_list)
    if not cases:
        print("no case selected", flush=True)
        return 0
    # the pinned routes' one-time wrapper gate runs on the first declared
    # case; with workers the cases go largest first (the long poles
    # start earliest) -- the store orders the rows, not the run
    first = cases[0]
    order = list(cases)
    if args.jobs > 1:
        order = sorted(cases, key=lambda c: -_size_of(c.tag))
    tasks = [(case, [r for r in ROUTES if r in set(args.route_list)],
              args.ground_only, case is first) for case in order]

    def finish(res) -> None:
        print(res["text"], end="", flush=True)
        for _name, row in res["rows"]:
            drop_row(row)
        print("  ledger (static / E(2) / total against ROHF):\n    "
              f"{_ledger(res['tag'])}", flush=True)
        v.adopt(res["verdict"])

    t0 = time.time()
    print(f"===== the ladder: {len(cases)} case(s); routes "
          f"{','.join(args.route_list)}; {args.jobs} job(s) x "
          f"{os.environ.get('OMP_NUM_THREADS', '?')} threads =====", flush=True)
    if args.jobs == 1:
        _start()
        for task in tasks:
            finish(_one_case(task))
    else:
        signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))
        pool = multiprocessing.get_context("spawn").Pool(
            args.jobs, initializer=_start_worker)
        try:
            for res in pool.imap_unordered(_one_case, tasks):
                finish(res)
            pool.close()
        except BaseException:                  # a killed leg takes its workers
            pool.terminate()
            raise
        finally:
            pool.join()
    print(f"\n===== the ledger ({time.time() - t0:.0f} s) =====", flush=True)
    _ladder_deltas()
    if args.base:
        print()
        _increments(args.base)
    ok = v.report()
    return 0 if ok else 1


def _size_of(tag) -> int:
    """A tag's active-space size (the long poles first)."""
    m = re.match(r"\((\d+),(\d+)\)", tag)
    return int(m.group(2)) if m else 0


def _start_worker() -> None:
    _start()


if __name__ == "__main__":
    sys.exit(main())
