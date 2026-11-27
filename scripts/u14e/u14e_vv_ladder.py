"""The CAS-ladder leg -- THE WRITER of the table's method rows.

Every case of ``vv.cases`` goes through the routes the table keeps:
the free trio (PySCF-direct, on the purified core) and the pinned
state average (the production wiring, on the same window and the
same core).  Each solution is priced (strongly contracted NEVPT2),
metered (<L²>, <S²>), checked and written.  The pinned route first
runs its frozen CI (print-only: the gate's reading and the ledger's
frozen point -- the table keeps only the pinned SA).

WHAT GATES (constructions and declared records, never a physics
outcome): each window is S-orthonormal and rotation-closed; every
solution passes its birth certificate (orbitals orthonormal; CI
stationary in its own orbitals -- ``vv.routes.assert_consistent``);
the frozen ground carries a definite L; and each case's declared
expectations (``vv.cases``: the frozen term, the SA's term identity,
the committed E(2) anchors, the pin-on-pur gate at (4,12)@r1, the
certified component class at (6,13), the s channel's pair at
(6,13)@r1).

WHAT IS DISCLOSED: the per-root tables, the E(2) mean over the
term's roots (the free routes' table convention is root 0; a pinned
row is the component mean -- the means are the like-for-like
reading), the free-vs-pinned energy difference, each solution's
core leak, and the static/dynamic ledger of every row written.

A window its certificate refuses is logged ``[refused]`` and skipped
(the declared refusal class); nothing is stored.

Run from the repo root:
    OMP_NUM_THREADS=64 python scripts/u14e/u14e_vv_ladder.py \\
        [--cases "(4,12)" "(6,13)"] [--routes casci,casscf,sa,pinned] \\
        [--ground-only] [--base "(4,12)@r1"]  |  tee logs/u14e_vv_ladder.log
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
from pyscf import ao2mo
from pyscf.mcscf import addons as mcaddons

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from atomterm.irreps import orb_irreps  # noqa: E402

from vv import cases as ladder  # noqa: E402
from vv.reference import Substrate, anchor, rohf_energy  # noqa: E402
from vv.routes import (free_casci, free_casscf, free_sa,  # noqa: E402
                       pinned_casci, pinned_sa)
from vv.symmetry import symmetrize_columns  # noqa: E402
from vv.table import Route, drop_row, load_rows  # noqa: E402
from vv.verdict import Verdict  # noqa: E402
from vv.windows import s_column  # noqa: E402

ROUTES = ("casci", "casscf", "sa", "pinned")

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
    closure=1e-10,  # a window's rotation-closure meters (observed ~1e-13)
    pair=1.5,       # the s channel's occupation in the frozen ground
                    # (the pair forms at ~1.99; the as-built read 3e-03)
    spin=1e-6,      # Eh-free: PySCF's spin_square against the meters
                    # (observed <= 1e-8)
    wrapper=1e-9,   # Eh: the window's ERIs on the D2h wrapper against
                    # the C1 molecule's (observed 0.0)
    pure=1e-6,      # the active columns' irrep impurity on the wrapper
                    # (the labeling seam refuses past it)
)


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
                   help="price root 0 only on the free SA")
    p.add_argument("--base", default=None,
                   help="a case tag: the ledger prints every SA row's "
                        "increments against it")
    a = p.parse_args(argv)
    bad = [r for r in a.routes.split(",") if r not in ROUTES]
    if bad:
        p.error(f"unknown route(s) {bad}; choose from {ROUTES}")
    a.route_list = [r for r in a.routes.split(",") if r]
    a.case_list = list(a.cases)
    return a


def _print_solution(name, sol, priced, seconds) -> None:
    c = sol.certificate
    print(f"  [{name}] {seconds:.1f} s | conv {sol.conv} | certificate: "
          f"orbitals {c['ortho']:.1e}, eigenpairs {c['residual']:.1e}, "
          f"core leak {c['core_leak']:.1e}"
          + (f", readout shift {c['readout_shift']:.1e}"
             if "readout_shift" in c else "")
          + (f", orbitals as returned {c['ortho_in']:.1e}"
             if "ortho_in" in c else "")
          + (f", component mean vs ensemble {c['ensemble_dev']:.1e}"
             if "ensemble_dev" in c else ""), flush=True)
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


def _run(name, go, roots=None):
    t0 = time.time()
    sol = go()
    priced = sol.price(roots=roots)
    _print_solution(name, sol, priced, time.time() - t0)
    return priced


def _window(v, case, built, sub) -> None:
    space = built.space
    for line in built.report:
        print(f"  {line}", flush=True)
    dev = space.certify(sub.mf)
    cl = space.closure(sub.mol)
    v.check(f"{case.tag}: the window is S-orthonormal and rotation-closed "
            "(the definite-L certificate)",
            cl["lz"] <= BAND["closure"] and cl["lperp"] <= BAND["closure"],
            f"{space.ncas} active columns on {space.ncore} core pairs, "
            f"sector {space.sector}; max|C^T S C - I| {dev:.1e}; "
            f"||[P,Lz]|| {cl['lz']:.2e}, ||[P,Lperp]|| {cl['lperp']:.2e} "
            f"(band {BAND['closure']:g})")


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
                    occ >= BAND["pair"], f"occupation {occ:.4f} (band >= {BAND['pair']})")
    if "sa" in sols:
        sa = sols["sa"]
        if case.sa_term is not None:
            dev = max(max(abs(x - case.sa_term[0]) for x in sa.l2),
                      max(abs(x - case.sa_term[1]) for x in sa.s2))
            v.check(f"{case.tag}: every free-SA root carries the term",
                    dev <= BAND["l2"],
                    f"{len(sa.states)} roots, max dev {dev:.2e} (band "
                    f"{BAND['l2']:g}); <L2> in "
                    f"[{min(sa.l2):.4f}, {max(sa.l2):.4f}]")
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


def _ledger(tag) -> str:
    """The static/dynamic split of every row written for a tag (mEh
    against the ROHF row)."""
    e0 = rohf_energy()
    parts = []
    for key, r in load_rows().items():
        if r.cas == tag and not r.reference and not r.withdrawn:
            parts.append(f"{r.method}{' pinned' if r.route is Route.PINNED else ''} "
                         f"{(r.e_method - e0) * 1e3:+.2f}/{r.e2 * 1e3:+.2f}/"
                         f"{(r.e_total - e0) * 1e3:+.2f}")
    return "  ".join(parts) if parts else "no rows"


def _increments(base) -> None:
    rows = load_rows()
    b = rows.get(f"sa-casscf|{Route.FREE.value}|{base}")
    if b is None or b.withdrawn:
        print(f"  (no priced free SA row for the base {base})", flush=True)
        return
    print(f"  the SA ladder against {base} "
          f"(static/dynamic/total increments, mEh):", flush=True)
    for key, r in rows.items():
        if (r.method == "sa-casscf" and not r.withdrawn
                and (r.route is Route.FREE or r.cas == base)):
            inc = ((r.e_method - b.e_method) * 1e3, (r.e2 - b.e2) * 1e3,
                   (r.e_total - b.e_total) * 1e3)
            print(f"    {r.route.value:10s} {r.cas:10s} "
                  f"{inc[0]:+8.2f}/{inc[1]:+8.2f}/{inc[2]:+8.2f}", flush=True)


def main(argv=None) -> int:
    args = _parse(sys.argv[1:] if argv is None else argv)
    v = Verdict()
    sub = Substrate.load()
    want = set(args.route_list)
    wrapper_checked = False
    for case in ladder.select_cases(args.case_list):
        print(f"===== {case.tag} =====", flush=True)
        try:
            built = case.window(sub)
        except ValueError as e:
            print(f"  [refused] the window: {e}", flush=True)
            v.note(f"{case.tag}: the window refused by its certificate",
                   str(e)[:240])
            continue
        _window(v, case, built, sub)
        space = built.space
        if "pinned" in want:
            if not wrapper_checked:
                _wrapper(v, sub, space)
                wrapper_checked = True
            _purity(v, sub, space, case.tag)
        t = ladder.TERM
        sols = {}
        if "casci" in want:
            sols["casci"] = _run("casci", lambda: free_casci(sub, space))
        if "casscf" in want:
            sols["casscf"] = _run("casscf", lambda: free_casscf(sub, space))
        if "sa" in want:
            sols["sa"] = _run("sa", lambda: free_sa(sub, space, ladder.NROOTS),
                              roots=[0] if args.ground_only else None)
        if "pinned" in want:
            sols["pinned-casci"] = _run(
                "pinned-casci (print-only)", lambda: pinned_casci(sub, space, t))
            sols["pinned"] = _run("pinned-sa", lambda: pinned_sa(sub, space, t))
        _expectations(v, case, sub, space, sols)
        for name in ("casci", "casscf", "sa", "pinned"):
            if name in sols:
                drop_row(sols[name].row(case.tag))
        print(f"  ledger (static/dynamic/total, mEh vs ROHF): {_ledger(case.tag)}",
              flush=True)
    print("\n===== the ledger =====", flush=True)
    if args.base:
        _increments(args.base)
    ok = v.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
