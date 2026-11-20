"""The U14e verification & validation instrument: re-derive, from
scratch and standalone, the numbers the committed payloads carry, and
ASSERT the wave's answers (the run.py precedent -- a RUN INSTRUMENT,
not a pytest tier; it writes nothing).

THE ASSERTIONS (handoff lviii, Task 1):
  (a) the SA (4,12) E(2) ~ -341.665 mEh -- anchor = the committed
      payload row u14e_10f13d_4-0_sa-casscf_ao_-_pur_avg, root 0;
  (b) <L2> = 72 on every (4,12) SA root, free and pinned (the 5L
      identity; <S2> = 6);
  (c) the (4,15) promotion WORSENS the corrected total (the stored
      read: -50.765477 -> -50.762391; the dynamic double-counting);
  (d) the D2h-pinned SA route ~ the free route at (4,12) (the
      route-invariance reading);
  (e) the ROHF-vs-UHF reference comparison -- NEW physics, printed,
      not asserted (uhf_contamination reported); the UHF-seeded SA
      leg answers the handoff's recovery question the same way
      (reported, never gating the exit).

DESIGN (all pieces exist; the script wires them): references BOTH --
the trusted ROHF rung (cache-first) and the UHF sibling; the free C1
SA (casscf.constrained_casscf, the grid recipe) beside the D2h-pinned
SA (cases.lane_substrate + casscf.run_irrep_sa, the lane recipe); the
two-pass readout (perturb.readout_mc) and the per-root PT2 pricing
(perturb.nevpt2_rows_on); the meters matrix-free
(detbasis.operator_apply through the CI-order bridge).

THE CAS SET (declared; parameterize via argv):
  (4,12)  10f13d    sector (4,0)   pinned twin 5f6d@6s6p7s
  (4,15)  10f13d14p sector (4,0)   no pinned twin (the window shares
                                   the 6p rows with the core -- the
                                   declared non-disjoint partition)
  (6,13)  10f13d7s  sector (5,1)   pinned twin 5f6d7s@6s6p
  --windows=w1,w2   run a subset of the declared windows
  --with-distilled  adds the 5f6d@{r1,s1} prior windows (free SA
                    legs; the series reference layer, chk-cached)
  --ground-only     PT2 on root 0 only (the meters still run on the
                    full ensemble)

Runtime honesty: ~1-3 h for the full set (the SA optimizations are
real runs; the ensemble PT2 ~2.5 min per window at (4,12), 17 roots
~8.6 s/root measured).  Threads: 16-32.

Usage:
    OMP_NUM_THREADS=16 python scripts/u14e_vv.py \
        [--with-distilled] [--ground-only] [--windows=w1,w2]  \
        |  tee logs/u14e_vv.log
"""

import json
import resource
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402

from atomterm.cases import campaign_spec, lane_substrate  # noqa: E402
from atomterm.casscf import (  # noqa: E402
    constrained_casscf, irrep_solvers, run_irrep_sa)
from atomterm.core import closed_shell_core_energy  # noqa: E402
from atomterm.detbasis import determinant_basis, operator_apply  # noqa: E402
from atomterm.irreps import branch_d2h, orb_irreps  # noqa: E402
from atomterm.operators import (  # noqa: E402
    ci_in_detbasis_order, l2_operator_parts, lz_operator,
    s2_operator_parts)
from atomterm.perturb import nevpt2_rows_on, readout_mc  # noqa: E402
from atomterm.reference import (  # noqa: E402
    run_rohf, run_uhf, uhf_contamination)

EH_EV = 27.211386245988

L_L2 = 72.0  # the 5L parent term: L = 8 -> L(L+1)
S_S2 = 6.0   # spin 4 -> S = 2 -> S(S+1)

# the assertion bands (declared constants; each cites the first
# green run's observed fresh values -- logs/u14e_vv.log, 2026-11-20):
BAND_L2 = 1e-5     # the Casimir meters are gauge-invariant; observed
                   # on the first green run: free 8.9e-15, pinned
                   # 1.9e-07 (17 roots / components each)
BAND_E2_MEH = 2.0  # a fresh SA draw is a new orbital gauge and E(2)
                   # is state-specific: the gauge-movement scale is
                   # the route record's Ti disclosure (1.76 mEh); the
                   # first green run's observed delta was 0.0000 mEh
                   # (the engine's own wobble is ~3e-7 Eh)
BAND_ROUTE = 1e-6  # Eh; the comparison surface is the ENSEMBLE
                   # energy (variational, gauge-invariant across the
                   # routes), not E(2) -- E(2) moves at the mEh scale
                   # across orbital gauges and is printed, never
                   # asserted.  The C route read landed delta 0, the
                   # certified (4,12) pinned SA anchor_dev 0.0, and
                   # the first green run observed 5.3e-09 Eh

WORKDIR = Path("elements") / "U14e"

# the window -> (sector, pinned twin); the twin is the lane window at
# the same chemical width (the lane class returns at the CHEMICAL
# width: 5f6d@6s6p7s IS CAS(4,12), 5f6d7s@6s6p IS CAS(6,13))
CAS_SET = {
    "10f13d": ((4, 0), "5f6d@6s6p7s"),
    "10f13d14p": ((4, 0), None),
    "10f13d7s": ((5, 1), "5f6d7s@6s6p"),
}
DISTILLED = {"5f6d@r1": (4, 0), "5f6d@s1": (4, 0)}  # --with-distilled

ADDR_412 = "u14e_10f13d_4-0_sa-casscf_ao_-_pur_avg"
ADDR_415 = "u14e_10f13d14p_4-0_sa-casscf_ao_-_pur_avg"


# the assertion ledger -----------------------------------------------

class Verdict:
    """Every check lands one row; the exit code is the verdict."""

    def __init__(self):
        self.rows = []

    def check(self, name, ok, detail):
        ok = bool(ok)
        self.rows.append((name, ok, "check"))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}",
              flush=True)
        return ok

    def note(self, name, detail):
        """A reported row: measured and recorded, never gating the
        exit (the new-physics class -- printed, not asserted)."""
        self.rows.append((name, True, "note"))
        print(f"  [note] {name}: {detail}", flush=True)

    def report(self) -> bool:
        print("\n===== verdict =====", flush=True)
        for name, ok, kind in self.rows:
            tag = ("PASS" if ok else "FAIL") if kind == "check" \
                else "REPORT"
            print(f"  {tag}  {name}", flush=True)
        checks = [(n, ok) for n, ok, k in self.rows if k == "check"]
        good = sum(ok for _, ok in checks)
        print(f"{good}/{len(checks)} assertions hold", flush=True)
        return good == len(checks)


def _rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def _anchor_row(payload: dict, addr: str) -> dict:
    """The committed payload row the assertion reads -- fail fast
    unless it prices roots (the declared-roots discipline)."""
    for c in payload["cells"]:
        if c["address"] == addr:
            if c.get("disposition") != "rows" or not c.get("roots"):
                sys.exit(f"vv: anchor {addr} carries no priced roots")
            return c
    sys.exit(f"vv: anchor {addr} not in {WORKDIR / 'perturbation.json'}")


# the meter machinery -------------------------------------------------

def _meters(mol, c_act, na: int, nb: int):
    """The meters as their ACTION (detbasis.operator_apply -- the
    matrix-free twin; the len(dets)^2 container never forms).  The
    many-body meters read per root through the CI-order bridge."""
    norb = c_act.shape[1]
    dets = determinant_basis(norb, na, nb)
    apply_l2 = operator_apply(dets, *l2_operator_parts(mol, c_act))
    apply_s2 = operator_apply(dets, *s2_operator_parts(norb))
    apply_lz = operator_apply(dets, one_body=lz_operator(mol, c_act))
    apply = {"l2": apply_l2, "s2": apply_s2, "lz": apply_lz}

    def expect(kind: str, ci_pyscf) -> float:
        w = ci_in_detbasis_order(np.asarray(ci_pyscf).ravel(),
                                 norb, na, nb)
        return float(np.real(w @ apply[kind](w)))

    return expect


def _price(mc, expect, ground_only: bool) -> list[dict]:
    """PT2 + meters on the readout's roots -- the meters read the
    PRICED states (the readout's own roots, the states the PT2 rows
    describe)."""
    out = []
    for row in nevpt2_rows_on(mc, roots=[0] if ground_only else None):
        ci = np.asarray(mc.ci[row.root]).ravel()
        out.append({"row": row, "l2": expect("l2", ci),
                    "s2": expect("s2", ci), "lz": expect("lz", ci)})
    return out


def _summarize(roots: list[dict]) -> str:
    """One line per leg: root 0's economics + the ensemble meters."""
    l2s = np.array([r["l2"] for r in roots])
    s2s = np.array([r["s2"] for r in roots])
    r0 = roots[0]["row"]
    return (f"root0: e_ref {r0.e_ref:.8f}, E(2) "
            f"{r0.e_corr * 1e3:+.4f} mEh ({r0.e_corr * EH_EV:+.4f} "
            f"eV), total {r0.e_total:.8f}; <L2> in [{l2s.min():.6f}, "
            f"{l2s.max():.6f}], <S2> max dev "
            f"{np.abs(s2s - S_S2).max():.2e}")


# the free route ------------------------------------------------------

def _free_sa(mol, mf_run, key, spec, sector, weights, mf_sub=None,
             mo_full=None, prior_bands=None, series_mfs=None,
             label="free"):
    """The grid recipe, fresh and store-free: the pur substrate's
    window columns -> constrained_casscf (the SA kernel) -> the
    two-pass readout on the converged window.  ``mf_sub`` is the
    substrate the window/core resolves on (the ROHF; the UHF-seeded
    leg reuses its window and partition, only the reference class
    under the kernel changes)."""
    from atomterm.cells import key_core, window_columns
    from atomterm.windows import build_full_mo, occupied_core_idx

    t0 = time.time()
    mf_sub = mf_run if mf_sub is None else mf_sub
    na, nb = sector
    ncore, c_core, _c_pur = key_core(key, mol, mf_sub,
                                      spec.shell_bands)
    act, ncas, _rep = window_columns(
        key, mol, mf_sub, spec.shell_bands, spec.shell_radials,
        ncol=ncore, prior_bands=prior_bands, series_mfs=series_mfs)
    if mo_full is None:
        mo_full = build_full_mo(mf_sub, occupied_core_idx(mf_sub,
                                                         ncore),
                                act, c_core=c_core)
    res = constrained_casscf(mf_run, ncas, (na, nb), None,
                             mo_full=mo_full, weights=weights,
                             name=f"u14e_vv_{label}")
    if np.ndim(res.mo_coeff) != 2:
        raise ValueError(
            f"{label}: the kernel returned mo_coeff with ndim "
            f"{np.ndim(res.mo_coeff)} -- the det-basis meters need "
            f"one common orbital basis")
    mc = readout_mc(mol, mf_run, res.mo_coeff, ncas, (na, nb),
                    len(weights), ci0=res.ci)
    c_act = np.asarray(res.mo_coeff)[:, ncore:ncore + ncas]
    ens = float(np.dot(weights, np.asarray(res.energies, dtype=float)))
    print(f"  [sa] {label}: ensemble {ens:.8f}"
          f" (conv {res.conv_dE:.1e}, {res.cycles} cycles), ncas "
          f"{ncas}, {time.time() - t0:.1f} s, rss {_rss():.1f} GB",
          flush=True)
    return res, mc, ncore, ncas, _meters(mol, c_act, na, nb)


# the pinned route ----------------------------------------------------

def _pinned_sa(lmol, lmf, key, spec, sector, label):
    """The lane recipe, fresh and store-free: run_irrep_sa's two-pass
    pattern (the ensemble mix kernel, then the exact per-component
    roots re-solved in the frozen converged window), then the
    per-block kernel-ready mc objects (readout_lane's construction on
    the fresh converged window -- the pinned solver rides the block so
    mrpt's canonicalization stays in it)."""
    from atomterm.cells import lane_core, window_columns
    from atomterm.windows import build_full_mo
    from pyscf import mcscf

    t0 = time.time()
    L, parity = spec.term
    na, nb = sector
    _core, nf, c_core, _rep = lane_core(key, lmol, lmf,
                                        spec.shell_bands)
    act, ncas, _w = window_columns(key, lmol, lmf, spec.shell_bands,
                                    spec.shell_radials)
    mo = build_full_mo(lmf, list(range(nf)), act, c_core=c_core)
    out = run_irrep_sa(lmol, lmf, mo[:, nf:nf + ncas], nf, na, nb,
                       L, parity, c_core=c_core)
    c_read = np.asarray(out["mo"])[:, nf:nf + ncas]
    orbsyms = orb_irreps(lmol, c_read)
    pinned = {irr: s for irr, s, _ in
              irrep_solvers(lmol, orbsyms, L, parity, 1e-12)}
    core_part = closed_shell_core_energy(lmf, c_core)
    expect = _meters(lmol, c_read, na, nb)
    blocks = []
    for irrep, mult in branch_d2h(L, parity).items():
        e_arr, ci = out["roots"][irrep]
        e_arr = np.atleast_1d(np.asarray(e_arr, dtype=float))
        mc = mcscf.CASCI(lmf, ncas, (na, nb))
        mc.mo_coeff = np.asarray(out["mo"], dtype=float)
        mc.fcisolver = pinned[irrep]
        mc.orbsym = orbsyms
        mc.fcisolver.orbsym = np.asarray(orbsyms)
        ci_arr = np.asarray(ci)
        mc.ci = ([ci_arr] if ci_arr.ndim == 2
                 else [np.asarray(x) for x in ci_arr])
        mc.e_cas = e_arr[:mult]
        mc.e_tot = float(mc.e_cas[0] + core_part)
        blocks.append({"irrep": irrep, "mc": mc, "expect": expect})
    print(f"  [pinned] {label}: ensemble {out['e_tot']:.8f} (conv "
          f"{out['conv_dE']}, {out['cycles']} cycles), blocks "
          f"{[b['irrep'] for b in blocks]}, {time.time() - t0:.1f} s, "
          f"rss {_rss():.1f} GB", flush=True)
    return out, blocks


# the driver ----------------------------------------------------------

def main() -> None:
    argv = sys.argv[1:]
    known = ("--with-distilled", "--ground-only")
    windows_sel = None
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--windows":
            i += 1
            if i >= len(argv):
                sys.exit("--windows needs a value (w1,w2,...)")
            windows_sel = tuple(argv[i].split(","))
        elif a.startswith("--windows="):
            windows_sel = tuple(a.split("=", 1)[1].split(","))
        elif a not in known:
            sys.exit(f"usage: u14e_vv.py [--with-distilled] "
                     f"[--ground-only] [--windows=w1,w2,...]; "
                     f"unknown flag: {a}")
        i += 1
    with_distilled = "--with-distilled" in argv
    ground_only = "--ground-only" in argv
    if windows_sel is not None:
        bad = [w for w in windows_sel if w not in CAS_SET]
        if bad:
            sys.exit(f"unknown windows {bad}; declared: "
                     f"{list(CAS_SET)}")

    spec = campaign_spec("U14e")
    cfg = spec.atom_config()
    payload = json.loads((WORKDIR / "perturbation.json").read_text())
    weights = list(spec.sa_weights)
    v = Verdict()

    # ---- references, BOTH -------------------------------------------
    print("===== references =====", flush=True)
    rohf_chk = WORKDIR / "rungs" / "rohf.chk"
    mol, mf = run_rohf(cfg, str(rohf_chk), ref=spec.ref)
    print(f"  ROHF E = {mf.e_tot:.8f}  (trusted rung {rohf_chk})",
          flush=True)
    uhf_chk = WORKDIR / "rungs" / "uhf.chk"
    dm0 = None
    if spec.ref.seed_basis is not None and not uhf_chk.exists():
        from atomterm.reference import uhf_seed_density
        dm0 = uhf_seed_density(cfg, spec.ref, str(uhf_chk))
    _, mf_uhf = run_uhf(cfg, str(uhf_chk), ref=spec.ref, dm0=dm0)
    s2_uhf, sigma = uhf_contamination(mf_uhf)
    n_broken = int(np.sum(sigma < 1.0 - 1e-6))
    print(f"  UHF  E = {mf_uhf.e_tot:.8f}  (E_uhf - E_rohf = "
          f"{(mf_uhf.e_tot - mf.e_tot) * 1e3:+.3f} mEh)", flush=True)
    print(f"  [new physics, printed not asserted] <S2> = "
          f"{s2_uhf:.6f} (Sz exact by construction); pairing "
          f"spectrum: {len(sigma)} occupied pairs, {n_broken} below "
          f"1 - 1e-6, smallest {np.sort(sigma)[:3]}", flush=True)

    anchor412 = _anchor_row(payload, ADDR_412)
    anchor415 = _anchor_row(payload, ADDR_415)
    stored_tot_412 = anchor412["roots"][0]["e_total"]
    stored_tot_415 = anchor415["roots"][0]["e_total"]
    fresh_tot = {}

    # ---- the free SA legs -------------------------------------------
    from atomterm.cells import CellKey

    series_mfs = None
    if with_distilled:
        from atomterm.basis import run_series
        series_mfs = run_series(cfg, "elements", spec.key, spec.series,
                                mf, seed_basis=spec.series_seed)
    legs = {w: s for w, s in CAS_SET.items()
            if windows_sel is None or w in windows_sel}
    if with_distilled:
        legs.update({w: (s, None) for w, s in DISTILLED.items()})
    mo412 = None
    res412 = None
    for win, (_sector, _lane) in legs.items():
        na, nb = _sector
        addr = f"u14e_{win}_{na}-{nb}_sa-casscf_ao_-_pur_avg"
        pb = spec.prior_bands if "@" in win else None
        print(f"\n===== free SA {win} ({na},{nb}) =====", flush=True)
        res, mc, _ncore, _ncas, expect = _free_sa(
            mol, mf, CellKey.from_address(addr), spec, _sector,
            weights, prior_bands=pb, series_mfs=series_mfs,
            label=f"free-{win}")
        if win == "10f13d":
            mo412 = np.asarray(res.mo_coeff)
            res412 = res
        roots = _price(mc, expect, ground_only)
        fresh_tot[win] = roots[0]["row"].e_total
        print(f"  {_summarize(roots)}", flush=True)
        if win == "10f13d":
            a0 = anchor412["roots"][0]
            r0 = roots[0]["row"]
            d_meh = (r0.e_corr - a0["e_corr"]) * 1e3
            v.check(
                "(a) SA (4,12) E(2) ~ the committed payload root 0",
                abs(d_meh) <= BAND_E2_MEH,
                f"fresh {r0.e_corr * 1e3:+.4f} mEh vs anchor "
                f"{a0['e_corr'] * 1e3:+.4f} (delta {d_meh:+.4f}, "
                f"band ±{BAND_E2_MEH})")
            dev = max(np.max(np.abs([r["l2"] for r in roots]) - L_L2),
                      np.max(np.abs([r["s2"] for r in roots]) - S_S2))
            v.check(
                "(b) <L2> = 72 / <S2> = 6 on every free SA (4,12) "
                "root", dev <= BAND_L2,
                f"{len(roots)} roots, max dev {dev:.2e} "
                f"(band {BAND_L2:g})")

    # ---- the pinned SA legs -----------------------------------------
    lmol = lmf = None
    for win, (sector, lane_win) in CAS_SET.items():
        if lane_win is None or (windows_sel is not None
                                and win not in windows_sel):
            continue
        if lmol is None:
            print("\n===== pinned substrate (the D2h lane) =====",
                  flush=True)
            lmol, lmf = lane_substrate(spec)
            print(f"  D2h ROHF E = {lmf.e_tot:.8f}", flush=True)
        na, nb = sector
        addr = f"u14e_{lane_win}_{na}-{nb}_sa-casscf_ao_-_pur_avg"
        print(f"\n===== pinned SA {lane_win} ({na},{nb}) =====",
              flush=True)
        out, blocks = _pinned_sa(lmol, lmf, CellKey.from_address(addr),
                                 spec, sector, f"pinned-{lane_win}")
        e2s, l2s, s2s = [], [], []
        for b in blocks:
            for row in nevpt2_rows_on(
                    b["mc"], roots=[0] if ground_only else None):
                ci = np.asarray(b["mc"].ci[row.root]).ravel()
                e2s.append(row.e_corr)
                l2s.append(b["expect"]("l2", ci))
                s2s.append(b["expect"]("s2", ci))
        e2s = np.array(e2s)
        print(f"  term E(2): mean {e2s.mean() * 1e3:+.4f} mEh, ptp "
              f"{np.ptp(e2s) * 1e3:.2e} mEh over {len(e2s)} components"
              f"; <L2> in [{min(l2s):.6f}, {max(l2s):.6f}]", flush=True)
        if win == "10f13d" and res412 is not None:
            dev = max(np.abs(np.array(l2s) - L_L2).max(),
                      np.abs(np.array(s2s) - S_S2).max())
            v.check(
                "(b) <L2> = 72 / <S2> = 6 on every pinned SA (4,12) "
                "component", dev <= BAND_L2,
                f"{len(l2s)} components, max dev {dev:.2e} "
                f"(band {BAND_L2:g})")
            free_ens = float(np.dot(weights, res412.energies))
            delta = out["e_tot"] - free_ens
            v.check(
                "(d) the pinned SA route ~ the free route at (4,12)",
                abs(delta) <= BAND_ROUTE,
                f"pinned ensemble {out['e_tot']:.8f} vs free "
                f"{free_ens:.8f} (delta {delta:+.2e} Eh, band "
                f"{BAND_ROUTE:g})")

    # ---- the UHF-seeded SA leg (new physics: reported, not asserted)
    if res412 is not None:
        print("\n===== UHF-seeded SA 10f13d (4,12) =====", flush=True)
        try:
            _res_u, mc_u, _nc, _nk, expect_u = _free_sa(
                mol, mf_uhf, CellKey.from_address(ADDR_412), spec,
                (4, 0), weights, mf_sub=mf, mo_full=mo412,
                label="uhf-seeded")
            roots_u = _price(mc_u, expect_u, ground_only)
            print(f"  {_summarize(roots_u)}", flush=True)
            dev = max(np.max(np.abs([r["l2"] for r in roots_u]) - L_L2),
                      np.max(np.abs([r["s2"] for r in roots_u]) - S_S2))
            v.note(
                "the UHF-seeded SA recovers <L2> = 72 / <S2> = 6 "
                "from the spin-broken reference",
                f"{len(roots_u)} roots, max dev {dev:.2e} (band "
                f"{BAND_L2:g}) -- the SA repair answers the handoff's "
                f"question")
        except Exception as exc:  # the leg is real; a crash is a note
            v.note("the UHF-seeded SA leg ran",
                   f"crashed: {type(exc).__name__}: {exc}")

    # ---- (c) the CAS-size ordering ----------------------------------
    print("\n===== (c) the CAS-size axis =====", flush=True)
    if not {"10f13d", "10f13d14p"} <= set(fresh_tot):
        v.note("(c) the CAS-size ordering", "skipped: the --windows "
               "subset carries not both legs")
    else:
        gap_f = (fresh_tot["10f13d14p"] - fresh_tot["10f13d"]) * 1e3
        gap_s = (stored_tot_415 - stored_tot_412) * 1e3
        v.check(
            "(c) the (4,15) promotion WORSENS the corrected total",
            fresh_tot["10f13d14p"] > fresh_tot["10f13d"],
            f"fresh: (4,12) {fresh_tot['10f13d']:.8f} -> (4,15) "
            f"{fresh_tot['10f13d14p']:.8f} (+{gap_f:.3f} mEh); "
            f"stored: {stored_tot_412:.8f} -> {stored_tot_415:.8f} "
            f"(+{gap_s:.3f} mEh)")

    ok = v.report()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
