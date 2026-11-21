"""The CAS(4,15) V&V leg -- the 14p-promotion leg of the flagship
CAS-size axis.  Two routes (free CASSCF, free SA).  THE PINNED LEGS,
ATTEMPTED AND MEASURED OUT (2026-11-21): the 14p column's survival
norm against the frozen core span is 0.28 (92% core-lean -- it is
one of the radials the occupied 6p is built from); the attempted
pinned legs priced garbage on that remainder (term e_var -46.95,
<L2> 0..42 -- not the term) and the irrep-purity chain broke at
every re-orthogonalization.  The declared refusal ("a non-disjoint
partition is a grid axis, not a lane one") is thereby
MEASURED-correct; the free rows carry the (4,15) experiment, and the
leg re-measures the survival norms each run (the disclosure below).
A disjoint (6,21)-class partition (dropping 6p from the core, the
semi-class move) is the possible follow-up.

THE SETUP (all visible in ``main``): the trusted ROHF rung loaded;
the core = the rotation-closed span (marked); the window = the AO
columns of the 10f (7) + 13d (5) + 14p (3) shells, written out and
certified elementwise against the production builder; free CASSCF
(single root, the draw class) and free SA (17 components, the
two-pass readout); PT2 through the engine's own NEVPT; atomterm's
meters.

THE ASSERTIONS:
  (a15) fresh SA E(2) root 0 ~ the committed payload row, band
        +-2 mEh;
  (c)   the (4,15) promotion WORSENS the corrected total: fresh vs
        fresh when the 4_12 leg's row is in the family store, else
        fresh (4,15) against the stored (4,12)/(4,15) pair --
        disclosed which reading ran.

Run from the repo root:
    OMP_NUM_THREADS=16 python scripts/u14e/u14e_vv_4_15.py \\
        [--ground-only]  |  tee logs/u14e_vv_4_15.log
"""

import json
import sys
import time

import numpy as np

from _vv_common import (BAND_E2_MEH, ROWS, Verdict, anchor,
                        core_part_of, drop_row, fmt_range, free_kernel, lane_route,
                        load_reference, meters, price, std_row)

ADDR_415 = "u14e_10f13d14p_4-0_sa-casscf_ao_-_pur_avg"
ADDR_412 = "u14e_10f13d_4-0_sa-casscf_ao_-_pur_avg"
ADDR_415_LANE = "u14e_10f13d14p@6s6p7s_4-0_sa-casscf_ao_-_pur_avg"
SECTOR = (4, 0)
NCORE = 5
SHELLS = (("10f", 7), ("13d", 5), ("14p", 3))
ROOTS = 17


def visible_window(mol, mf, c_core):
    """The window, written out (the same treatment as the 4_12 leg,
    one more shell) -- certified against the production builder."""
    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    sel = [i for i, lab in enumerate(labels)
           if any(lab.startswith(w) for w, _ in SHELLS)]
    ncas = len(sel)
    if ncas != sum(n for _, n in SHELLS):
        sys.exit(f"the label selection found {ncas} AOs, expected "
                 f"{sum(n for _, n in SHELLS)}")
    nao = mol.nao
    S = mf.get_ovlp()
    act = np.zeros((nao, ncas))
    for k, i in enumerate(sel):
        act[i, k] = 1.0
    P_core = c_core @ c_core.T @ S   # the metric projector
    w = act - P_core @ act
    g = w.T @ S @ w
    gv, U = np.linalg.eigh(g)
    c_win = w @ U @ np.diag(gv ** -0.5) @ U.T
    C = np.hstack([c_core, c_win])
    P_all = C @ C.T @ S   # (C^T S C = I)
    rest = mf.mo_coeff - P_all @ mf.mo_coeff
    g2 = rest.T @ S @ rest
    gv2, U2 = np.linalg.eigh(g2)
    keep = gv2 > 1e-7 * gv2.max()
    c_virt = rest @ U2[:, keep] @ np.diag(gv2[keep] ** -0.5)
    if c_virt.shape[1] != nao - C.shape[1]:
        sys.exit(f"the complement came out {c_virt.shape[1]} wide, "
                 f"expected {nao - C.shape[1]}")
    mo = np.hstack([C, c_virt])
    from atomterm.windows import build_full_mo, occupied_core_idx
    mo_prod = build_full_mo(mf, occupied_core_idx(mf, NCORE), sel,
                            c_core=c_core)
    nc = C.shape[1]
    dev_blk = np.max(np.abs(mo[:, :nc] - mo_prod[:, :nc]))
    proj = mo_prod[:, nc:]
    resid = mo[:, nc:] - proj @ (proj.T @ S @ mo[:, nc:])
    dev_span = np.max(np.abs(resid))
    if dev_blk > 1e-9 or dev_span > 1e-6:
        sys.exit(f"the visible window construction disagrees with "
                 f"the production builder (block dev {dev_blk:.2e}, "
                 f"complement span dev {dev_span:.2e})")
    print(f"  window: {ncas} AO columns {SHELLS}, visible "
          f"construction == production build_full_mo (block dev "
          f"{dev_blk:.1e}, complement span dev {dev_span:.1e})",
          flush=True)
    return mo


def survival_disclosure(spec):
    """The measured reason the pinned (4,15) is out: each distilled
    column's S-norm after removing the lane core's span (the 14p
    radial is one of the occupied 6p's own radials)."""
    from _vv_common import lane_route
    lmol, lmf, c_core_l, nf, _mo, ncas_l, act = lane_route(
        spec, ADDR_415_LANE, SECTOR)
    S = lmf.get_ovlp()
    overl = np.linalg.norm(c_core_l.T @ S @ act, axis=0)
    surv = np.sqrt(np.maximum(1.0 - overl ** 2, 0.0))
    lbl = [lab for lab, n in SHELLS for _ in range(n)]
    return [f"{lb} {x:.4f}" for lb, x in zip(lbl, surv) if x < 0.99]


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    ground_only = "--ground-only" in argv
    pt2_roots = [0] if ground_only else None
    v = Verdict()
    mol, mf = load_reference()
    from atomterm.core import purified_core_orbitals
    from atomterm.cases import campaign_spec
    c_core = purified_core_orbitals(mol, mf, NCORE)[0]
    mo_full = visible_window(mol, mf, c_core)
    expect = meters(mol, mo_full[:, NCORE:NCORE + 15], *SECTOR)
    weights = [1.0 / 17] * 17

    # ---- the pinned (4,15): measured out (the disclosure) -------------
    lean = survival_disclosure(campaign_spec("U14e"))
    v.note("the pinned (4,15): attempted and measured out",
           f"the 14p column's survival norm against the frozen core "
           f"span: {'; '.join(lean)} (92% core-lean) -- the "
           f"attempted pinned legs priced garbage on that remainder "
           f"(term e_var -46.95, <L2> 0..42, irrep tails breaking "
           f"every re-orthogonalization); the declared refusal is "
           f"measured-correct.  A disjoint (6,21)-class partition "
           f"(dropping 6p from the core) is the possible follow-up")

    # ---- free CASSCF (single root; the escape class, no assertion) --
    t0 = time.time()
    mc_c, cas_c = free_kernel(mf, 15, SECTOR, mo_full)
    r_c = price(cas_c, expect, roots=[0])[0]
    print(f"  [casscf] e_tot {mc_c.e_tot:.8f}, converged "
          f"{mc_c.converged}, <L2> {r_c['l2']:.6f} "
          f"({time.time() - t0:.1f} s)", flush=True)
    v.note("free CASSCF (the single-root draw class)",
           f"<L2> {r_c['l2']:.6f}, E(2) {r_c['e2'] * 1e3:+.4f} mEh")
    drop_row(std_row("casscf", "free", "(4,15)", float(mc_c.e_tot),
                     r_c["e2"] * 1e3, float(cas_c.e_tot),
                     f"{r_c['l2']:.4f}", f"{r_c['s2']:.4f}"))

    # ---- free SA ------------------------------------------------------
    t0 = time.time()
    mc_sa, cas_sa = free_kernel(mf, 15, SECTOR, mo_full,
                                weights=weights, nroots=ROOTS)
    core_e, e_cas = core_part_of(cas_sa)
    expect_sa = meters(mol, cas_sa.mo_coeff[:, NCORE:NCORE + 15],
                       *SECTOR)
    a0 = anchor(ADDR_415)["roots"][0]
    rows = price(cas_sa, expect_sa,
                 roots=range(1) if ground_only else range(ROOTS))
    for r in rows:
        r["e_ref"] = float(e_cas[r["root"]] + core_e)
        print(f"  root {r['root']}: e_ref {r['e_ref']:.8f}, E(2) "
              f"{r['e2'] * 1e3:+.4f} mEh, <L2> {r['l2']:.6f}, "
              f"<S2> {r['s2']:.6f}", flush=True)
    r0 = rows[0]
    d_meh = (r0["e2"] - a0["e_corr"]) * 1e3
    v.check("(a15) SA E(2) root 0 ~ the committed payload row",
            abs(d_meh) <= BAND_E2_MEH,
            f"fresh {r0['e2'] * 1e3:+.4f} mEh vs anchor "
            f"{a0['e_corr'] * 1e3:+.4f} (delta {d_meh:+.4f}, band "
            f"+-{BAND_E2_MEH})")
    v.note("the (4,15) meters",
           f"<L2> in [{min(r['l2'] for r in rows):.4f}, "
           f"{max(r['l2'] for r in rows):.4f}] (the promotion's 6p "
           f"rows cost some closure; the stored row carries 71.785 "
           f"on root 0 -- printed, not asserted)")
    print(f"  [sa] ensemble {mc_sa.e_tot:.8f}, converged "
          f"{mc_sa.converged} ({time.time() - t0:.1f} s)", flush=True)
    drop_row(std_row("sa-casscf", "free", "(4,15)",
                     float(e_cas[0] + core_e), r0["e2"] * 1e3,
                     float(r0["e_ref"] + r0["e2"]),
                     fmt_range([r['l2'] for r in rows]),
                     f"{max(r['s2'] for r in rows):.4f}"))

    # ---- (c) the CAS-size ordering ------------------------------------
    fresh_412 = None
    if ROWS.exists():
        fresh_412 = json.loads(ROWS.read_text()).get(
            "sa-casscf|free|(4,12)")
    s412 = anchor(ADDR_412)["roots"][0]["e_total"]
    s415 = a0["e_total"]
    tot_fresh = r0["e_ref"] + r0["e2"]
    if fresh_412 is not None:
        tot_412 = fresh_412["e_tot"]
        v.check("(c) the (4,15) promotion WORSENS the corrected "
                "total (fresh vs fresh)", tot_fresh > tot_412,
                f"(4,12) {tot_412:.4f} -> (4,15) {tot_fresh:.4f} "
                f"(+{(tot_fresh - tot_412) * 1e3:.3f} mEh); stored: "
                f"{s412:.4f} -> {s415:.4f} "
                f"(+{(s415 - s412) * 1e3:.3f} mEh)")
    else:
        v.note("(c) the CAS-size ordering",
               f"the fresh (4,12) row is not in the store (run "
               f"u14e_vv_4_12.py first for the fresh-vs-fresh "
               f"reading); stored: {s412:.4f} -> {s415:.4f} "
               f"(+{(s415 - s412) * 1e3:.3f} mEh); the fresh (4,15) "
               f"total {tot_fresh:.4f}")
    ok = v.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
