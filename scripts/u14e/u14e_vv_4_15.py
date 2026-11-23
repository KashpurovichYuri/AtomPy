"""The CAS(4,15) V&V leg -- the 14p-promotion leg of the flagship
CAS-size axis.  THE PINNED ROUTE AT (4,15): ATTEMPTED AND MEASURED
(2026-11-21/22, under the corrected symmetrizer -- the first
attempt's garbage pricing was untrusted, its symmetrizer measured
impurity on unnormalized columns and destroyed one).  The leg
discloses, per run: each distilled column's SURVIVAL norm against
the frozen core span (the 14p column reads ~0.28 -- 92% core-lean:
it is one of the occupied 6p's own radials), and the irrep tails
the projection's conditioning leaves (re-projected, the disclosed
cleanup; observed ~3.5e-06).  The pinned casci then RUNS -- and its
block grounds come out NOT the 5L: their L classes read 0..42
(5S/5I classes; the same route at (4,12) reads 17 components at
72.000000).  The promoted direction destabilizes the term identity
in the pinned route -- the same class as the certified (6,13)
block-ground shift.  THE PINNED SA (the user's call, attempted on
the certified intersection window): MEASURED OUT -- the 14p columns
vanish into the pinned route's frozen core (post-subtraction
S-norm ~7e-09, the typed VanishedColumnError, re-measured every
run).  THE MECHANISM, MEASURED (the 2026-11-22 second validation,
the block-content ranks): the declared core span of @6s6p7s
realizes 6p as THREE bands (13p/14p/15p -- ``shell_bands``); the
pinned route's frozen core is the C1 reference's m-mixed occupied
6p, S-projected onto that span -- it contributes THREE independent
directions per p block, SPANNING THE 14p SLOTS ENTIRELY (the
instrumented ranks: B2u 3, B3u 3).  The promoted columns are
core-span content BY CONSTRUCTION for the pinned route: no
core-orthogonal remainder exists for them, and the pinned (4,15)
is the pinned (4,12) at this partition.  THE FREE ROUTE'S CORE
(the pur substrate) does NOT contain the slots (the survival
norms 0.28) and hosts the promotion -- its measured cost +3.09
mEh (the (c) check) is the experiment's answer.  VERDICT: the
project-then-validate promotion into lane_core is DECLINED (src
unchanged); no (6,21)-class partition is needed for now.

THE SETUP (all visible in ``main``): the trusted ROHF rung loaded;
the core = the rotation-closed span (marked); the window = the AO
columns of the 10f (7) + 13d (5) + 14p (3) shells, written out and
certified elementwise against the production builder; free CASCI /
free CASSCF / free SA on the visible window; pinned CI on the
BAND-NAMED lane window (``10f13d14p@6s6p7s`` -- the chemical-named
lane key is the refused form, the shell ``6p`` would sit in both
the window and its core; the band realizations may share rows per
``lane_core``'s own contract); PT2 through the engine's own NEVPT;
atomterm's meters.

THE ASSERTIONS:
  (a15) fresh SA E(2) root 0 ~ the committed payload row, band
        +-2 mEh;
  (c)   the (4,15) promotion WORSENS the corrected total: fresh vs
        fresh when the 4_12 leg's row is in the family store, else
        fresh (4,15) against the stored pair -- disclosed which
        reading ran.
The free CASSCF row carries NO assertion (the draw class); the
casci meters are disclosed (the certified INPUT-window ground is
the 5I -- the SA promotes it to the 5L).

Run from the repo root:
    OMP_NUM_THREADS=16 python scripts/u14e/u14e_vv_4_15.py \\
        [--ground-only]  |  tee logs/u14e_vv_4_15.log
"""

import json
import sys
import time

import numpy as np
from pyscf import mcscf

from _vv_common import (BAND_E2_MEH, ROWS, VanishedColumnError, Verdict,
                        anchor, core_part_of, drop_row, fmt_range,
                        free_kernel, lane_route, load_reference,
                        make_blocks, meters, meters_all, nevpt_rows,
                        pinned_ci_run, pinned_sa_run, price, purify_window,
                        std_row, symmetrize_columns)

ADDR_415 = "u14e_10f13d14p_4-0_sa-casscf_ao_-_pur_avg"
ADDR_412 = "u14e_10f13d_4-0_sa-casscf_ao_-_pur_avg"
ADDR_415_LANE = "u14e_10f13d14p@6s6p7s_4-0_sa-casscf_ao_-_pur_avg"
SECTOR = (4, 0)
NCORE = 5
SHELLS = (("10f", 7), ("13d", 5), ("14p", 3))
L, PARITY = 8, "u"
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

    # ---- free CASCI (the frozen input window; the pricing class) -----
    t0 = time.time()
    cas_f = mcscf.CASCI(mf, 15, SECTOR)
    cas_f.canonicalization = False
    cas_f.kernel(mo_full)
    r_f = price(cas_f, expect, roots=[0])[0]
    l_f = r_f["l2"]
    l_def = min(abs(l_f - lx * (lx + 1)) for lx in range(16))
    print(f"  [casci] e_tot {cas_f.e_tot:.8f}, <L2> {l_f:.6f} "
          f"({time.time() - t0:.1f} s)", flush=True)
    v.note("free CASCI (the frozen input window)",
           f"<L2> {l_f:.6f} -- closest L(L+1) within {l_def:.2e}; "
           f"the certified INPUT-window ground class (the 5I, L = 6, "
           f"the same reading as the (4,12) leg) -- the SA promotes "
           f"it to the 5L; E(2) {r_f['e2'] * 1e3:+.4f} mEh")
    drop_row(std_row("casci", "free", "(4,15)", float(cas_f.e_tot),
                     r_f["e2"] * 1e3,
                     float(cas_f.e_tot) + r_f["e2"],
                     f"{l_f:.4f}", f"{r_f['s2']:.4f}",
                     "5I input ground"))

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
                     r_c["e2"] * 1e3,
                     float(mc_c.e_tot) + r_c["e2"],
                     f"{r_c['l2']:.4f}", f"{r_c['s2']:.4f}",
                     "escape draw"))

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
                     f"{max(r['s2'] for r in rows):.4f}",
                     "5L, closure-costed"))

    # ---- (c) the CAS-size ordering (the flagship read) ----------------
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

    # ---- the pinned substrate (built once: disclosures + two legs) ----
    lmol, lmf, c_core_l, nf, mo_l, ncas_l, act_l = lane_route(
        campaign_spec("U14e"), ADDR_415_LANE, SECTOR)
    S_l = lmf.get_ovlp()
    # the survival norms: the PRE-projection distilled columns vs the
    # frozen core span (mo's block is already projected + renormalized
    # -- measuring it would be vacuous)
    overl = np.linalg.norm(c_core_l.T @ S_l @ act_l, axis=0)
    surv = np.sqrt(np.maximum(1.0 - overl ** 2, 0.0))
    lbl = [lab for lab, n in SHELLS for _ in range(n)]
    lean = [f"{lb} {x:.4f}" for lb, x in zip(lbl, surv) if x < 0.99]
    print(f"  [pinned window] {ncas_l} columns; survival norms vs the "
          f"frozen core span: "
          f"{'; '.join(lean) if lean else 'all at 1.0'} "
          f"(the 14p remainder is the promoted direction)", flush=True)
    # the frozen-core projection leaves tiny irrep tails (label_orb_symm
    # refuses past ~1e-6); re-project onto the dominant D2h irrep --
    # the disclosed CLEANUP (the core projection is physics, not a knob)
    c_read_l = mo_l[:, nf:nf + ncas_l]
    c_read_l, imp = symmetrize_columns(lmol, S_l, c_read_l)
    lean_in = float(np.linalg.norm(c_core_l.T @ S_l @ c_read_l,
                                   axis=0).max())
    print(f"  [pinned window] irrep tails re-projected (max "
          f"{max(imp):.1e}); the solvers' window's residual core lean "
          f"{lean_in:.1e} (the 14p slots sit inside the pinned core's "
          f"span -- the block projection cannot remove them; "
          f"disclosed)", flush=True)

    # ---- pinned casci (per-block solvers, no optimizer) ---------------
    t0 = time.time()
    roots_ci, _ = pinned_ci_run(lmol, lmf, mo_l, nf, ncas_l, SECTOR,
                                L, PARITY, c_core_l, c_read=c_read_l)
    blocks, expect_p, _ = make_blocks(lmol, lmf, mo_l, nf, ncas_l,
                                      SECTOR, L, PARITY, c_core_l,
                                      roots_ci, c_read=c_read_l)
    l2s, s2s = meters_all(blocks, expect_p)
    e2s = [row.e_corr for _, mc in blocks
           for row in nevpt_rows(mc, pt2_roots)]
    term_means = []
    for _, mc in blocks:
        ec = np.atleast_1d(mc.e_cas).ravel()
        cb = float(mc.e_tot) - float(ec[0])
        term_means.append(float(np.mean(ec + cb)))
    e_m_ci = float(np.mean(term_means))
    print(f"  [pinned casci] term e_method {e_m_ci:.8f}, E(2) mean "
          f"{np.mean(e2s) * 1e3:+.4f} mEh, ptp "
          f"{np.ptp(e2s) * 1e3:.2e} mEh, <L2> "
          f"[{min(l2s):.6f}, {max(l2s):.6f}] "
          f"({time.time() - t0:.1f} s)", flush=True)
    v.note("pinned casci (the band-named lane window, corrected "
           "symmetrizer)",
           f"<L2> in [{min(l2s):.4f}, {max(l2s):.4f}], <S2> max "
           f"{max(s2s):.4f}; E(2) mean {np.mean(e2s) * 1e3:+.4f} mEh "
           f"(the block grounds' L classes ARE the finding: the "
           f"clean-5L reading is 17 components at 72)")
    drop_row(std_row("casci", "pinned", "(4,15)", e_m_ci,
                     float(np.mean(e2s) * 1e3),
                     float(e_m_ci + np.mean(e2s)),
                     fmt_range(l2s), f"{max(s2s):.4f}",
                     "grounds not the 5L (5S/5I)"))

    # ---- pinned SA: the intersection window (the user's (a)) ----------
    # the pinned casci's block grounds are NOT the 5L (L classes
    # 0..42 -- 5S/5I classes; the same route at (4,12) reads 17 x
    # 72.000000), so the pinned SA's ensemble would be that ground
    # ensemble -- the mix kernel has no root-selection knob.  The
    # user's call: attempt it anyway on the certified intersection
    # window (purify_window) and disclose what it prices.
    # MEASURED (2026-11-22, the second validation): the promoted
    # columns are core-span content BY CONSTRUCTION for the pinned
    # route -- the declared core span of @6s6p7s realizes 6p as
    # THREE bands (13p/14p/15p), and the C1 reference's m-mixed
    # occupied 6p, S-projected onto that span, contributes THREE
    # independent directions per p block (the instrumented ranks:
    # B2u 3, B3u 3), SPANNING THE 14p SLOTS ENTIRELY: the columns
    # vanish (post-subtraction S-norm ~7e-09).
    try:
        c_sa_l, rep = purify_window(lmol, S_l, c_core_l,
                                    mo_l[:, nf:nf + ncas_l])
    except VanishedColumnError as e:
        v.note("the pinned SA at (4,15): MEASURED OUT -- the "
               "promoted band is core-span content for the pinned "
               "route",
               f"purify_window refuses: {e}.  The mechanism: the "
               f"declared core span of @6s6p7s realizes 6p as three "
               f"bands (13p/14p/15p, shell_bands), and the C1 "
               f"reference's m-mixed occupied 6p, S-projected onto "
               f"that span, spans the 14p slots entirely (the block "
               f"ranks 3/3 measured) -- no core-orthogonal remainder "
               f"exists for the promoted columns; the pinned (4,15) "
               f"is the pinned (4,12) at this partition.  The free "
               f"route's core (pur) does NOT contain the slots "
               f"(survival 0.28) and hosts the promotion: +3.09 mEh "
               f"(the (c) check).  VERDICT: the project-then-"
               f"validate promotion into lane_core is DECLINED (src "
               f"unchanged); no (6,21)-class partition is needed for "
               f"now")
        ok = v.report()
        return 0 if ok else 1
    print(f"  [pinned window] the core/irrep intersection converged "
          f"in {rep['passes']} pass(es) (core residue "
          f"{rep['core_res']:.1e}, irrep impurity "
          f"{rep['irrep_imp']:.1e}, orthonormality "
          f"{rep['ortho_dev']:.1e})", flush=True)
    t0 = time.time()
    out_sa, _ = pinned_sa_run(lmol, lmf, mo_l, nf, ncas_l, SECTOR,
                              L, PARITY, c_core_l, c_read=c_sa_l)
    blocks_s, expect_s, _ = make_blocks(lmol, lmf,
                                        np.asarray(out_sa["mo"]),
                                        nf, ncas_l, SECTOR, L, PARITY,
                                        c_core_l, out_sa["roots"])
    l2s_s, s2s_s = meters_all(blocks_s, expect_s)
    e2s_s = [row.e_corr for _, mc in blocks_s
             for row in nevpt_rows(mc, pt2_roots)]
    print(f"  [pinned sa] ensemble {out_sa['e_tot']:.8f}, term E(2) "
          f"mean {np.mean(e2s_s) * 1e3:+.4f} mEh, <L2> "
          f"[{min(l2s_s):.6f}, {max(l2s_s):.6f}], <S2> max "
          f"{max(s2s_s):.4f} ({time.time() - t0:.1f} s)", flush=True)
    v.check("(a15b) the pinned SA converges on the intersection "
            "window (the certified core/irrep fixed point; the "
            "optimized block stays irrep-pure -- the readout's "
            "labeling passed)",
            bool(out_sa["converged"]),
            f"converged {out_sa['converged']}, cycles "
            f"{out_sa['cycles']}, conv_dE {out_sa['conv_dE']:.1e}")
    delta = float(out_sa["e_tot"]) - float(mc_sa.e_tot)
    v.note("(d15) the pinned SA ensemble vs the free SA ensemble",
           f"pinned {out_sa['e_tot']:.8f} vs free {mc_sa.e_tot:.8f} "
           f"(delta {delta:+.2e} Eh -- the ensemble-surface reading; "
           f"at (4,12) the same reading held at 1.5e-09)")
    v.note("the pinned (4,15) SA prices the BLOCK-GROUND ensemble",
           f"the optimized components read <L2> in "
           f"[{min(l2s_s):.4f}, {max(l2s_s):.4f}] -- the ground "
           f"classes, NOT the 5L (the SA mix kernel has no "
           f"root-selection knob: its ensemble IS the block grounds); "
           f"the row is stored with this disclosure")
    drop_row(std_row("sa-casscf", "pinned", "(4,15)",
                     float(out_sa["e_tot"]),
                     float(np.mean(e2s_s) * 1e3),
                     float(out_sa["e_tot"] + np.mean(e2s_s)),
                     fmt_range(l2s_s), f"{max(s2s_s):.4f}",
                     "ground ensemble, not the 5L"))
    ok = v.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
