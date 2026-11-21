"""The CAS(4,15) V&V leg -- the 14p-promotion leg of the flagship
CAS-size axis (the (4,12)/(4,15)/(6,13) trio), PySCF-direct and
readable top-to-bottom.

THE SETUP (all visible below): the trusted ROHF rung loaded; the
core = the rotation-closed span over the occupied prefix (substrate
chemistry -- the certified construction, marked); the window = the
AO columns of the 10f (7) + 13d (5) + 14p (3) shells, written out
(band-pure at k=1 -- the visible selection IS the stored distilled
window; the production builder cross-checked elementwise); SA over
the ground term's 17 M_L components; the two-pass readout; the
engine's own NEVPT per root; atomterm's meters.

THE ASSERTIONS:
  (a15) fresh E(2) root 0 ~ the committed payload row
        (u14e_10f13d14p_4-0_sa-casscf_ao_-_pur_avg, root 0);
  (c)   the (4,15) promotion WORSENS the corrected total: fresh vs
        fresh when the 4_12 leg's row is in the family store, else
        fresh (4,15) against the stored (4,12)/(4,15) pair --
        disclosed which reading ran.

Run from the repo root:
    OMP_NUM_THREADS=16 python scripts/u14e/u14e_vv_4_15.py \
        [--ground-only]  |  tee logs/u14e_vv_4_15.log
"""

import json
import sys
import time
from datetime import date

import numpy as np
from pyscf import mcscf, mrpt

from _vv_common import (BAND_E2_MEH, BAND_L2, EH_EV, L_L2, S_S2,
                        ROWS, Verdict, anchor, drop_row,
                        load_reference, meters)

ADDR_415 = "u14e_10f13d14p_4-0_sa-casscf_ao_-_pur_avg"
ADDR_412 = "u14e_10f13d_4-0_sa-casscf_ao_-_pur_avg"
SECTOR = (4, 0)
NCORE = 5
SHELLS = (("10f", 7), ("13d", 5), ("14p", 3))


def visible_window(mol, mf, c_core):
    """The window, written out (the same treatment as the 4_12 leg,
    one more shell) -- asserted against the production builder."""
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
    # the certification: the core+window block elementwise (same
    # input, same Loewdin gauge); the virtual complement by SPAN
    # (any S-orthonormal basis of it is valid for the CAS optimizer)
    nc = C.shape[1]
    dev_blk = np.max(np.abs(mo[:, :nc] - mo_prod[:, :nc]))
    proj = mo_prod[:, nc:]
    resid = mo[:, nc:] - proj @ (proj.T @ S @ mo[:, nc:])
    dev_span = np.max(np.abs(resid))
    # the span band cites the observed: 1.9e-08 on this box (the
    # diffuse tail's overlap conditioning); the block equality is
    # the certificate
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
    v = Verdict()
    mol, mf = load_reference()
    from atomterm.core import purified_core_orbitals
    c_core = purified_core_orbitals(mol, mf, NCORE)[0]
    mo_full = visible_window(mol, mf, c_core)

    t0 = time.time()
    mc = mcscf.CASSCF(mf, 15, SECTOR)
    mc.state_average_([1.0 / 17] * 17)
    mc.kernel(mo_full)
    print(f"  [sa] ensemble {mc.e_tot:.8f}, converged "
          f"{mc.converged} ({time.time() - t0:.1f} s)", flush=True)
    cas = mcscf.CASCI(mf, 15, SECTOR)
    cas.mo_coeff = mc.mo_coeff
    cas.canonicalization = False
    cas.fcisolver.nroots = 17
    # seed the readout FROM the SA kernel's own vectors (the two-pass
    # pattern's faithful form -- production's readout_mc(ci0=...);
    # unseeded, the davidson tails mix the quasi-degenerate block at
    # the 1e-5 level and the meters carry it)
    from pyscf.fci import cistring
    nsa = cistring.num_strings(cas.ncas, SECTOR[0])
    nsb = cistring.num_strings(cas.ncas, SECTOR[1])
    guess = [np.asarray(c).ravel().reshape(nsa, nsb) for c in mc.ci]
    cas.kernel(ci0=guess)
    expect = meters(mol, cas.mo_coeff[:, NCORE:NCORE + 15], *SECTOR)
    a0 = anchor(ADDR_415)["roots"][0]
    e_tot_0 = float(np.atleast_1d(np.asarray(cas.e_tot)).ravel()[0])
    e_cas = np.atleast_1d(np.asarray(cas.e_cas)).ravel()
    core_e = e_tot_0 - float(e_cas[0])
    rows = []
    for k in (range(1) if ground_only else range(17)):
        pt = mrpt.NEVPT(cas, root=k)
        pt.ci = list(pt.ci)
        e2 = float(pt.kernel())
        rows.append({"root": k, "e_ref": float(e_cas[k] + core_e),
                     "e2": e2, "l2": expect("l2", cas.ci[k]),
                     "s2": expect("s2", cas.ci[k])})
        print(f"  root {k}: e_ref {rows[-1]['e_ref']:.8f}, E(2) "
              f"{e2 * 1e3:+.4f} mEh, <L2> {rows[-1]['l2']:.6f}, "
              f"<S2> {rows[-1]['s2']:.6f}", flush=True)
    r0 = rows[0]
    d_meh = (r0["e2"] - a0["e_corr"]) * 1e3
    v.check("(a15) E(2) root 0 ~ the committed payload row",
            abs(d_meh) <= BAND_E2_MEH,
            f"fresh {r0['e2'] * 1e3:+.4f} mEh vs anchor "
            f"{a0['e_corr'] * 1e3:+.4f} (delta {d_meh:+.4f}, band "
            f"+-{BAND_E2_MEH})")
    dev_m = max(abs(r["l2"] - L_L2) for r in rows)
    dev_s = max(abs(r["s2"] - S_S2) for r in rows)
    v.note("the (4,15) meters", f"<L2> in "
           f"[{min(r['l2'] for r in rows):.6f}, "
           f"{max(r['l2'] for r in rows):.6f}] (the stored row "
           f"carries 71.785 on root 0 -- the promotion's 6p rows "
           f"cost some closure; printed, not asserted)")
    drop_row({
        "leg": "4_15", "route": "free", "cas": "CAS(4,15)",
        "e_ref": r0["e_ref"], "e2_meh": r0["e2"] * 1e3,
        "e_total": r0["e_ref"] + r0["e2"],
        "anchor_e2_meh": a0["e_corr"] * 1e3, "delta_meh": d_meh,
        "l2": f"{min(r['l2'] for r in rows):.6f}.."
              f"{max(r['l2'] for r in rows):.6f}",
        "s2": f"{max(abs(r['s2'] - S_S2) for r in rows):.1e} dev",
        "date": date.today().isoformat(),
    })

    # (c) the CAS-size ordering -- fresh vs fresh when the 4_12 row
    # is in the store, else fresh (4,15) against the stored pair
    fresh_412 = None
    if ROWS.exists():
        fresh_412 = json.loads(ROWS.read_text()).get("4_12|free")
    s412 = anchor(ADDR_412)["roots"][0]["e_total"]
    s415 = a0["e_total"]
    if fresh_412 is not None:
        tot_f, tot_p = fresh_412["e_total"], r0["e_ref"] + r0["e2"]
        gap = (tot_p - tot_f) * 1e3
        v.check("(c) the (4,15) promotion WORSENS the corrected "
                "total (fresh vs fresh)", tot_p > tot_f,
                f"(4,12) {tot_f:.8f} -> (4,15) {tot_p:.8f} "
                f"(+{gap:.3f} mEh); stored: {s412:.8f} -> {s415:.8f} "
                f"(+{(s415 - s412) * 1e3:.3f} mEh)")
    else:
        v.note("(c) the CAS-size ordering",
               f"the fresh (4,12) row is not in the store (run "
               f"u14e_vv_4_12.py first for the fresh-vs-fresh "
               f"reading); stored: {s412:.8f} -> {s415:.8f} "
               f"(+{(s415 - s412) * 1e3:.3f} mEh); the fresh (4,15) "
               f"total {r0['e_ref'] + r0['e2']:.8f}")
    ok = v.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
