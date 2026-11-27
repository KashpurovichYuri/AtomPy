"""The (6,13) active-space identity probe -- do the routes and the
windows reach the SAME active space?  The energies cannot answer on a
quasi-degenerate manifold (the class spread ~0.1 mEh degenerates every
partition); the identity instruments:

  1. the SPAN PRINCIPAL ANGLES -- the converged active subspaces'
     cross-Gram singular values under the AO overlap (the gauge-
     invariant subspace comparison: identical spans -> all cosines
     1.0);  the C1 and D2h substrates share the AO basis and overlap
     (asserted before use);
  2. the OCCUPATION FINGERPRINTS -- the ensemble-average 1-RDM
     diagonals per active column (equal weight over the 17 roots on
     both sides), compared as sorted profiles (the column gauges
     differ across routes; the profile is the gauge-free witness);
  3. the <L2> PARTITION WALKS -- the per-root <L2> per run, snapped
     to the L(L+1) grid (the term-partition identity: the free's
     L=9 partition vs the pinned's counts).

THE THREE RUNS (the two-axis comparison):
  the route axis at the fixed as-built window:
     free SA (C1)  vs  pinned SA (pinned@pur);
  the window axis at the fixed pinned@pur route:
     as-built window  vs  @r1 window (act_s_mo_of's mixed-source
     s channel -- the leg's own construction, reused verbatim).

Prints only (no store writes).  Run from the repo root:
    OMP_NUM_THREADS=64 python scripts/probes/probe_u14e_span_compare.py \\
        | tee logs/probe_u14e_span_compare.log
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "u14e"))

import numpy as np  # noqa: E402

from vv.meters import meters  # noqa: E402
from vv.reference import load_reference  # noqa: E402
from vv.routes import free_kernel, pinned_sa_run  # noqa: E402

SECTOR = (5, 1)
NCORE = 4
ROOTS = 17
BUDGET = {"max_cycle_macro": 200}
ADDR_613 = "u14e_10f13d7s_5-1_sa-casscf_ao_-_pur_avg"


def root_rdm_diag(cis, ncas, na, nb):
    """The ensemble-average 1-RDM diagonal over the given CI roots
    (equal root weights -- the free ensemble's [1/17] weighting).
    ``direct_spin1.make_rdm1`` returns the 1-RDM matrix itself."""
    from pyscf.fci import direct_spin1
    dms = []
    for c in cis:
        out = direct_spin1.make_rdm1(np.asarray(c), ncas, (na, nb))
        dm = out[0] if isinstance(out, tuple) else out
        dms.append(np.asarray(dm))
    return np.diag(np.mean(dms, axis=0)).real


def block_cis(cii):
    """The CI list from a run_irrep_casci roots entry: a single 2-D
    root stays a one-element list; a multi-root entry is the list."""
    arr = np.asarray(cii)
    if arr.ndim == 2:
        return [arr]
    return [np.asarray(x) for x in cii]


def snap_counts(vals, tol=1e-3):
    """The {L(L+1): count} census of the <L2> values (the nearest
    grid point within tol; 'other' otherwise)."""
    counts = {}
    for v in vals:
        ll = round((np.sqrt(1 + 4 * v) - 1) / 2)
        hit = ll * (ll + 1)
        key = str(int(hit)) if abs(v - hit) < tol else "other"
        counts[key] = counts.get(key, 0) + 1
    return counts


def main() -> int:
    t00 = time.time()
    mol, mf = load_reference()
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    import u14e_vv_unification as uni
    mol_d2h, mf_d2h = uni.wrapper_substrate(mol, mf, spec)

    # the shared-AO-overlap assertion (the cross-Gram's precondition)
    s1 = mol.intor_symmetric("int1e_ovlp")
    s2 = mol_d2h.intor_symmetric("int1e_ovlp")
    print(f"[overlap] |S_c1 - S_d2h|_max "
          f"{np.max(np.abs(s1 - s2)):.2e} (the cross-Gram precondition)",
          flush=True)

    from atomterm.windows import build_full_mo, occupied_core_idx
    from atomterm.cells import CellKey, key_core, window_columns
    key = CellKey.from_address(ADDR_613)
    nc, _cc, c_pur = key_core(key, mol, mf, spec.shell_bands)
    act, ncas, _ = window_columns(key, mol, mf, spec.shell_bands,
                                  spec.shell_radials, ncol=nc)
    weights = [1.0 / ROOTS] * ROOTS

    runs = {}

    # ---- run 1: the free SA (C1, the as-built window, free@pur) ------
    print("===== 1. the free SA (C1, the as-built window) =====",
          flush=True)
    t0 = time.time()
    mo_f = build_full_mo(mf, occupied_core_idx(mf, NCORE), act,
                         c_core=c_pur)
    mc1, cas1 = free_kernel(mf, ncas, SECTOR, mo_f, weights=weights,
                            nroots=ROOTS, track_iters=True)
    print(f"  E {mc1.e_tot:.8f}, converged {mc1.converged} "
          f"({time.time() - t0:.0f} s)", flush=True)
    runs["free-C1"] = {
        "span": np.asarray(cas1.mo_coeff)[:, NCORE:NCORE + ncas],
        "ci": [np.asarray(c) for c in cas1.ci], "ncas": ncas,
        "c_act": cas1.mo_coeff[:, NCORE:NCORE + ncas], "mol": mol}

    # ---- run 2: the pinned SA (pinned@pur, the as-built window) ------
    print("===== 2. the pinned SA (pinned@pur, the as-built window) "
          "=====", flush=True)
    t0 = time.time()
    mo_d2h = build_full_mo(mf_d2h, occupied_core_idx(mf_d2h, NCORE),
                           act, c_core=c_pur)
    out2, _ = pinned_sa_run(mol_d2h, mf_d2h, mo_d2h, NCORE, ncas,
                            SECTOR, 8, "u", c_pur, mc_opts=BUDGET)
    print(f"  ensemble {out2['e_tot']:.8f}, converged "
          f"{out2['converged']} (cycles {out2['cycles']}, "
          f"{time.time() - t0:.0f} s)", flush=True)
    runs["pinned-asbuilt"] = {
        "span": np.asarray(out2["mo"])[:, NCORE:NCORE + ncas],
        "roots": out2["roots"], "ncas": ncas, "mol": mol_d2h}

    # ---- run 3: the pinned SA (pinned@pur, the @r1 window) -----------
    print("===== 3. the pinned SA (pinned@pur, the @r1 window) =====",
          flush=True)
    t0 = time.time()
    act_s, ncas_s, _purity = uni.act_s_mo_of(mol, mf, mol_d2h, spec)
    mo_s = build_full_mo(mf_d2h, occupied_core_idx(mf_d2h, NCORE),
                         act_s, c_core=c_pur)
    out3, _ = pinned_sa_run(mol_d2h, mf_d2h, mo_s, NCORE, ncas_s,
                            SECTOR, 8, "u", c_pur, mc_opts=BUDGET)
    print(f"  ensemble {out3['e_tot']:.8f}, converged "
          f"{out3['converged']} (cycles {out3['cycles']}, "
          f"{time.time() - t0:.0f} s)", flush=True)
    runs["pinned-r1"] = {
        "span": np.asarray(out3["mo"])[:, NCORE:NCORE + ncas_s],
        "roots": out3["roots"], "ncas": ncas_s, "mol": mol_d2h}

    # ---- 1. the span principal angles (the pairwise) ------------------
    print("\n===== A. the span principal angles (the cosines; 1.0 = "
          "the shared direction) =====", flush=True)
    names = list(runs)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = runs[names[i]], runs[names[j]]
            g = a["span"].T @ s1 @ b["span"]
            sv = np.linalg.svd(g, compute_uv=False)
            shared = int(np.sum(sv > 0.99))
            print(f"  {names[i]:16s} vs {names[j]:16s}: cosines "
                  f"[min {sv.min():.4f} / mean {sv.mean():.4f} / max "
                  f"{sv.max():.4f}], the shared directions "
                  f"{shared}/{ncas}", flush=True)

    # ---- 2. the occupation fingerprints -------------------------------
    print("\n===== B. the occupation fingerprints (the ensemble "
          "1-RDM diagonals, equal root weights) =====", flush=True)
    occs = {}
    for name in names:
        r = runs[name]
        if "ci" in r:
            cis = r["ci"]
        else:
            cis = [c for _e, cii in r["roots"].values()
                   for c in block_cis(cii)]
        occs[name] = root_rdm_diag(cis, r["ncas"], *SECTOR)
        print(f"  {name:16s}: "
              + "  ".join(f"{float(o):.2f}" for o in occs[name]),
              flush=True)
    print("  -- the sorted-profile deviations (the gauge-free "
          "comparison) --", flush=True)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            p1 = np.sort(occs[names[i]])[::-1]
            p2 = np.sort(occs[names[j]])[::-1]
            n = min(len(p1), len(p2))
            dev = float(np.max(np.abs(p1[:n] - p2[:n])))
            print(f"  {names[i]:16s} vs {names[j]:16s}: max |d occ| "
                  f"(sorted profiles) {dev:.3f}", flush=True)

    # ---- 3. the <L2> partition walks ----------------------------------
    print("\n===== C. the <L2> partition walks =====", flush=True)
    for name in names:
        r = runs[name]
        expect = meters(r["mol"], r["span"], *SECTOR)
        vals = []
        if "ci" in r:
            for c in r["ci"]:
                vals.append(float(expect("l2", np.asarray(c).ravel())))
            print(f"  {name}:", "  ".join(f"{v:.1f}" for v in vals),
                  flush=True)
        else:
            for irr, (_e, cii) in r["roots"].items():
                cl = block_cis(cii)
                block = [float(expect("l2", np.asarray(c).ravel()))
                         for c in cl]
                vals.extend(block)
                print(f"  {name} [{irr}]: "
                      + "  ".join(f"{v:.1f}" for v in block), flush=True)
        print(f"  {name} the counts: {snap_counts(vals)}", flush=True)

    # ---- 4. the input windows vs the converged spans -------------------
    print("\n===== D. the INPUT windows vs the converged spans (did "
          "the optimizer move the span?) =====", flush=True)
    sv_f = np.linalg.svd(act.T @ s1 @ runs["free-C1"]["span"],
                         compute_uv=False)
    sv_p = np.linalg.svd(act.T @ s1 @ runs["pinned-asbuilt"]["span"],
                         compute_uv=False)
    sv_s = np.linalg.svd(act_s.T @ s1 @ runs["pinned-r1"]["span"],
                         compute_uv=False)
    sv_cross = np.linalg.svd(act_s.T @ s1 @ act, compute_uv=False)
    print(f"  input as-built act vs converged free-C1 span: "
          f"[min {sv_f.min():.4f} / max {sv_f.max():.4f}]", flush=True)
    print(f"  input as-built act vs converged pinned-asbuilt span: "
          f"[min {sv_p.min():.4f} / max {sv_p.max():.4f}]", flush=True)
    print(f"  input @r1 act_s vs converged pinned-r1 span: "
          f"[min {sv_s.min():.4f} / max {sv_s.max():.4f}], the shared "
          f"directions {int(np.sum(sv_s > 0.99))}/{ncas_s}", flush=True)
    print(f"  (the INPUT windows themselves: act_s vs act: "
          f"[min {sv_cross.min():.4f} / max {sv_cross.max():.4f}], "
          f"the shared directions {int(np.sum(sv_cross > 0.99))}/"
          f"{ncas})", flush=True)

    print(f"\n({time.time() - t00:.0f} s total)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
