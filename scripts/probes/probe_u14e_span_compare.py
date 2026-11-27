"""The (6,13) active-space identity probe -- do the routes and the
windows reach the SAME active space?  The energies cannot answer on a
quasi-degenerate manifold (the class spread ~0.1 mEh degenerates every
partition); the identity instruments:

  1. the SPAN PRINCIPAL ANGLES -- the converged active subspaces'
     cross-Gram singular values under the AO overlap (the gauge-
     invariant subspace comparison: identical spans -> all cosines
     1.0);  the C1 and D2h substrates share the AO basis and overlap
     (certified where the wrapper is built: ``vv.reference``);
  2. the OCCUPATION FINGERPRINTS -- the ensemble-average 1-RDM
     diagonals per active column (equal weight over the 17 roots on
     both sides), compared as sorted profiles (the column gauges
     differ across routes; the profile is the gauge-free witness);
  3. the <L2> PARTITION WALKS -- the per-root <L2> per run, snapped
     to the L(L+1) grid (the term-partition identity: the free's
     L=9 partition vs the pinned's counts).

THE THREE RUNS (the two-axis comparison), each a certified solution
of ``vv.routes`` (its states carry the orbitals their CI vectors were
solved in):
  the route axis at the fixed as-built window:
     free SA (C1)  vs  pinned SA (pinned@pur);
  the window axis at the fixed pinned@pur route:
     as-built window  vs  @r1 window (``vv.windows.window_6_13_r1``,
     the mixed-source s channel).

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

from vv.reference import Substrate  # noqa: E402
from vv.routes import free_sa, pinned_sa  # noqa: E402
from vv.windows import window_6_13_r1, window_from_key  # noqa: E402

SECTOR = (5, 1)
ROOTS = 17
TERM = (8, "u")
ADDR_613 = "u14e_10f13d7s_5-1_sa-casscf_ao_-_pur_avg"
ADDR_613_R1 = "u14e_5f6d7s@r1_5-1_casci_ao_-_pur_avg"


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


def _run(name, go, t00) -> dict:
    t0 = time.time()
    sol = go()
    print(f"  E {sol.ensemble:.8f}, conv {sol.conv} ({time.time() - t0:.0f} s); "
          f"certificate {sol.certificate}", flush=True)
    st = sol.states[0]
    return {"sol": sol, "span": st.active, "ncas": st.ncas,
            "mol": st.mol, "name": name}


def main() -> int:
    t00 = time.time()
    sub = Substrate.load()
    s1 = sub.mf.get_ovlp()
    built_a = window_from_key(sub, ADDR_613)       # the as-built window
    built_r = window_6_13_r1(sub, ADDR_613_R1)     # the @r1 window
    act, act_s = built_a.space.active, built_r.space.active
    ncas, ncas_s = built_a.space.ncas, built_r.space.ncas

    runs = {}
    print("===== 1. the free SA (C1, the as-built window) =====", flush=True)
    runs["free-C1"] = _run("free-C1", lambda: free_sa(sub, built_a.space, ROOTS), t00)

    print("===== 2. the pinned SA (pinned@pur, the as-built window) =====",
          flush=True)
    runs["pinned-asbuilt"] = _run("pinned-asbuilt", lambda: pinned_sa(sub, built_a.space, TERM), t00)

    print("===== 3. the pinned SA (pinned@pur, the @r1 window) =====", flush=True)
    runs["pinned-r1"] = _run("pinned-r1", lambda: pinned_sa(sub, built_r.space, TERM), t00)
    names = list(runs)

    # ---- 1. the span principal angles (the pairwise) ------------------
    print("\n===== A. the span principal angles (the cosines; 1.0 = "
          "the shared direction) =====", flush=True)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = runs[names[i]], runs[names[j]]
            sv = np.linalg.svd(a["span"].T @ s1 @ b["span"], compute_uv=False)
            shared = int(np.sum(sv > 0.99))
            print(f"  {names[i]:16s} vs {names[j]:16s}: cosines "
                  f"[min {sv.min():.4f} / mean {sv.mean():.4f} / max "
                  f"{sv.max():.4f}], the shared directions "
                  f"{shared}/{min(a['ncas'], b['ncas'])}", flush=True)

    # ---- 2. the occupation fingerprints -------------------------------
    print("\n===== B. the occupation fingerprints (the ensemble "
          "1-RDM diagonals, equal root weights) =====", flush=True)
    occs = {}
    for name in names:
        r = runs[name]
        occs[name] = root_rdm_diag([st.ci for st in r["sol"].states],
                                   r["ncas"], *SECTOR)
        print(f"  {name:16s}: "
              + "  ".join(f"{float(o):.2f}" for o in occs[name]), flush=True)
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
        sol = runs[name]["sol"]
        meters = sol.meters()
        vals = []
        for irrep in dict.fromkeys(st.irrep for st in sol.states):
            block = [round(l2, 1) for st, (l2, _s2) in zip(sol.states, meters)
                     if st.irrep == irrep]
            vals.extend(block)
            label = name if irrep is None else f"{name} [{irrep}]"
            print(f"  {label}: " + "  ".join(f"{v:.1f}" for v in block),
                  flush=True)
        print(f"  {name} the counts: {snap_counts(vals)}", flush=True)

    # ---- 4. the input windows vs the converged spans -------------------
    print("\n===== D. the INPUT windows vs the converged spans (did "
          "the optimizer move the span?) =====", flush=True)
    for label, a, b in (
            ("input as-built act vs converged free-C1 span", act,
             runs["free-C1"]["span"]),
            ("input as-built act vs converged pinned-asbuilt span", act,
             runs["pinned-asbuilt"]["span"]),
            ("input @r1 act_s vs converged pinned-r1 span", act_s,
             runs["pinned-r1"]["span"]),
            ("the INPUT windows themselves: act_s vs act", act_s, act)):
        sv = np.linalg.svd(a.T @ s1 @ b, compute_uv=False)
        print(f"  {label}: [min {sv.min():.4f} / max {sv.max():.4f}], "
              f"the shared directions {int(np.sum(sv > 0.99))}/"
              f"{min(a.shape[1], b.shape[1])}", flush=True)

    print(f"\n({time.time() - t00:.0f} s total)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
