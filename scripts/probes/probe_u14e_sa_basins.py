"""The (6,13)@r1 SA/CASSCF basin probe -- the valley map, printed
only (NO store writes; the user evaluates before any table/leg
change).

THE QUESTION (the user's review of the sa-casscf (free) block): the
single-root CASSCF descends the QUINTET valley (-50.4623, <L2> 72.08,
<S2> 6) while the free SA descends the SEPTET valley (-50.4529,
<L2> 90, <S2> 12 -- the exchange-maximal six-parallel class).  Which
wins -- the ensemble OBJECTIVE or the SEED?  Five runs:

  A  CASSCF from the input window      the quintet valley record
  B  free SA x17 from the input window the septet landing + the
                                       occupation FINGERPRINT of the
                                       L=9 root (the prediction:
                                       5f3(m=3,2,1) 6d2(m=2,1) 7s1(m=0),
                                       Sum m = 9; the alternative f3d3
                                       reads the s column ~0)
  C  free SA from A's orbitals         objective vs seed
  D  CASSCF from B's orbitals          does the single root escape to
                                       the quintet valley?
  E1 the pinned SA (the term's own     the 5L-block orbitals
     wfnsym-anchored ensemble)
  E2 free SA from E1's orbitals        the gauge-anchoring test: does
                                       the free SA STAY on the 5L
                                       block when started from the
                                       term's orbitals?

Every run is a certified solution of ``vv.routes`` (a restart rides
``ActiveSpace.on``: the same partition on a converged solution's
orbitals) and prints: E, <L2>/<S2> of the ground root, the
converged-window CLOSURE meters (the (a) measurement: the input
windows close at 1e-13; the defect is created by the optimization --
here measured, not inferred), and the active OCCUPATIONS (the class
fingerprint).

Run from the repo root:
    OMP_NUM_THREADS=64 python scripts/probes/probe_u14e_sa_basins.py \\
        | tee logs/probe_u14e_sa_basins.log
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "u14e"))

import numpy as np  # noqa: E402
from pyscf.fci import direct_spin1  # noqa: E402
from pyscf.mcscf import addons as mcaddons  # noqa: E402

from atomterm.operators import window_closure  # noqa: E402

from vv.cases import ADDR, NROOTS, TERM  # noqa: E402
from vv.reference import Substrate  # noqa: E402
from vv.routes import free_casscf, free_sa, pinned_sa  # noqa: E402
from vv.windows import s_column, window_6_13_r1  # noqa: E402


def _occupations(sol, root=0):
    """The active occupations of one root -- the class fingerprint."""
    st = sol.states[root].state
    return np.diag(direct_spin1.make_rdm1(st.ci, st.ncas, st.nelecas)).real


def _occ_print(tag, sol, root=0):
    occ = _occupations(sol, root)
    print(f"  [{tag}] occupations:", flush=True)
    print("    " + "  ".join(f"{i}:{o:.2f}" for i, o in enumerate(occ)),
          flush=True)
    return occ


def _ground_report(tag, sol, root=0):
    """<L2>/<S2> of one root + the converged-window closure."""
    st = sol.states[root].state
    (l2, s2) = sol.meters()[root]
    cl = window_closure(st.mol, st.active)
    ss, mult = mcaddons.spin_square(st.mc, ci=st.ci)
    print(f"  [{tag}] ground root: <L2> {l2:.6f}, <S2> {s2:.6f} "
          f"(spin_square: {ss:.6f}, {mult:.1f}); the converged "
          f"window's closure: ||[P,Lz]|| {cl['lz']:.2e}, "
          f"||[P,Lperp]|| {cl['lperp']:.2e}", flush=True)
    return l2, s2, cl


def _run(tag, go, what):
    t0 = time.time()
    sol = go()
    energy = (f"E {sol.states[0].e_ref:.8f}" if sol.ensemble is None
              else f"ensemble {sol.ensemble:.8f}")
    print(f"  [{tag}] {energy}, conv {sol.conv} ({time.time() - t0:.0f} s); "
          f"certificate {sol.certificate} -- {what}", flush=True)
    return sol


def main() -> int:
    t00 = time.time()
    sub = Substrate.load()
    built = window_6_13_r1(sub, ADDR["(6,13)@r1"])
    space = built.space
    for line in built.report:
        print(f"  {line}", flush=True)
    cl_in = space.closure(sub.mol)
    j_s = s_column(space.active, sub.mol)
    print(f"  [window] the input closure ||[P,Lz]|| {cl_in['lz']:.2e}, "
          f"||[P,Lperp]|| {cl_in['lperp']:.2e}; the s column at {j_s}",
          flush=True)

    print("===== A: CASSCF, input window =====", flush=True)
    a = _run("A", lambda: free_casscf(sub, space), "the input seed")
    _ground_report("A", a)
    _occ_print("A", a)

    print("===== B: free SA x17, input window =====", flush=True)
    b = _run("B", lambda: free_sa(sub, space, NROOTS), "the input seed")
    _ground_report("B", b)
    occ_b = _occ_print("B", b)
    print(f"  [B] the fingerprint read: s column occ "
          f"{occ_b[j_s]:.3f}; f block sum {occ_b[:7].sum():.2f}; "
          f"d block sum {occ_b[7:12].sum():.2f}", flush=True)

    print("===== C: free SA x17, started from A's orbitals =====", flush=True)
    c = _run("C", lambda: free_sa(sub, space.on(a.states[0].mo_coeff, sub.mf), NROOTS),
             "A's seed")
    _ground_report("C", c)

    print("===== D: CASSCF, started from B's orbitals =====", flush=True)
    d = _run("D", lambda: free_casscf(sub, space.on(b.states[0].mo_coeff, sub.mf)),
             "B's seed")
    _ground_report("D", d)

    print("===== E1: the pinned SA on the @r1 window =====", flush=True)
    e1 = _run("E1", lambda: pinned_sa(sub, space, TERM),
              "the term's own wfnsym-anchored ensemble")
    l2_e1 = e1.meters()
    print(f"  [E1] <L2> [{min(l2 for l2, _ in l2_e1):.6f}, "
          f"{max(l2 for l2, _ in l2_e1):.6f}], <S2> max "
          f"{max(s2 for _, s2 in l2_e1):.6f}", flush=True)

    print("===== E2: free SA x17, started from E1's orbitals =====", flush=True)
    e2 = _run("E2", lambda: free_sa(sub, space.on(e1.states[0].mo_coeff, sub.mf), NROOTS),
              "E1's seed")
    _ground_report("E2", e2)
    _occ_print("E2", e2)

    # ---- the summary ---------------------------------------------------
    print("\n===== the valley map =====", flush=True)
    print(f"  A  casscf   (input seed)      E {a.states[0].e_ref:.8f}", flush=True)
    print(f"  D  casscf   (B's seed)        E {d.states[0].e_ref:.8f}", flush=True)
    print(f"  B  sa x17   (input seed)      ensemble {b.ensemble:.8f}", flush=True)
    print(f"  C  sa x17   (A's seed)        ensemble {c.ensemble:.8f}", flush=True)
    print(f"  E1 pinned sa                  ensemble {e1.ensemble:.8f}", flush=True)
    print(f"  E2 sa x17   (E1's seed)       ensemble {e2.ensemble:.8f}", flush=True)
    print(f"({time.time() - t00:.0f} s total)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
