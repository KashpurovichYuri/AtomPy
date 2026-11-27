"""The reference-comparison V&V leg -- ROHF vs UHF, plus the
spin-broken-start SA question.  THE NEW-PHYSICS LEG: everything is
printed/reported; nothing gates.  Also drops the table's two
REFERENCE rows (ROHF -- the zero of every dE -- and UHF).

  (e) the UHF sibling's spin honesty: E beside the ROHF rung, <S^2>
      by the OVERLAP FORMULA (written out below -- six lines, the
      analysis of interest), the corresponding-orbital pairing
      spectrum;
  the UHF-seeded SA: the visible free-route construction riding the
  UHF object (PySCF's CASSCF accepts it directly) -- does the SA
  repair recover <L2>/<S2> from a RAW spin-broken start?  Measured
  and reported (the committed orchestrator record holds the
  CONVERGED-start reading).

The UHF rung: loaded from ``rungs/uhf.chk`` when present (trusted
cache), else warmed through the production runner (``run_uhf`` with
the basin gate -- one marked call, disclosed).

Run from the repo root:
    OMP_NUM_THREADS=16 python scripts/u14e/u14e_vv_uhf.py \
        |  tee logs/u14e_vv_uhf.log
"""

import sys
import time
from pathlib import Path

import numpy as np

from vv.bands import L_L2, S_S2
from vv.meters import meters
from vv.reference import RUNGS, load_reference, load_rung
from vv.table import EH_EV, drop_row, std_row
from vv.verdict import Verdict

SECTOR = (4, 0)
NCORE = 5
SHELLS = (("10f", 7), ("13d", 5))


def uhf_rung():
    """The UHF sibling: the trusted chk, else the production runner
    (marked call -- the forced-occupation hooks + the basin gate are
    the reference layer's discipline, not under test here)."""
    chk = str(RUNGS / "uhf.chk")
    if Path(chk).exists():
        mol, mf = load_rung(chk)
        print(f"  UHF (trusted rung, loaded): E = {mf.e_tot:.8f}",
              flush=True)
        return mol, mf
    from atomterm.reference import run_uhf
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    mol, mf = run_uhf(spec.atom_config(), chk, ref=spec.ref)
    print(f"  UHF (production runner, marked): E = {mf.e_tot:.8f}",
          flush=True)
    return mol, mf


def contamination(mf):
    """<S^2> by the overlap formula -- written out: Sz(Sz+1) + N_beta
    - sum_ij |<phi_a_i|phi_b_j>|^2 over the occupied blocks; a
    spin-paired determinant makes the sum N_beta and <S^2> = S(S+1)
    exactly; any excess is spin contamination.  Returns
    (s2, sigma) -- sigma the occupied-overlap singular values
    (1 = a paired couple, 0 = fully broken)."""
    ca, cb = mf.mo_coeff
    occ = np.asarray(mf.mo_occ)
    a_idx = np.flatnonzero(occ[0] > 0)
    b_idx = np.flatnonzero(occ[1] > 0)
    na, nb = len(a_idx), len(b_idx)
    x = ca[:, a_idx].conj().T @ mf.get_ovlp() @ cb[:, b_idx]
    sz = (na - nb) / 2.0
    s2 = sz * (sz + 1.0) + nb - float(np.sum(np.abs(x) ** 2))
    return s2, np.linalg.svd(x, compute_uv=False)


def visible_window(mol, mf_rohf, c_core):
    """The same visible window as the 4_12 leg (band-pure shells)."""
    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    sel = [i for i, lab in enumerate(labels)
           if any(lab.startswith(w) for w, _ in SHELLS)]
    if len(sel) != 12:
        sys.exit(f"the label selection found {len(sel)} AOs, "
                 "expected 12")
    nao = mol.nao
    S = mf_rohf.get_ovlp()
    act = np.zeros((nao, len(sel)))
    for k, i in enumerate(sel):
        act[i, k] = 1.0
    P_core = c_core @ c_core.T @ S   # the metric projector
    w = act - P_core @ act
    g = w.T @ S @ w
    gv, U = np.linalg.eigh(g)
    c_win = w @ U @ np.diag(gv ** -0.5) @ U.T
    C = np.hstack([c_core, c_win])
    P_all = C @ C.T @ S
    rest = mf_rohf.mo_coeff - P_all @ mf_rohf.mo_coeff
    g2 = rest.T @ S @ rest
    gv2, U2 = np.linalg.eigh(g2)
    keep = gv2 > 1e-7 * gv2.max()
    c_virt = rest @ U2[:, keep] @ np.diag(gv2[keep] ** -0.5)
    if c_virt.shape[1] != nao - C.shape[1]:
        sys.exit(f"the complement came out {c_virt.shape[1]} wide, "
                 f"expected {nao - C.shape[1]}")
    return np.hstack([C, c_virt])


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    v = Verdict()
    mol, mf_rohf = load_reference()
    e_rohf = float(mf_rohf.e_tot)
    drop_row(std_row("rohf", "-", "-", e_rohf, None, e_rohf,
                     "-", "-"))
    mol_u, mf_uhf = uhf_rung()
    e_gap = float(mf_uhf.e_tot) - e_rohf
    s2, sigma = contamination(mf_uhf)
    n_broken = int(np.sum(sigma < 1.0 - 1e-6))
    v.note("(e) the UHF sibling's spin honesty",
           f"<S2> = {s2:.4f} (Sz exact by construction; the exact "
           f"5-fold would read {S_S2:.3f}); pairing spectrum: "
           f"{len(sigma)} occupied pairs, {n_broken} below "
           f"1 - 1e-6, smallest {np.sort(sigma)[:3]}")
    v.note("(e) the ROHF-vs-UHF reference comparison",
           f"UHF {e_gap:+.4f} mEh ({e_gap * EH_EV:+.4f} eV) below "
           f"the ROHF basin; the cores assume one closed-shell "
           f"density -- the UHF object is a reference-class "
           f"product, not a cells axis")
    drop_row(std_row("uhf", "-", "-", float(mf_uhf.e_tot), None,
                     float(mf_uhf.e_tot), "-", f"{s2:.4f} (det)"))

    # the UHF-seeded SA: the visible free-route construction riding
    # the UHF object -- from a RAW start, does the SA repair recover
    # the term identity?  free_kernel = the family's faithful two-pass
    # (the readout seeded from the SA kernel's own vectors) + the conv
    # capture (track_iters)
    from atomterm.core import purified_core_orbitals
    from vv.routes import free_kernel
    c_core = purified_core_orbitals(mol, mf_rohf, NCORE)[0]
    mo_full = visible_window(mol, mf_rohf, c_core)
    t0 = time.time()
    mc, cas = free_kernel(mf_uhf, 12, SECTOR, mo_full,
                          weights=[1.0 / 17] * 17, nroots=17,
                          track_iters=True)
    print(f"  [sa] UHF-seeded ensemble {mc.e_tot:.8f}, converged "
          f"{mc.converged}, conv {mc.final_de:.1e}/{mc.n_macro} "
          f"({time.time() - t0:.1f} s)", flush=True)
    expect = meters(mol_u, cas.mo_coeff[:, NCORE:NCORE + 12],
                    *SECTOR)
    l2s = [expect("l2", cas.ci[k]) for k in range(17)]
    dev = max(abs(x - L_L2) for x in l2s)
    if dev > 1.0:
        v.note("the UHF-seeded SA from the RAW start descends broken "
               "(this draw)",
               f"<L2> dev {dev:.2e} (<L2> in [{min(l2s):.4f}, "
               f"{max(l2s):.4f}], ensemble {mc.e_tot:.4f}, conv "
               f"{mc.final_de:.1e}/{mc.n_macro}"
               f"{'' if mc.converged else '*'} -- an "
               f"energy of the UHF-CASCI partition, not comparable "
               f"to the ROHF-route optimum)")
    else:
        v.note("the UHF-seeded SA recovers <L2> = 72 from the RAW "
               "start (this draw)",
               f"17 roots, max dev {dev:.2e}, ensemble {mc.e_tot:.8f}, "
               f"conv {mc.final_de:.1e}/{mc.n_macro}"
               f"{'' if mc.converged else '*'} "
               f"= the ROHF-route optimum.  DRAW-DEPENDENT class: "
               f"the family's earlier draw descended broken (<L2> "
               f"dev 7.2e+01, that log in git history at the "
               f"824b5d8 commit), and the committed converged-start "
               f"reading holds regardless (dev 8.9e-15) -- the SA "
               f"repair holds the converged gauge; the raw "
               f"spin-broken start is a coin toss")
    ok = v.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
