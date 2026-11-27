"""The reference-comparison leg -- the table's two REFERENCE rows (ROHF,
the zero of every dE; UHF) and the UHF sibling's spin honesty, plus
the UHF-seeded SA question (a note; nothing here gates).

  (e) the UHF sibling's spin honesty: E beside the ROHF rung, <S²>
      by the OVERLAP FORMULA (written out below -- six lines, the
      analysis of interest), the corresponding-orbital pairing
      spectrum;
  the UHF-seeded SA: the free SA (``vv.routes.free_sa``) riding the
      UHF object -- PySCF's CASSCF takes a UHF mean-field spin-free
      on given orbitals -- from a RAW spin-broken start, does the SA
      repair recover the term identity?  DRAW-DEPENDENT: the family's
      earlier draw descended broken, a later one recovered (both
      readings are in the committed logs); reported, never asserted.

The UHF rung: loaded from the trusted chk when present, else warmed
through the production runner (``run_uhf`` with the basin gate -- one
marked call, disclosed).

Run from the repo root:
    OMP_NUM_THREADS=16 python scripts/u14e/u14e_vv_uhf.py \\
        |  tee logs/u14e_vv_uhf.log
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from vv import cases as ladder  # noqa: E402
from vv.reference import Substrate, load_rung, rung_path  # noqa: E402
from vv.routes import free_sa  # noqa: E402
from vv.table import EH_EV, Route, Row, drop_row  # noqa: E402
from vv.verdict import Verdict  # noqa: E402


def uhf_rung(sub: Substrate) -> Substrate:
    """The UHF sibling: the trusted chk, else the production runner
    (a marked call -- the forced-occupation hooks and the basin gate
    are the reference layer's discipline, not under test here)."""
    chk = rung_path("U14e", "uhf.chk")
    if chk.exists():
        mol, mf = load_rung(chk)
        print(f"  UHF (trusted rung, loaded): E = {mf.e_tot:.8f}",
              flush=True)
        return Substrate(mol, mf, sub.spec)
    from atomterm.reference import run_uhf
    mol, mf = run_uhf(sub.spec.atom_config(), str(chk), ref=sub.spec.ref)
    print(f"  UHF (production runner, marked): E = {mf.e_tot:.8f}",
          flush=True)
    return Substrate(mol, mf, sub.spec)


def contamination(mf) -> tuple[float, np.ndarray]:
    """<S²> by the overlap formula:
    Sz(Sz+1) + N_beta − Σ_ij |<φ_i^α|φ_j^β>|² over the occupied blocks.
    A spin-paired determinant makes the sum N_beta and <S²> = S(S+1)
    exactly; any excess is spin contamination.

    Returns
    -------
    (s2, sigma) -- ``sigma`` the occupied-overlap singular values
    (1 = a paired couple, 0 = fully broken).
    """
    ca, cb = mf.mo_coeff
    occ = np.asarray(mf.mo_occ)
    a_idx = np.flatnonzero(occ[0] > 0)
    b_idx = np.flatnonzero(occ[1] > 0)
    na, nb = len(a_idx), len(b_idx)
    x = ca[:, a_idx].conj().T @ mf.get_ovlp() @ cb[:, b_idx]
    sz = (na - nb) / 2.0
    s2 = sz * (sz + 1.0) + nb - float(np.sum(np.abs(x) ** 2))
    return s2, np.linalg.svd(x, compute_uv=False)


def main(argv=None) -> int:
    argparse.ArgumentParser(
        description="the U14e V&V reference leg: the two reference "
                    "rows + the UHF sibling's spin honesty (notes only)"
    ).parse_args(sys.argv[1:] if argv is None else argv)
    v = Verdict()
    sub = Substrate.load()
    e_rohf = float(sub.mf.e_tot)
    drop_row(Row("rohf", Route.REFERENCE, "-", e_rohf, None, "-", "-"))

    uhf = uhf_rung(sub)
    e_uhf = float(uhf.mf.e_tot)
    e_gap = e_uhf - e_rohf
    s2, sigma = contamination(uhf.mf)
    n_broken = int(np.sum(sigma < 1.0 - 1e-6))
    v.note("(e) the UHF sibling's spin honesty",
           f"<S2> = {s2:.4f} (Sz exact by construction; the exact "
           f"5-fold would read {ladder.S2_QUINTET:.3f}); pairing "
           f"spectrum: {len(sigma)} occupied pairs, {n_broken} below "
           f"1 - 1e-6, smallest {np.sort(sigma)[:3]}")
    v.note("(e) the ROHF-vs-UHF reference comparison",
           f"UHF {e_gap * 1e3:+.4f} mEh ({e_gap * EH_EV:+.4f} eV) below "
           f"the ROHF basin; the cores assume one closed-shell "
           f"density -- the UHF object is a reference-class product, "
           f"not a cells axis")
    drop_row(Row("uhf", Route.REFERENCE, "-", e_uhf, None, "-",
                 f"{s2:.4f} (det)"))

    # the UHF-seeded SA: the free route riding the UHF object, on the
    # (4,12) visible window built from the ROHF reference
    space = ladder.visible_window(sub, ladder.SHELLS_412, ladder.SECTOR_4).space
    sol = free_sa(uhf, space, ladder.NROOTS)
    meters = sol.meters()
    dev = max(abs(l2 - ladder.L2_TERM) for l2, _ in meters)
    cls = "recovers <L2> = 72 from the RAW start" if dev <= 1.0 else \
          "descends broken from the RAW start"
    v.note(f"the UHF-seeded SA {cls} (this draw; a coin toss across draws)",
           f"{len(meters)} roots, <L2> in "
           f"[{min(l2 for l2, _ in meters):.4f}, {max(l2 for l2, _ in meters):.4f}] "
           f"(max dev {dev:.2e}), <S2> up to {max(s for _, s in meters):.4f}; "
           f"ensemble {sol.ensemble:.8f}, conv {sol.conv}; certificate "
           f"{sol.certificate}")
    v.dump()          # a notes-only leg: nothing here gates
    return 0


if __name__ == "__main__":
    sys.exit(main())
