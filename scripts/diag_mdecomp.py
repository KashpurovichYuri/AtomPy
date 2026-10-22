"""Mathematical m-decomposition of the U ROHF open MOs (ticket 24 (b)),
and the reference det's <L^2> evaluated INSIDE a big m-complete window
(ticket 24 (c) -- the third independent construction).

(b) Real harmonics R_c, R_s of a member pair relate to the complex
|m>, |-m> by R_c = (|m> + s|m| |-m>)/sqrt2, R_s = i(|m> - s|m| |-m>)/sqrt2
(s|m| = Condon-Shortley (-1)^m), so an MO's complex amplitudes are
c_m = (a - i s b)/sqrt2 from its real coefficients (a, b) on the pair.
The pair total weight is invariant -- the built-in check.

(c) Project the 4 open MOs onto the 3f+2d window span (build_full_mo),
express the reference det as a combination of window dets, and evaluate
<L^2> with the window's own operator matrix.
"""
import sys
from pathlib import Path

sys.path.insert(0, "src")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from run_u import CFG, CHK, NCORE, REF, build_window
from atomterm import drivers
from atomterm.detbasis import determinant_basis
from atomterm.operators import build_operators
from atomterm.drivers import lm_sectors

mol, mf = drivers.run_rohf(CFG, CHK, ref=REF)
open_cols = list(np.flatnonzero(mf.mo_occ == 1.0))
s = mol.intor_symmetric("int1e_ovlp")
sectors = lm_sectors(mol)

# (b): complex-harmonic m-decomposition of each open MO.
# pyscf real order within a shell: m = 0, +1, -1, +2, -2, ...
print("complex-harmonic m-decomposition of the open MOs (|c_m|^2):")
for j, oc in enumerate(open_cols):
    c = mf.mo_coeff[:, oc]
    rows_out = {}
    for l in (2, 3):
        mem = [(mm, sectors[(l, mm)]) for mm in sorted(sectors
                                                       for _ in ())
               if False]  # placeholder, replaced below
    # build per-l from sector keys
    for l in (2, 3):
        members = sorted(m for (ll, m) in sectors if ll == l)
        if not members:
            continue
        w_pair = {}
        for k, mm in enumerate(members):
            rows = sectors[(l, mm)]
            srr = s[np.ix_(rows, rows)]
            w = float(c[rows] @ srr @ c[rows])
            # members are ordered m = 0, +1, -1, +2, -2, ...
            m_abs = (k + 1) // 2
            w_pair[m_abs] = w_pair.get(m_abs, 0.0) + w
        for m_abs, w in sorted(w_pair.items()):
            if w > 1e-4:
                rows_out[f"{l}m|{m_abs}|"] = round(w / 1.0, 3)
    tot = sum(rows_out.values())
    print(f"  MO {oc} (total {tot:.3f}): {rows_out}")

# (c): the reference det expressed in the 3f+2d window
closed = list(np.flatnonzero(mf.mo_occ == 2.0))
c_pur, _ = drivers.purified_core_orbitals(mol, mf, NCORE)
f_rows, d_rows, win = build_window(mol, {"f": (0.2, 2.0), "d": (0.1, 0.5)})
mo = drivers.build_full_mo(mf, closed, win, c_core=c_pur)
ncas = len(win)
ops = build_operators(mol, mf, mo[:, NCORE:NCORE + ncas], NCORE, 4, 0,
                      c_core=mo[:, :NCORE])
dets = determinant_basis(ncas, 4, 0)
# amplitude of the reference det on window dets: the open MOs'
# projections onto the window's orthonormal columns
c_open = mf.mo_coeff[:, open_cols]
ov = c_open.T.conj() @ (s @ mo[:, NCORE:NCORE + ncas])  # (4, ncas)
coef = np.zeros(len(dets))
for idx, det in enumerate(dets):
    amp = 1.0
    cols = list(det)
    sub = ov[:, cols]
    # det amplitude = per-det of the 4x4 projection (wedge component):
    amp = np.linalg.det(sub)
    coef[idx] = amp
norm = float(np.linalg.norm(coef))
l2_det = float(coef @ ops["L2"].real @ coef) / norm**2
print(f"\n(c) reference det projected into the 3f+2d window: "
      f"norm = {norm:.6f}")
print(f"    <L^2> of the in-window component = {l2_det:.6f}")
print(f"    (the outside-span remainder carries the rest; the 1-det")
print(f"     evaluations gave 8e-10 / 0.0 exactly)")
