"""Why is the gated U ROHF determinant a 5S-deg state? (ticket 24)

Three looks at the same fact:
1. m-composition of the 4 open MOs (which m-set did the SCF pick?);
2. <L^2> of the reference det via the SECOND-QUANTIZED construction
   (independent of the Slater-Condon path that measured 8e-10);
3. the ensemble <L^2>: the occupation-number-weighted average of L^2
   over ALL 4-open-electron dets of the open space (the AOC reading --
   is the 5S-deg a property of the specific det or of the ensemble?).
"""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import numpy as np

from run_u import CFG, CHK, NCORE, REF
from atomterm import drivers
from atomterm.detbasis import determinant_basis
from atomterm.operators import build_operators
from atomterm.drivers import lm_sectors

mol, mf = drivers.run_rohf(CFG, CHK, ref=REF)
open_cols = list(np.flatnonzero(mf.mo_occ == 1.0))
c_pur, _ = drivers.purified_core_orbitals(mol, mf, NCORE)

# 1. m-composition of each open MO
sectors = lm_sectors(mol)
s = mol.intor_symmetric("int1e_ovlp")
print("open-MO S-weights per (l, m) sector:")
for j, oc in enumerate(open_cols):
    c = mf.mo_coeff[:, oc]
    parts = {}
    for (l, m), rows in sorted(sectors.items()):
        srr = s[np.ix_(rows, rows)]
        w = float(c[rows] @ srr @ c[rows])
        if w > 1e-4:
            parts[(l, m)] = round(w, 4)
    print(f"  MO {oc}: {parts}")

# 2. second-quantized <L^2> of the reference det (1-det space)
mo = np.hstack([c_pur, mf.mo_coeff[:, open_cols]])
ops = build_operators(mol, mf, mo[:, NCORE:], NCORE, 4, 0,
                      c_core=mo[:, :NCORE], method="secondquant")
ev, vec = np.linalg.eigh(ops["H"].real)
l2 = float(vec[:, 0] @ ops["L2"] @ vec[:, 0])
l2m = np.linalg.eigvalsh(ops["L2"].real)
print(f"secondquant: E = {ev[0]:.6f}  <L^2>_ground = {l2:.3e}"
      f"  L2 eigenvalues in space: {np.round(l2m, 6)}")

# 3. ensemble average: the occupation (AOC) density's <L^2>
nf = np.identity(4)
print("ensemble diag of L^2 over the open space:",
      np.round(np.diag(ops["L2"].real), 4))
