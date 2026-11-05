"""One-shot probe (removed after use): the U14e shell-band identities.

1. the shell-slot inventory (the labels the window/fold tags must name);
2. per-orbital shell decomposition (Mulliken gross populations over the
   shell spans) for the closed core, the open shells, and the lowest
   virtuals -- which bands carry 6s/7s/6p and the valence f/d;
3. CASCI(4,12) ranking over the candidate 1f+1d AO pairs (the ROHF
   closed columns as the frozen core) -- the pair the minimal lane
   window declares.

Results -> elements/U14e/logs/probe_bands.log (provenance).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pyscf import mcscf  # noqa: E402

from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.system import build_mol  # noqa: E402
from atomterm.reference import run_rohf  # noqa: E402

import re  # noqa: E402

cfg = CAMPAIGN["U14e"]
mol = build_mol(cfg.atom_config())
_, mf = run_rohf(cfg.atom_config(), "elements/U14e/rungs/rohf.chk",
                 ref=cfg.ref)
S = mol.intor("int1e_ovlp")
labels = [line.split()[-1] for line in mol.ao_labels()]

# shell spans: group AO rows by the m-suffix-stripped label
_M_SUFFIX = re.compile(r"(x2-y2|xy|xz|yz|z\^2|x|y|z|[+-]\d+)$")


def base(lbl):
    return _M_SUFFIX.sub("", lbl)


shells: dict[str, list[int]] = {}
for i, lbl in enumerate(labels):
    shells.setdefault(base(lbl), []).append(i)
print(f"n AO {mol.nao}  n shells {len(shells)}")
by_l = {}
for sh in shells:
    by_l.setdefault(sh[-1], []).append(sh)
for l, names in by_l.items():
    print(f"  {l}: {names}")

# Mulliken gross population of orbital c on a shell span:
# w = c_band @ (S @ c)[band]  (sums to 1 across all bands)
def weights(c):
    Sc = S @ c
    return {sh: float(c[idx] @ Sc[idx]) for sh, idx in shells.items()}

def show(name, c):
    w = weights(c)
    top = sorted(w.items(), key=lambda kv: -kv[1])[:3]
    tops = "  ".join(f"{n}:{v:.3f}" for n, v in top if v > 0.005)
    print(f"  {name:12s} eps={mf.mo_energy[i]: .4f}  {tops}")

occ_cols = [i for i in range(mol.nao) if mf.mo_occ[i] == 2]
open_cols = [i for i in range(mol.nao) if 0 < mf.mo_occ[i] < 2]
print(f"\nclosed columns {occ_cols}  open columns {open_cols}")
mo = mf.mo_coeff
print("\n-- closed core (the fold) --")
for i in occ_cols:
    show(f"mo[{i}]", mo[:, i])
print("-- open shells --")
for i in open_cols:
    show(f"mo[{i}]", mo[:, i])
print("-- lowest virtuals (by eps) --")
virt = sorted((i for i in range(mol.nao) if mf.mo_occ[i] == 0),
              key=lambda i: mf.mo_energy[i])
for i in virt[:15]:
    show(f"mo[{i}]", mo[:, i])

# CASCI(4,12) ranking over candidate 1f+1d pairs -- through the
# PROJECT's own anchor construction (run_cell: the same fold, E_core,
# closure meters the declared anchors carry; the raw PySCF-CASCI
# permutation probe mis-decoded the ROHF columns and was discarded)
from atomterm.cells import CellKey, run_cell  # noqa: E402

core_idx = occ_cols
print("\n-- run_cell ranking (anchor construction), 1f+1d pairs --")
import os  # noqa: E402

cand_f = ["8f", "9f", "10f", "11f"]
cand_d = ["12d", "13d"]
for f in cand_f:
    for d in cand_d:
        key = CellKey("u14e", f"{f}{d}", (4, 0), method="casci")
        facts = run_cell(key, mf, workdir="elements")
        st0 = facts.analysis.states[0]
        print(f"  {key.address}: E = {st0.energy:.6f}  "
              f"<L2> = {st0.l2:.4f}  term = {st0.term}  "
              f"|[P,Lz]| = {facts.win_lz:.1e}  "
              f"|[P,Lp]| = {facts.win_lperp:.1e}")
        os.unlink(Path("elements") / "U14e" / "rungs" / f"{key.address}.chk")
print("\nscratch chks removed")
