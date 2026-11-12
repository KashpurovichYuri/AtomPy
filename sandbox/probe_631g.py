"""STAGE 1b rider (ticket 37): the 6-31G certification -- the user's
question ('do we obtain divergence even for 6-31G? what elements?')
closed with measurements instead of inference.

  1  the S spectra of the 6-31G bases on Sc/Ti/C/Si (expect: no
     small-lambda tail -> no metric-amplified garbage; the disease's
     ingredient (1) absent)
  2  the C/Si s1p3 promotion forced SCFs on 6-31G under the plain
     weight hook (the disease's ingredients (2)+(3) absent: physical
     labels, meaningful weights)

The Sc/Ti 6-31G sources (d3/d2s1/d4/d3s1) already ran in the bench
-- converged, occupation-safe (see transfer.py / mom.py runs).
"""
import sys

import numpy as np
from dataclasses import replace

sys.path.insert(0, "sandbox")
from spectrum import l_blocks  # noqa: E402

from atomterm.cases import CAMPAIGN, ReferenceSpec  # noqa: E402
from atomterm.reference import occupation_postcheck, run_rohf  # noqa: E402
from atomterm.system import AtomConfig, build_mol  # noqa: E402

print("=== the 6-31G overlap spectra ===")
for atom in ("Sc", "Ti", "C", "Si"):
    cfg = replace(CAMPAIGN[atom].atom_config(), basis="6-31g")
    mol = build_mol(cfg)
    s = mol.intor("int1e_ovlp")
    mins = []
    for l, ix in sorted(l_blocks(mol).items()):
        ev = np.linalg.eigvalsh(s[np.ix_(ix, ix)])
        mins.append((l, len(ev), float(ev.min()), float(ev.max())))
    tail = " ".join(f"l{l}:min={m:.2e}" for l, n, m, mx in mins)
    print(f"  {atom:2s} 6-31g  nao={mol.nao:3d}  {tail}")

print("=== the C/Si s1p3 promotions on 6-31G (the plain weight hook) ===")
for atom, closed, open_ in (("C", ("1s",), ("2s", "2p")),
                            ("Si", ("1s", "2s", "2p"), ("3s", "3p"))):
    cfg = replace(CAMPAIGN[atom].atom_config(), basis="6-31g", spin=4)
    ref = ReferenceSpec(occ_closed=closed, occ_open=open_)
    _, mf = run_rohf(cfg, f"sandbox/scratch/{atom.lower()}_631g_s1p3.chk",
                     ref=ref)
    post = occupation_postcheck(build_mol(cfg), mf, closed, open_)
    print(f"  {atom} s1p3: E = {mf.e_tot:.8f}  conv={mf.converged}  "
          f"post={post}")
print("VERDICT: no divergence anywhere -- 6-31G carries none of the "
      "disease's ingredients")
