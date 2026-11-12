"""STAGE 1b probe 2 (ticket 37): WHY does the character-of-occupied
post-check fail runs whose energies/occupations are right?  Inspect
the converged MOM run's occupied columns -- weights, norms -- and the
span selections on the unc basis (label semantics)."""
import pathlib
import sys

import numpy as np
from dataclasses import replace

sys.path.insert(0, "sandbox")
from project import project_orbitals  # noqa: E402
from mom import seeded_run  # noqa: E402

from atomterm.cases import CAMPAIGN, ReferenceSpec  # noqa: E402
from atomterm.reference import run_rohf  # noqa: E402
from atomterm.system import build_mol  # noqa: E402
from atomterm.windows import select_span, span_projector  # noqa: E402

key, contracted, tag = "Sc-unc", "Sc", "d3"
spec_u = CAMPAIGN[key]
member = next(s for s in CAMPAIGN[contracted].series if s.tag == tag)
mol_u = build_mol(replace(spec_u.atom_config(), spin=member.spin))

labels = mol_u.ao_labels()
d_idx = [i for i, l in enumerate(labels) if l.endswith("d")]
print("the d-family AO labels:", sorted({labels[i][2:] if len(labels[i]) > 2
                                         else labels[i] for i in d_idx}))
print("select_span('3d'):", len(select_span(mol_u, "3d")), "AOs;",
      "select_span('4s'):", len(select_span(mol_u, "4s")))

cfg_s = replace(CAMPAIGN[contracted].atom_config(), basis="6-31g",
                spin=member.spin)
chk = pathlib.Path(f"sandbox/scratch/{key.lower()}_6-31g_{tag}.chk")
chk.unlink(missing_ok=True)
mol_s, mf_s = run_rohf(cfg_s, str(chk),
                       ref=ReferenceSpec(occ_closed=member.closed,
                                         occ_open=member.open))
occ = np.asarray(mf_s.mo_occ)
c_b, rep = project_orbitals(mol_u, mol_s, mf_s.mo_coeff[:, occ > 0],
                            thresh=1e-8)
mf = seeded_run(mol_u, c_b, occ[occ > 0], accel="diis", max_cycle=80)

s = mol_u.intor_symmetric("int1e_ovlp")
p_closed = span_projector(s, sum((select_span(mol_u, p)
                                   for p in member.closed), []))
p_open = span_projector(s, sum((select_span(mol_u, p)
                                for p in member.open), []))
occ_r = np.asarray(mf.mo_occ)
print(f"\nE = {mf.e_tot:.8f} conv={mf.converged}")
for kind, w_proj in (("closed", p_closed), ("open", p_open)):
    idx = np.flatnonzero(occ_r == (2.0 if kind == "closed" else 1.0))
    ws = [float(np.real(mf.mo_coeff[:, i].conj() @ w_proj
                        @ mf.mo_coeff[:, i])) for i in idx]
    ns = [float(np.linalg.norm(mf.mo_coeff[:, i])) for i in idx]
    print(f"{kind}: n={len(idx)} weights="
          f"{np.array2string(np.array(ws), precision=3)}")
    print(f"       norms={np.array2string(np.array(ns), precision=2)}")
