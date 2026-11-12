"""STAGE 1b probe (ticket 37): is F[dm_trusted] genuinely pathological,
and what exactly do the MOM runs + the post-check disagree on?

  1  the trusted record's stored mo_energy range vs eig(F[dm], S)
  2  ||dm_trusted||_2 (max eigenvalue) and the veff eigen-range
  3  the P1 legs with FULL reporting (E/conv/post per accel)
  4  the post-check's assignment vs MOM's on the converged MOs
"""
import sys

import numpy as np
from dataclasses import replace

sys.path.insert(0, "sandbox")
from project import project_orbitals  # noqa: E402
from mom import mom_hook, seeded_run  # noqa: E402

from atomterm.cases import CAMPAIGN, ReferenceSpec  # noqa: E402
from atomterm.reference import (_forced_rohf_occ, occupation_postcheck,
                                run_rohf)
from atomterm.system import build_mol  # noqa: E402
from pyscf.scf import hf as hf  # noqa: E402

key = "Sc-unc"
spec = CAMPAIGN[key]
head = spec.series[0]
mol = build_mol(spec.atom_config())
s = mol.intor_symmetric("int1e_ovlp")
_, mf = run_rohf(spec.atom_config(), f"elements/{key}/rungs/rohf.chk")
dm = mf.make_rdm1()

print(f"stored mo_energy range: [{mf.mo_energy.min():.2f}, "
      f"{mf.mo_energy.max():.2f}]")
dm2 = np.asarray(dm)
if dm2.ndim == 3:                       # (alpha, beta) form
    dm2 = dm2[0] + dm2[1]
ev_dm = np.linalg.eigvalsh(dm2)
print(f"||dm|| eigen-range: [{float(ev_dm.min()):.3f}, "
      f"{float(ev_dm.max()):.3f}] "
      f"(trace {float(np.trace(dm2 @ s)):.1f})")

for tag, f in (("get_fock(dm)", mf.get_fock(dm=dm)),
               ("h+veff(dm)", mf.get_hcore() + mf.get_veff(mol, dm))):
    e, c = hf.eig(f, s)
    n_occ = int((np.asarray(mf.mo_occ) > 0).sum())
    print(f"{tag}: occupied eps [{e[:n_occ].min():.2f}, "
          f"{e[:n_occ].max():.2f}]; full range [{e.min():.2f}, "
          f"{e.max():.2f}]")

veff = mf.get_veff(mol, dm)
print(f"veff eigen-range: [{np.linalg.eigvalsh(veff).min():.3f}, "
      f"{np.linalg.eigvalsh(veff).max():.3f}]")

# ---- the MOM run + the post-check disagreement, on the record
occ_s = np.asarray(mf.mo_occ)
c_seed = mf.mo_coeff[:, occ_s > 0]
m = seeded_run(mol, c_seed, occ_s[occ_s > 0], max_cycle=40)
hook = _forced_rohf_occ(mol, list(head.closed), list(head.open))
occ_hook = hook(m.mo_energy, m.mo_coeff)
occ_run = np.asarray(m.mo_occ)
print(f"MOM rerun of the trusted seed: E = {m.e_tot:.8f} "
      f"conv={m.converged}")
print(f"  run occ  (nonzero idx): {np.flatnonzero(occ_run > 0)}")
print(f"  hook occ (nonzero idx): {np.flatnonzero(occ_hook > 0)}")
print(f"  hook's closed picks' eps: "
      f"{np.sort(m.mo_energy[occ_hook == 2.0])}")
norms = np.linalg.norm(m.mo_coeff, axis=0)
print(f"  hook's closed picks' ||c||_2: "
      f"{np.sort(norms[occ_hook == 2.0])[::-1]}")
