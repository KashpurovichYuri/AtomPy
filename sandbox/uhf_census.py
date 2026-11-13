"""The rich UHF convergence series on Ti-unc (the user's directive:
"just as for ROHF"), reshaped by the Stage-1b findings:

  the WEIGHT-hook arms (guess x accelerator) reproduce the disease --
  every forced-hook path fills the metric-amplified garbage on these
  bases, whatever the start (the old census's divergence, expected);
  the MOM-anchored arms (the per-spin overlap contest following the
  seed's columns -- the certified repair, now for UHF) from two
  seeds: the trusted ROHF chk's density (alpha=beta=ROHF) and the
  projected contracted-Ti head.

The basin table decides whether the UHF class can reach the ROHF's
deep basin (the cached anchor sits +6.5 mEh above it), and whether
MOM earns its UHF promotion.
"""
import sys

import numpy as np
from pyscf import lib, scf

sys.path.insert(0, "sandbox")
from transfer import projected_dm0  # noqa: E402

from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.reference import run_rohf  # noqa: E402
from atomterm.system import build_mol  # noqa: E402

spec = CAMPAIGN["Ti-unc"]
cfg = spec.atom_config()
mol = build_mol(cfg)
_, mf_rohf = run_rohf(cfg, "elements/Ti-unc/rungs/rohf.chk")
s = mol.intor_symmetric("int1e_ovlp")

# the seeds ---------------------------------------------------------
dm = np.asarray(mf_rohf.make_rdm1())
if dm.ndim == 3:
    dm = dm[0] + dm[1]
dm_rohf = np.asarray(dm.real, dtype=np.float64)
occ_rohf = np.asarray(mf_rohf.mo_occ)
occ_alpha_seed = np.where(occ_rohf == 2.0, 1.0, occ_rohf)   # 12 ones
occ_beta_seed = np.where(occ_rohf == 2.0, 1.0, 0.0)          # 10 ones

_, mf_c = run_rohf(CAMPAIGN["Ti"].atom_config(),
                   "elements/Ti/rungs/rohf.chk")
mol_c = build_mol(CAMPAIGN["Ti"].atom_config())
dm_proj, rep = projected_dm0(mol_c, mf_c, mol)
dm_proj = np.asarray(dm_proj.real, dtype=np.float64)


def mom_uhf_hook(mol, mo_seed, occ_alpha, occ_beta):
    """The MOM-anchored per-spin get_occ: occupations FOLLOW the
    seed's occupied columns by S-overlap (the certified ROHF form's
    per-spin twin; the contest state resets every call)."""
    seed = np.asarray(mo_seed)

    def _occ(mo_energy, mo_coeff=None):
        if mo_coeff is None:
            raise ValueError("MOM occupations need the MOs")
        occ = np.zeros((2, mo_coeff[0].shape[1]))
        for spin, seed_occ in ((0, occ_alpha), (1, occ_beta)):
            taken = np.zeros(mo_coeff[0].shape[1], dtype=bool)
            ov = np.abs(mo_coeff[spin].conj().T @ s @ seed)
            for j in np.flatnonzero(seed_occ > 0):
                free = np.flatnonzero(~taken)
                bi = free[np.argmax(ov[free, j])]
                occ[spin][bi] = seed_occ[j]
                taken[bi] = True
        return occ

    return _occ


def honest_occ(m, occ_closed, occ_open):
    """The character-of-occupied check (the honest spans)."""
    from atomterm.windows import select_span, span_projector
    p_c = span_projector(s, sum((select_span(mol, p)
                                 for p in occ_closed), []))
    p_o = span_projector(s, sum((select_span(mol, p)
                                 for p in occ_open), []))
    ok = True
    for spin in (0, 1):
        for i in np.flatnonzero(m.mo_occ[spin] > 0):
            w_c = float(np.real(m.mo_coeff[spin][:, i].conj()
                                @ p_c @ m.mo_coeff[spin][:, i]))
            w_o = float(np.real(m.mo_coeff[spin][:, i].conj()
                                @ p_o @ m.mo_coeff[spin][:, i]))
            if w_c < 0.5 and w_o < 0.5:
                ok = False
    return ok


def arm(label, hook=None, guess=None, accel="diis", dm0=None,
        max_cycle=120):
    m = scf.UHF(mol)
    if hook is not None:
        m.get_occ = hook
    m.conv_tol = 1e-10
    m.max_cycle = max_cycle
    m.verbose = 0
    if guess is not None:
        m.init_guess = guess
    if accel == "newton":
        m = m.newton()
    try:
        m.run(dm0=dm0)
    except Exception as e:
        print(f"  {label:36s} RAISED {type(e).__name__}: {str(e)[:36]}")
        return
    ss = m.spin_square()
    hc = honest_occ(m, ("l0", "l1"), ("l2",))
    print(f"  {label:36s} conv={str(m.converged):5s} "
          f"E={m.e_tot:.6f}  <S2>={ss[0]:.4f}  occ_ok={hc}")


print("reference readings: ROHF deep basin -848.407443 (the trusted")
print("cache) | the cached UHF anchor -848.400969 (the +6.5 mEh gap)")
print(f"seeds: ||d(dm_rohf, dm_proj)||/||dm_rohf|| = "
      f"{np.linalg.norm(dm_rohf - dm_proj) / np.linalg.norm(dm_rohf):.3f}")
print()
print("=== the weight-hook arms (the disease's reproduction) ===")
arm("weight hook, default guess diis", guess=None, accel="diis")
arm("weight hook, atom guess newton", guess="atom", accel="newton")
print("=== the MOM-anchored arms (the repair class) ===")
mom_rohf = mom_uhf_hook(mol, mf_rohf.mo_coeff, occ_alpha_seed,
                        occ_beta_seed)
arm("MOM, ROHF-seed dm0, diis", hook=mom_rohf, dm0=dm_rohf)
arm("MOM, ROHF-seed dm0, newton", hook=mom_rohf, accel="newton",
    dm0=dm_rohf)
mom_proj = mom_uhf_hook(mol, mf_rohf.mo_coeff, occ_alpha_seed,
                        occ_beta_seed)
arm("MOM, projected-seed dm0, diis", hook=mom_proj, dm0=dm_proj)
arm("MOM, projected-seed dm0, newton", hook=mom_proj, accel="newton",
    dm0=dm_proj)
