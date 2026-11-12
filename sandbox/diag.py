"""STAGE 1b diagnostics (ticket 37): WHY does a capture-1.0 projected
seed still diverge?  The map's 2026-11-08 record claims DIIS from the
trusted unc record rolls back to its own anchor -- if that holds but
the projected seed diverges, the difference lives in dm0; if both
diverge, the SCF step itself meets the disease regardless of seed.

Arms (Sc-unc head, d1s2):
  A  dm0 = the TRUSTED unc record's own density (same basis)
  B  dm0 = the projected contracted head (the transfer seed)
  C  arm B WITHOUT the forced-occupation hook (plain aufbau ROHF)
Each prints the energy trajectory's first/last few cycles.
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, "sandbox")
from project import project_orbitals  # noqa: E402
from transfer import projected_dm0, seeded_rohf  # noqa: E402

from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.reference import run_rohf  # noqa: E402
from atomterm.system import build_mol  # noqa: E402
from pyscf import scf  # noqa: E402

key, contracted = "Sc-unc", "Sc"
spec = CAMPAIGN[key]
head = spec.series[0]
mol_u = build_mol(spec.atom_config())

# the trusted unc record (cached chk -- trusted, no SCF)
_, mf_trusted = run_rohf(spec.atom_config(),
                         f"elements/{key}/rungs/rohf.chk")
print(f"trusted {key} head: E = {mf_trusted.e_tot:.8f} "
      f"(anchor {head.e_ref:.8f})")

# arm A: the trusted density as dm0 (S-project the mo columns = they
# are native, so dm = the record's own density)
occ = np.asarray(mf_trusted.mo_occ)
dm_a = (mf_trusted.mo_coeff[:, occ > 0] * occ[occ > 0]) \
    @ mf_trusted.mo_coeff[:, occ > 0].conj().T

# arm B: the projected contracted head
_, mf_c = run_rohf(CAMPAIGN[contracted].atom_config(),
                   f"elements/{contracted}/rungs/rohf.chk")
dm_b, rep = projected_dm0(mol_c := build_mol(CAMPAIGN[contracted]
                                              .atom_config()),
                          mf_c, mol_u)
print(f"arm B capture_min = {rep['capture'].min():.6f}")
print(f"||dm_a - dm_b||_F / ||dm_a||_F = "
      f"{np.linalg.norm(dm_a - dm_b) / np.linalg.norm(dm_a):.4f}")


def trace(dm0, hook=True, accel="diis", n=6):
    m = scf.ROHF(mol_u)
    if hook:
        from atomterm.reference import _forced_rohf_occ
        m.get_occ = _forced_rohf_occ(mol_u, list(head.closed),
                                     list(head.open))
    m.conv_tol = 1e-12
    m.max_cycle = 30
    m.verbose = 0
    es = []

    def cb(locals):
        es.append(locals["e_tot"] if "e_tot" in locals else np.nan)
    m.callback = cb
    if accel == "newton":
        m = m.newton()
    m.run(dm0=dm0.real.astype(np.float64))
    tail = [f"{e:.4f}" for e in es[:n]]
    end = [f"{e:.4f}" for e in es[-3:]]
    print(f"  conv={m.converged} E_final={m.e_tot:.6f} "
          f"cycles={len(es)} first={tail} last={end}")
    return m


print("arm A  (trusted native dm0, hook, diis):")
trace(dm_a)
print("arm B  (projected contracted dm0, hook, diis):")
trace(dm_b)
print("arm C  (projected dm0, NO hook, diis):")
trace(dm_b, hook=False)
print("arm D  (projected dm0, hook, newton):")
trace(dm_b, accel="newton")
