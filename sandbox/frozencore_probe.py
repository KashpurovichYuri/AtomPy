"""The frozencore optimizer probe (the user's variety extension, the
(b) hardening): hold the CORE block in the optimizer's own terms --
every rotation pair touching a core orbital blocked, window+virtuals
free -- and measure whether the swap basin still lands on the twin
leg where the free optimizer found it (sc-unc 3d4s free ao pol_avg).

Expected if the mechanism is right: no swap (the semi-core lives in
the now-immobile block), in-window convergence, the energy between
the free-in-window CI reading and the swap's.
"""
import sys

import numpy as np

sys.path.insert(0, "sandbox")
from project import project_orbitals  # noqa: F401  (the bench import)

from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.cells import CellKey, _key_core, _window_columns  # noqa: E402
from atomterm.casscf import constrained_casscf  # noqa: E402
from atomterm.reference import run_rohf  # noqa: E402
from atomterm.stages import stage_facts  # noqa: E402
from atomterm.system import build_mol, sector_na_nb  # noqa: E402
from atomterm.windows import (build_full_mo,  # noqa: E402
                              occupied_core_idx)


def mask_frozen_core(self, nmo, ncore, ncas, frozen):
    """Allowed: every rotation among (active+virtual); blocked: any
    pair touching the core -- the core is immobile, the rest free.
    The ``frozen`` mask's complement in scope (window frozen there,
    core relaxed here)."""
    m = np.zeros((nmo, nmo), dtype=bool)
    m[ncore:, ncore:] = True
    return m


spec = CAMPAIGN["Sc-unc"]
cfg = spec.atom_config()
mol = build_mol(cfg)
_, mf = run_rohf(cfg, "elements/Sc-unc/rungs/rohf.chk")
key = CellKey.from_address("sc-unc_3d4s_2-1_casscf_ao_-_pol_avg")
na, nb = key.sector
ncore, c_core, c_pur = _key_core(key, mol, mf, spec.shell_bands)
act, ncas, rep = _window_columns(key, mol, mf, spec.shell_bands,
                                 spec.shell_radials, ncol=ncore,
                                 prior_bands=spec.prior_bands)
core_idx = occupied_core_idx(mf, ncore)
mo_full = build_full_mo(mf, core_idx, act, c_core=c_core)

chk = "sandbox/scratch/sc_unc_3d4s_frozencore.chk"
res = constrained_casscf(mf, ncas, (na, nb), chk, mo_full=mo_full,
                         mask=mask_frozen_core, name="frozencore")
print(f"frozencore: E = {res.energies[0]:.8f}  cycles={res.cycles} "
      f"conv_dE={res.conv_dE}")

s = mol.intor_symmetric("int1e_ovlp")


def span_ov(c, cref):
    def proj(cols):
        return cols @ np.linalg.inv(cols.conj().T @ s @ cols) \
            @ cols.conj().T
    return float(np.trace(proj(c) @ s @ proj(cref) @ s).real
                 / c.shape[1])


core_ref = mf.mo_coeff[:, occupied_core_idx(mf, ncore)]
c_out = res.mo_coeff
print("core span overlap (out vs ref):",
      round(span_ov(c_out[:, :ncore], core_ref), 6),
      "(1 = immobile)")
print("active span overlap (out vs input window):",
      round(span_ov(c_out[:, ncore:ncore + ncas], act), 6))

facts = stage_facts(mol, mf, res, c_out[:, ncore:ncore + ncas],
                    ncore, na, nb, tag="frozencore")
print("swap_basin:", facts.swap_basin,
      "| core_span_overlap:", facts.core_span_overlap)
print("e_ground:", f"{facts.e_ground:.8f}",
      "| max_de_vs_pyscf:", facts.max_de_vs_pyscf)
print("ground term:", facts.ground_term,
      "| <L2>:", round(facts.ground_l2, 4),
      "| comm_norms:", facts.comm_norm_h_lz,
      facts.analysis.comm_norm_h_l2)
print()
print("reference readings: CI in-window -759.7342 (+3.2) | free CASSCF")
print("swap -759.7593 (-21.9) | frozen-window -759.7343 (+3.1)")
