"""Per-rung ground-root <Lz^2> re-analysis from the cached chkfiles.

No kernel re-runs: the ROHF / frac-occ references load from their cached
chkfiles (trusted), every rung's window loads from ``mcscf/mo_coeff``,
the fold is rebuilt exactly as ``treatments._rung_payload`` builds it
(same mf, same spherical_core flag, core columns from the rung's own
mo_coeff), and the ground eigenvector of the folded H gives
<v|Lz@Lz|v>.  The independent v|L^2|v is printed against the stored
ground_l2 as the faithfulness check (identical up to degenerate-cluster
gauge).

Usage: OPENBLAS_NUM_THREADS=64 python scripts/lz2_analysis.py [Ti]
"""

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
from atomterm.reference import ReferenceSpec, run_frac_occ_scf, run_rohf
from atomterm.system import AtomConfig, ncore_of, sector_occupations
from atomterm.windows import select_valence_ao, build_full_mo, occupied_core_idx
from atomterm.core import purified_core_orbitals
from atomterm.operators import build_operators
from atomterm.treatments import _registry

# (cfg, shell patterns, nelec_act, sa_weights, frac_nelec_open, frac_n_open, ref)
CASES = {
    "Ti": (
        AtomConfig("Ti", "ccpvdz", spin=2),
        ["3d", "4s"], 4, [1 / 7] * 7, 4, 6,
        ReferenceSpec(occ_closed=("1s", "2s", "3s", "4s", "2p", "3p"),
                      occ_open=("3d",), occ_e_ref=-848.40648875),
    ),
}


def _load_mo(tag):
    """Window-substrate MO matrix: mcscf/mo_coeff for optimized rungs,
    rebuilt from the reference for the CASCI rungs (their chkfiles store
    only energies+ci; the substrate is deterministic)."""
    from pyscf.lib import chkfile
    mo = chkfile.load(chk(tag), "mcscf/mo_coeff")
    if mo is not None:
        return mo
    if tag in ("casci_raw", "casci_avg"):
        return build_full_mo(mf, occupied_core_idx(mf, ncore), act_idx)
    if tag == "casci_frac":
        return build_full_mo(mf_f, occupied_core_idx(mf_f, ncore), act_idx,
                             c_core=mf_f.mo_coeff[:, :ncore])
    if tag == "casci_avg_pur":
        return build_full_mo(mf, occupied_core_idx(mf, ncore), act_idx,
                             c_core=c_pur)
    return None


only = sys.argv[1:] or ["Ti"]
for name in only:
    cfg, patterns, nelec, sa_w, frac_nelec, frac_n, ref = CASES[name]
    out = Path("elements") / name / "ladder"
    label, spin = cfg.label, cfg.spin
    chk = lambda tag: str(out / f"{label}_{tag}.chk")  # noqa: E731

    mol, mf = run_rohf(cfg, chk("rohf"), ref=ref)
    act_idx = select_valence_ao(mol, patterns)
    mol_f, mf_f = run_frac_occ_scf(cfg, chk("fracocc"), frac_nelec, frac_n,
                                   open_idx=act_idx)
    ncas = len(act_idx)
    na, nb = sector_occupations(nelec, spin)
    ncore = ncore_of(mol.nelectron, nelec)
    ndet = math.comb(ncas, na) * math.comb(ncas, nb)
    c_pur, _ = purified_core_orbitals(mol, mf, ncore)
    registry = _registry(mf, mf_f, ncore, ncas, (na, nb), act_idx,
                         ndet, sa_w, chk, c_pur=c_pur)

    payload = json.load(open(out.parent / "ladder.json"))
    print(f"\n===== {name}: ground-root <Lz^2> in the rung's own fold =====")
    print(f"{'rung':<20} {'<Lz^2>':>12} {'v|L2|v':>12} "
          f"{'stored L2':>12} {'diff':>9}")
    for rung in registry:
        tag = rung.tag
        mo = _load_mo(tag)
        if mo is None:
            print(f"{tag:<20}  (no chk)")
            continue
        c_act = mo[:, ncore:ncore + ncas]
        ops = build_operators(
            mol, rung.fold_mf(), c_act, ncore, na, nb,
            spherical_core=rung.spherical_core, c_core=mo[:, :ncore])
        _, v = np.linalg.eigh(ops["H"])
        g = v[:, 0]
        lz2 = float(np.real(g @ (ops["Lz"] @ (ops["Lz"] @ g))))
        l2 = float(np.real(g @ (ops["L2"] @ g)))
        stored = payload[tag]["ground_l2"]
        print(f"{tag:<20} {lz2:>12.6f} {l2:>12.6f} {stored:>12.6f} "
              f"{abs(l2 - stored):>9.1e}")
