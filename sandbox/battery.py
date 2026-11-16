"""The subspace-ROHF battery (ticket 37 Stage 3, anchor 7).

Legs:
  exact -- the whole-AO-space degenerate anchor: subspace == free ROHF
           (family + the full 43-AO valence = everything), per filling.
  ion   -- the closed-shell cross-check: the valence-empty solve on the
           family span vs confined_core_orbitals's own Roothaan loop
           (two independent constructions, one confined problem).
  conf  -- one confined spot leg per filling (the construction's first
           contact with the real windows; the full table is table.py).

Run: OMP_NUM_THREADS=16 PYTHONPATH=src python sandbox/battery.py Sc-unc
"""

import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from atomterm.basis import series_window
from atomterm.cases import CAMPAIGN
from atomterm.core import closed_shell_core_energy, confined_core_orbitals
from atomterm.system import build_mol
from atomterm.windows import lowdin, select_valence_ao

from subspace import ao_family, load_substrate, subspace_rohf

CORE_SHELLS = ("1s", "2s", "2p", "3s", "3p")


def twin_layout(key: str):
    """The twin's assembly inputs: the case, mols, the trusted head and
    member substrates (each with its declared filling patterns)."""
    cs = CAMPAIGN[key]
    cfg = cs.atom_config()
    mol = build_mol(cfg)
    head = load_substrate(f"elements/{key}/rungs/rohf.chk")
    members = [(spec, load_substrate(
        f"elements/{key}/rungs/series/{spec.tag}.chk"))
        for spec in cs.series[1:]]
    return cs, cfg, mol, head, members


def fillings(cs, head, members):
    """The declared (patterns, anchor, substrate) rows per filling."""
    rows = [(*_spec_fields(cs.ref, head, head_ref=True), "head")]
    for spec, mf in members:
        rows.append((*_spec_fields(spec, mf), spec.tag))
    return rows


def _spec_fields(spec, mf, head_ref: bool = False):
    if head_ref:
        return (spec.occ_closed, spec.occ_open, spec.occ_e_ref, mf)
    return (spec.closed, spec.open, spec.e_ref, mf)


def main(key: str) -> None:
    cs, cfg, mol, head, members = twin_layout(key)
    fam = ao_family(mol, CORE_SHELLS, cs.shell_bands)
    s = mol.intor_symmetric("int1e_ovlp")
    rows = fillings(cs, head, members)

    # -- the exactness legs (whole space == free ROHF) ------------------
    # the TRUE whole space: the family plus EVERY remaining AO -- the
    # declared valence band families miss 14 undeclared columns on the
    # twins, and the raw-AO gauge of the short window carries the
    # metric-amplified garbage directions (eps +30-class) whose
    # energy-ascending picks destroy the filling; in the Loewdin gauge
    # of the full space the amplification is normalized away
    print(f"== {key}: the whole-space anchor legs ==", flush=True)
    n = mol.nao_nr()
    fam_rows = np.flatnonzero(np.any(np.abs(fam) > 1e-12, axis=1))
    rest_rows = [i for i in range(n) if i not in set(fam_rows.tolist())]
    eye_n = np.eye(n)
    worst = 0.0
    for closed, open_, anchor, mf, tag in rows:
        val = eye_n[:, rest_rows]
        t0 = time.time()
        r = subspace_rohf(mf, fam, cs.ncore, val, closed, open_,
                          conv_accelerator="newton")
        d = (r["e_tot"] - anchor) * 1e3
        worst = max(worst, abs(d))
        print(f"  [exact] {tag}: E = {r['e_tot']:.9f}  d_anchor = "
              f"{d:+.6f} mEh  ({time.time() - t0:.0f}s)", flush=True)
    # the Newton's conv_tol_grad = 1e-4 delivers energies to ~1e-9 Eh
    # class on the seed legs (measured: 2e-6 mEh); 1e-3 mEh is the
    # honest pin
    assert worst < 1e-3, f"the exactness leg failed: {worst} mEh"

    # -- the closed-shell cross-check -----------------------------------
    # the twins' family span is the ill-conditioned class: the
    # confined-core Roothaan floors at ~5e-10 past its own 100-cycle
    # budget (measured) -- the cross-check then runs where the span
    # conditioning allows it (the carbon-class 2-column family); the
    # twins' refusal is the recorded declared class
    print(f"== {key}: the closed-shell cross-check ==", flush=True)
    ion_charge = {"Sc": 3, "Ti": 4}[cs.symbol]
    mol_ion = build_mol(replace(cfg, charge=ion_charge, spin=0))
    s_ion = mol_ion.intor_symmetric("int1e_ovlp")
    fam_ion_raw = ao_family(mol_ion, CORE_SHELLS, cs.shell_bands)
    fam_ion = lowdin(fam_ion_raw, s_ion)
    mf_ion = __import__("pyscf").scf.ROHF(mol_ion)
    try:
        c_conf, rep = confined_core_orbitals(mf_ion, fam_ion, cs.ncore)
    except ValueError as e:
        print(f"  [ion] confined_core_orbitals REFUSES on the twins' "
              f"family span (the declared class): {e}", flush=True)
        c_conf = rep = None
    if c_conf is not None:
        e_loop = closed_shell_core_energy(mf_ion, c_conf)
        t0 = time.time()
        r = subspace_rohf(mf_ion, fam_ion_raw, cs.ncore,
                          np.zeros((mol_ion.nao_nr(), 0)), ("l0", "l1"), (),
                          conv_accelerator="newton")
        d = (r["e_tot"] - e_loop) * 1e3
        print(f"  [ion] loop = {e_loop:.9f}  subspace = {r['e_tot']:.9f}  "
              f"d = {d:+.6f} mEh  ({time.time() - t0:.0f}s)  "
              f"(loop: {rep['cycles']} cycles, comm {rep['comm']:.1e})",
              flush=True)
        assert abs(d) < 1e-3, f"the ion cross-check failed: {d} mEh"

    # -- one confined spot leg per filling (the @r1 3d4s window) --------
    print(f"== {key}: the confined spot legs ==", flush=True)
    for closed, open_, anchor, mf, tag in rows:
        win, _rep_w = series_window(mol, [head], s, ("3d", "4s"),
                                    cs.prior_bands, {"3d": 1, "4s": 1},
                                    cs.ncore)
        t0 = time.time()
        r = subspace_rohf(mf, fam, cs.ncore, win, closed, open_,
                          conv_accelerator="newton")
        d = (r["e_tot"] - anchor) * 1e3
        print(f"  [conf] 3d4s@r1 {tag}: E = {r['e_tot']:.9f}  "
              f"dE = {d:+.3f} mEh  ({time.time() - t0:.0f}s)  "
              f"capture_min {min(r['report']['capture']):.3f}  "
              f"splits {r['report']['window_splits']}", flush=True)

    print(f"== {key}: battery legs all green ==", flush=True)


def carbon_cross_check() -> None:
    """The closed-shell cross-check on a well-conditioned span: the
    carbon fixture's [1s,2s] family, one confined pair, valence empty
    -- the shimmed kernel vs confined_core_orbitals's own Roothaan
    loop (two independent constructions, one confined problem)."""
    import tests.carbon_reference as CR
    from pyscf import scf as pscf_sc

    cfg = CAMPAIGN["C"].atom_config()
    # the C4+ ion (2e = the one confined pair): both constructions run
    # on it (the nuclear field is atom-blind to the charge; the
    # counting follows mol.nelectron)
    mol = build_mol(replace(cfg, charge=4, spin=0))
    n = mol.nao_nr()
    s = mol.intor_symmetric("int1e_ovlp")
    core1 = select_valence_ao(mol, ["1s"])
    core2 = select_valence_ao(mol, ["2s"])
    fam_rows = core1 + core2
    eye = np.eye(n)
    fam_raw = eye[:, fam_rows]
    fam = lowdin(fam_raw, s)
    mf = pscf_sc.ROHF(mol)   # the ion-like substrate (reads no orbitals)
    c_conf, rep = confined_core_orbitals(mf, fam, 1)
    e_loop = closed_shell_core_energy(mf, c_conf)
    t0 = time.time()
    r = subspace_rohf(mf, fam_raw, 1, np.zeros((n, 0)), ("l0", "l1"), (),
                      conv_accelerator="diis")
    d = (r["e_tot"] - e_loop) * 1e3
    print(f"[carbon ion] loop = {e_loop:.9f}  subspace = {r['e_tot']:.9f}  "
          f"d = {d:+.6f} mEh  ({time.time() - t0:.0f}s)  "
          f"(loop: {rep['cycles']} cycles, comm {rep['comm']:.1e})",
          flush=True)
    assert abs(d) < 1e-3, f"the carbon ion cross-check failed: {d} mEh"


if __name__ == "__main__":
    key = sys.argv[1] if len(sys.argv) > 1 else "Sc-unc"
    if key == "carbon":
        carbon_cross_check()
    else:
        main(key)
