"""The irrep-SCF probe for U14e -- the D2h irrep-occupation approach
(the PySCF ``irrep_nelec`` device, as in the mrpt example
``41-for_state_average.py``) tried on our uranium, compared against
the trusted C1 rung (the l-span character-forced ROHF).  Prints only.

STAGE 1  the irrep-identification table: which AO bands land in which
         D2h irrep (the "problem with irrep identification" made
         concrete: Ag hosts 6s+7s+two d's at once; each u-irrep hosts
         a 6p pair AND two f functions).
STAGE 2  the references: (i) the C1 rung (trusted); (ii) the D2h ROHF
         with ``irrep_nelec`` = {Ag 5, Au 1, B1u 3, B2u 3, B3u 2}
         (the closed shells by aufbau-within-block, the open
         f3/d1 m-class FORCED); (iii) the D2h ROHF with OUR
         l-span forced filling (the rung's vocabulary on the D2h
         framework).  Compared: energies, the l-span character
         postcheck, the open-shell m-structure, the valence-MO
         overlaps (gauge distances).
STAGE 3  the windows from the irrep-SCF (the user's "MOs directly?"
         question, measured): (i) the MO-direct valence slice (the
         13 columns dominated by f/d/s bands); (ii) the @r1 device
         distilled from the irrep reference (+ its own pair-MO s
         channel).  Frozen CASCI (5,1) on both: closure, <L2>,
         energy -- against the rung-based records (+140.5 mEh,
         <L2> 72.0000).

STAGE 4 (deferred): the optimized levels on the winning window.

Run from the repo root:
    OMP_NUM_THREADS=64 python scripts/probes/probe_u14e_irrep_scf.py \
        | tee logs/probe_u14e_irrep_scf.log
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "u14e"))

import numpy as np  # noqa: E402

SECTOR = (5, 1)
NCORE = 4
FILLING = {"s": 4, "p": 6, "d": 1, "f": 3}   # 6s2+7s2 / 6p6 / 6d1 / 5f3
DICT_CHOICE = {"Ag": 5, "B1g": 0, "B2g": 0, "B3g": 0,
               "Au": 1, "B1u": 3, "B2u": 3, "B3u": 2}   # sum 14


def _ao_irrep_table(mol):
    """Stage 1: per irrep, the AO band labels; per band, the irrep."""
    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    rows = {}
    for k, name in enumerate(mol.irrep_name):
        d = np.asarray(mol.symm_orb[k])
        w = np.abs(d).sum(axis=1)
        for i in np.flatnonzero(w > 1e-10):
            rows.setdefault(i, name)
    table = {}
    for i, lab in enumerate(labels):
        table.setdefault(lab, set()).add(rows.get(i, "?"))
    print("  [stage1] the band -> irrep table:", flush=True)
    by_irrep = {}
    for lab, irs in sorted(table.items()):
        for ir in irs:
            by_irrep.setdefault(ir, []).append(lab)
    for ir, labs in sorted(by_irrep.items()):
        print(f"    {ir:4s}: {', '.join(sorted(set(labs), key=len))}",
              flush=True)
    return rows, by_irrep


def _l_characters(mol, mf, tag):
    """The occupied MOs' l-sector S-weights -- the filling identity."""
    s = mf.get_ovlp()
    c = np.asarray(mf.mo_coeff)
    occ = np.asarray(mf.mo_occ) > 0
    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    import re as _re
    sect = {}
    for i, lab in enumerate(labels):
        m = _re.match(r"\d+([spdf])", lab)
        sect.setdefault(m.group(1), []).append(i)
    tot = {}
    n_e = 0.0
    for j in np.flatnonzero(occ):
        v = c[:, j]
        n_e += float(mf.mo_occ[j])
        for l, idx in sect.items():
            tot[l] = tot.get(l, 0.0) + mf.mo_occ[j] * float(
                v[idx] @ s[np.ix_(idx, idx)] @ v[idx])
    print(f"  [{tag}] the occupied l-characters: "
          + ", ".join(f"{l}: {w:.3f}" for l, w in sorted(tot.items()))
          + f"  (electrons {n_e:.1f}; intended {FILLING})", flush=True)
    return tot


def _open_structure(mol, mf, tag):
    """The open columns' dominant AO labels (the m-structure)."""
    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    for j in np.flatnonzero(np.asarray(mf.mo_occ) > 0):
        if mf.mo_occ[j] >= 2.0:
            continue
        v = np.abs(np.asarray(mf.mo_coeff)[:, j])
        top = np.argsort(-v)[:2]
        print(f"  [{tag}] open column {j}: "
              + ", ".join(f"{labels[i]} {v[i]:.3f}" for i in top),
              flush=True)


def _gauge_distance(mol, mf1, mf2, tag):
    """The occupied-set S-overlap singular values (the gauge distance)."""
    s = mf2.get_ovlp()
    o1 = np.flatnonzero(np.asarray(mf1.mo_occ) > 0)
    o2 = np.flatnonzero(np.asarray(mf2.mo_occ) > 0)
    if len(o1) != len(o2):
        print(f"  [{tag}] the occupied counts differ "
              f"({len(o1)} vs {len(o2)}) -- no gauge distance",
              flush=True)
        return
    x = np.asarray(mf1.mo_coeff)[:, o1].conj().T @ s \
        @ np.asarray(mf2.mo_coeff)[:, o2]
    sv = np.linalg.svd(x, compute_uv=False)
    print(f"  [{tag}] the occupied-set gauge distance: SVD "
          f"min {sv.min():.4f}, mean {sv.mean():.4f}, max "
          f"{sv.max():.4f}", flush=True)


def main() -> int:
    t00 = time.time()
    from pyscf import gto, scf
    from atomterm.cases import campaign_spec
    from atomterm.reference import (_forced_rohf_occ,  # noqa: E402
                                    occupation_postcheck)
    from vv.reference import load_reference

    spec = campaign_spec("U14e")
    mol_d2h = gto.M(atom=f"{spec.symbol} 0 0 0", basis=spec.basis,
                    spin=spec.spin, ecp=spec.ecp, symmetry="D2h",
                    verbose=0)
    print(f"  [mol] the D2h U atom: nao {mol_d2h.nao}, spin "
          f"{mol_d2h.spin}, the irreps {mol_d2h.irrep_name}", flush=True)

    # ---- stage 1: the identification table -----------------------------
    _ao_irrep_table(mol_d2h)

    # ---- stage 2: the three references ---------------------------------
    mol_c1, mf_c1 = load_reference()
    _l_characters(mol_c1, mf_c1, "rung C1")

    print("===== (ii) the D2h irrep-nelec ROHF =====", flush=True)
    mf_irr = scf.ROHF(mol_d2h)
    mf_irr.irrep_nelec = dict(DICT_CHOICE)
    mf_irr.conv_tol = 1e-12
    t0 = time.time()
    mf_irr.kernel()
    print(f"  [ii] E {mf_irr.e_tot:.8f}, converged {mf_irr.converged} "
          f"({time.time() - t0:.0f} s)", flush=True)
    # the verification: the occupations per irrep, honored or not
    from pyscf.symm import label_orb_symm
    try:
        ids = label_orb_symm(mol_d2h, mol_d2h.irrep_id,
                             mol_d2h.symm_orb, mf_irr.mo_coeff)
        sums = {mol_d2h.irrep_name[k]: 0.0
                for k in range(len(mol_d2h.irrep_name))}
        for j, k in enumerate(ids):
            sums[mol_d2h.irrep_name[k]] += float(mf_irr.mo_occ[j])
        print(f"  [ii] the occupations per irrep: {sums} vs the dict "
              f"{DICT_CHOICE}", flush=True)
    except ValueError as e:
        print(f"  [ii] the labeling refused: {e}", flush=True)
    _l_characters(mol_d2h, mf_irr, "ii irrep")
    _open_structure(mol_d2h, mf_irr, "ii irrep")

    print("===== (iii) the D2h l-span-forced ROHF =====", flush=True)
    mf_d2f = scf.ROHF(mol_d2h)
    mf_d2f.get_occ = _forced_rohf_occ(mol_d2h, ("l0", "l1"),
                                      ("l2", "l3"))
    mf_d2f.conv_tol = 1e-12
    t0 = time.time()
    mf_d2f.kernel()
    print(f"  [iii] E {mf_d2f.e_tot:.8f}, converged {mf_d2f.converged} "
          f"({time.time() - t0:.0f} s); the character postcheck "
          f"{occupation_postcheck(mol_d2h, mf_d2f, ('l0', 'l1'), ('l2', 'l3'))}",
          flush=True)
    _l_characters(mol_d2h, mf_d2f, "iii d2f")
    _open_structure(mol_d2h, mf_d2f, "iii d2f")

    _gauge_distance(mol_c1, mf_c1, mf_irr, "rung vs irrep")
    _gauge_distance(mol_c1, mf_c1, mf_d2f, "rung vs d2f")
    _gauge_distance(mol_d2h, mf_d2f, mf_irr, "d2f vs irrep")

    # ---- stage 3: the windows from the irrep reference -----------------
    from atomterm.windows import build_full_mo, occupied_core_idx
    from atomterm.operators import window_closure
    from vv.meters import meters
    from pyscf import mcscf

    key = None   # the windows are built directly; no key needed here

    def _frozen(tag, mo_win_source, mf_ref, c_core):
        mo_prod = build_full_mo(mf_ref, occupied_core_idx(mf_ref, NCORE),
                                mo_win_source, c_core=c_core)
        mo_win = mo_prod[:, NCORE:NCORE + 13]
        cl = window_closure(mol_d2h, mo_win)
        cas = mcscf.CASCI(mf_ref, 13, SECTOR)
        cas.canonicalization = False
        cas.kernel(mo_prod)
        expect = meters(mol_d2h, mo_win, *SECTOR)
        ci = np.asarray(cas.ci[0] if isinstance(cas.ci, list)
                        else cas.ci).ravel()
        l2 = expect("l2", ci)
        s2 = expect("s2", ci)
        de_meh = (float(cas.e_tot) - mf_c1.e_tot) * 1e3
        print(f"  [{tag}] CASCI e_tot {cas.e_tot:.8f}, dE {de_meh:+.1f} "
              f"mEh, <L2> {l2:.6f}, <S2> {s2:.6f}; the window's closure "
              f"||[P,Lz]|| {cl['lz']:.2e}, ||[P,Lperp]|| "
              f"{cl['lperp']:.2e}", flush=True)
        return l2, cl

    print("===== stage 3: the windows from the irrep reference =====",
          flush=True)
    # (i) MOs directly: the family's moc-kind selector (the columns of
    # the window's character, ranked by the named bands' S-projection)
    labels = [lab.split()[-1] for lab in mol_d2h.ao_labels()]
    from atomterm.windows import select_valence_mo
    perm = select_valence_mo(mf_irr, ["10f", "13d", "7s"], NCORE)
    mo_direct = np.asarray(mf_irr.mo_coeff)[:, perm][:, NCORE:NCORE + 13]
    print(f"  [(i) MO-direct] the selector's window columns: "
          f"{list(perm[NCORE:NCORE + 13])}", flush=True)
    _frozen("(i) MO-direct", mo_direct, mf_irr, None)

    # (ii) the @r1 device on the irrep reference
    from atomterm.basis import series_window
    act, rep = series_window(mol_d2h, [mf_irr], mf_irr.get_ovlp(),
                             ["5f", "6d", "7s"], spec.shell_bands,
                             {sh: 1 for sh in ("5f", "6d", "7s")}, NCORE)
    act = np.array(act, dtype=float)
    j_raw = 12   # the s radial: the last shell's single column
    # the irrep reference's own pair MO, s-pure sector (the r1 pattern)
    from atomterm.windows import occupied_core_idx as _oci
    core_idx = set(_oci(mf_irr, NCORE))
    occ_j = [j for j in range(mf_irr.mo_coeff.shape[1])
             if mf_irr.mo_occ[j] > 0 and j not in core_idx
             and labels[int(np.argmax(np.abs(
                 np.asarray(mf_irr.mo_coeff)[:, j])))].endswith("s")]
    print(f"  [(ii)] the irrep reference's s-dominated occupied "
          f"columns: {occ_j}", flush=True)
    if occ_j:
        v = np.asarray(mf_irr.mo_coeff)[:, occ_j[0]]
        s_rows = [i for i, lab in enumerate(labels)
                  if _re.match(r"\d+s$", lab)]
        vp = np.zeros_like(v)
        vp[s_rows] = v[s_rows]
        vp = vp / np.sqrt(float(vp @ mf_irr.get_ovlp() @ vp))
        act[:, j_raw] = vp
    _frozen("(ii) @r1-from-irrep", act, mf_irr, None)

    print(f"({time.time() - t00:.0f} s total)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
