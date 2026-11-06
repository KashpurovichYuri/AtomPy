"""Probe: the ROHF recheck (scratch; the log is the record, the probe
chk removed at the end).

1. A FRESH from-scratch gated ROHF on a probe chk path (the trusted
   cache is untouched; the basin gate retries init guesses until
   |E - occ_e_ref| <= 1e-4): the energy vs the trusted value is a
   basin comparison, never bitwise (fact 10).
2. The character audit: occupations, the open columns' band
   identities and m-contents, the one-particle closure certificate
   of the open-shell span, the core's m-scalarity.
3. The reference determinant's exact full-space <L^2> -- the spin
   block one-body norm identity,
       <D|L^2|D> = sum_a [ |c_a^al| + SS_a^al + |c_a^be| + SS_a^be
                           + 2 c_a^al c_a^be ],
   c_sigma = sum_{i in occ_sigma} L_aa,ii,
   SS_sigma = sum_{i in occ_sigma, k not in occ_sigma} |L_aa,ki|^2,
   machine-certified against the detbasis L^2 matrix on a small
   random-mix case BEFORE use (the identity construction, not the
   number, is what the certification pins).
"""
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[0] / "src"))
from pyscf import lib  # noqa: E402

from atomterm.angular import l_matrices  # noqa: E402
from atomterm.cases import campaign_spec  # noqa: E402
from atomterm.detbasis import (determinant_basis,  # noqa: E402
                               slater_condon_matrix)
from atomterm.operators import (build_operators, core_anisotropy,  # noqa: E402
                                l2_operator_parts, window_closure)
from atomterm.reference import run_rohf  # noqa: E402

TRUSTED = -50.38185233


def band_ao_groups(mol):
    """The AO row indices per band token (the label's leading
    ``<n><l>``)."""
    groups = {}
    for mu, (_, label, _sym) in enumerate(mol.ao_labels(fmt=False)):
        m = re.match(r"(\d+[spdfgh])", label)
        if m is None:
            raise ValueError(f"AO label {label!r} has no band token")
        groups.setdefault(m.group(1), []).append(mu)
    return groups


def det_l2_diagonal(l_mo, occ_a, occ_b):
    """The exact full-space <D|L^2|D> by the spin block norm identity
    (the module docstring).  ``l_mo`` = {'x','y','z'} MO-basis
    complex Hermitian; the occupations are spatial column index
    lists (closed columns appear in BOTH lists)."""
    total = 0.0
    for a in "xyz":
        l = l_mo[a]
        c_a = float(np.real(np.sum(l[np.ix_(occ_a, occ_a)].diagonal())))
        c_b = float(np.real(np.sum(l[np.ix_(occ_b, occ_b)].diagonal())))
        virt = [k for k in range(l.shape[0]) if k not in occ_a]
        ss_a = float(np.sum(np.abs(l[np.ix_(virt, occ_a)]) ** 2))
        virt_b = [k for k in range(l.shape[0]) if k not in occ_b]
        ss_b = float(np.sum(np.abs(l[np.ix_(virt_b, occ_b)]) ** 2))
        total += c_a ** 2 + ss_a + c_b ** 2 + ss_b + 2.0 * c_a * c_b
    return total


def certify_l2_identity(mol, mf):
    """Machine-certify the identity against the detbasis L^2 matrix:
    a random unitary mix of six canonical columns, sector (3,1), the
    closed-pair + open-alpha occupation pattern."""
    rng = np.random.default_rng(7)
    q = np.linalg.qr(rng.normal(size=(mf.mo_coeff.shape[0], 6)))[0]
    dets = determinant_basis(6, 3, 2)
    parts = l2_operator_parts(mol, q)
    mat = slater_condon_matrix(dets, one_body=parts[0], two_body=parts[1])
    det = (0, 1, 2, 6, 7)  # alphas {0,1,2}, betas {0,1} (norb=6)
    idx = dets.index(det)
    exact = float(np.real(mat[idx, idx]))
    l_mo = {a: q.conj().T @ l_matrices(mol)[a] @ q for a in "xyz"}
    formula = det_l2_diagonal(l_mo, [0, 1, 2], [0, 1])
    print(f"  certification: matrix {exact:.12f} formula {formula:.12f} "
          f"|diff| {abs(exact - formula):.2e}", flush=True)
    if abs(exact - formula) > 1e-9:
        raise RuntimeError("the L^2 norm identity failed certification")
    return abs(exact - formula)


def main():
    spec = campaign_spec("U14e")
    cfg = spec.atom_config()
    probe_chk = "elements/U14e/rungs/rohf_probe_fresh.chk"

    print("== ROHF recheck (fresh from-scratch SCF, basin-gated) ==",
          flush=True)
    mol, mf = run_rohf(cfg, probe_chk, ref=spec.ref)
    print(f"  fresh E = {mf.e_tot:.8f}  trusted = {TRUSTED:.8f}  "
          f"|dE| = {abs(mf.e_tot - TRUSTED):.2e}", flush=True)
    print(f"  converged = {mf.converged}", flush=True)

    occ = np.asarray(mf.mo_occ)
    closed = list(np.flatnonzero(occ == 2.0))
    opens = list(np.flatnonzero((occ > 0) & (occ < 2)))
    print(f"  occupations: {len(closed)} closed pairs (occ=2), opens "
          f"{opens} (occ {[occ[i] for i in opens]})", flush=True)

    certify_l2_identity(mol, mf)

    bands = band_ao_groups(mol)
    s = mf.get_ovlp()
    g = mf.mo_coeff.conj().T @ s @ mf.mo_coeff  # rows sum to 1 per column

    def band_profile(j, n=3):
        pops = {b: float(np.sum(g[np.ix_([mu for mu in rows], [j])]))
                for b, rows in bands.items()}
        top = sorted(pops.items(), key=lambda kv: -kv[1])[:n]
        return "  ".join(f"{b}:{w:.3f}" for b, w in top), pops

    print("-- the closed core (the fold) --", flush=True)
    for j in closed:
        prof, _ = band_profile(j)
        print(f"  mo[{j:3d}]  eps={mf.mo_energy[j]:9.4f}  {prof}",
              flush=True)
    print("-- the open shells --", flush=True)
    lz_mo = mf.mo_coeff.conj().T @ l_matrices(mol)["z"] @ mf.mo_coeff
    for j in opens:
        prof, _ = band_profile(j)
        m_lz = float(np.real(lz_mo[j, j]))
        print(f"  mo[{j:3d}]  eps={mf.mo_energy[j]:9.4f}  <Lz>=m "
              f"{m_lz:+.4f}   {prof}", flush=True)

    cert = window_closure(mol, mf.mo_coeff[:, opens])
    print(f"-- open-shell span closure certificate --", flush=True)
    print(f"  ||[P,Lz]||  = {cert['lz']:.4f}", flush=True)
    print(f"  ||[P,Lperp]|| = {cert['lperp']:.4f}", flush=True)
    aniso = core_anisotropy(mol, mf, len(closed))
    print(f"-- core m-scalarity (anisotropy of the closed-shell "
          f"density) -- {aniso:.3e}", flush=True)

    l_mo_full = {a: mf.mo_coeff.conj().T @ l_matrices(mol)[a]
                 @ mf.mo_coeff for a in "xyz"}
    l2det = det_l2_diagonal(l_mo_full, closed + opens, closed)
    print("-- the reference determinant (Ms=2, 9 spatial columns) --",
          flush=True)
    print(f"  <L^2> = {l2det:.6f}   (a pure 5S component reads 0; "
          f"the 5L term reads 72)", flush=True)
    print("  <S^2> = 6 exactly (ROHF construction: S = 2)", flush=True)

    Path(probe_chk).unlink(missing_ok=True)
    print("probe chk removed", flush=True)


if __name__ == "__main__":
    main()
