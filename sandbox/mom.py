"""STAGE 1b (ticket 37), the mechanism + the repair prototype.

MEASURED (diag.py): the fresh-SCF disease on the twins is NOT the
guess and NOT the conditioning -- the CHARACTER HOOK destabilizes
even the trusted record's own density (+1.4e7 from cycle one), while
the same seed without the hook lands the anchor.

HYPOTHESIS (this file's arm 0): the hook's projector weight c'Pc is
S-bounded and BLIND to numerical garbage -- a near-null AO
combination (Euclidean norm ~ 1/sqrt(lambda_min) ~ 45 on these
bases, lambda_min ~ 5e-4 per the Stage-1a spectra) is S-normalized
and carries weight ~ 1 on the very span being forced.  The hook
assigns the electron to garbage; the density is S-bounded but
Euclidean-huge; the two-body energy explodes.

THE REPAIR (arm A'/B'/P1): the MOM-anchored occupation hook --
occupations FOLLOW the seed's occupied columns (max S-overlap), the
seed being the configuration declaration (the projected 6-31G /
contracted source).  run_folded's certified pick discipline, promoted
to the fresh runner.  The character spans stay as the POST-CHECK.
"""
import pathlib
import sys

import numpy as np
from dataclasses import replace

sys.path.insert(0, "sandbox")
from project import project_orbitals  # noqa: E402
from transfer import p1_decision, projected_dm0  # noqa: E402

from atomterm.cases import CAMPAIGN, ReferenceSpec  # noqa: E402
from atomterm.reference import (_forced_rohf_occ, occupation_postcheck,
                                run_rohf)
from atomterm.system import build_mol  # noqa: E402
from pyscf import scf  # noqa: E402
from pyscf.scf import hf  # noqa: E402


def honest_post(mol, mf, occ_closed, occ_open):
    """The CHARACTER-OF-OCCUPIED post-check (the well-posed form on
    ill-conditioned bases): verify the CONVERGED OCCUPIED columns
    carry the declared character -- never re-rank all columns (the
    virtual tail's metric-amplified garbage carries w ~ 1 too, the
    measured disease)."""
    from atomterm.windows import select_span, span_projector
    s = mol.intor_symmetric("int1e_ovlp")
    p_closed = span_projector(
        s, sum((select_span(mol, p) for p in occ_closed), []))
    p_open = span_projector(
        s, sum((select_span(mol, p) for p in occ_open), []))
    occ = np.asarray(mf.mo_occ)
    n_closed = int(((mol.nelectron - mol.spin) // 2))
    n_open = int(mol.spin)
    closed = np.flatnonzero(occ == 2.0)
    opens = np.flatnonzero(occ == 1.0)
    if len(closed) != n_closed or len(opens) != n_open:
        return False
    for i in closed:
        w = float(np.real(mf.mo_coeff[:, i].conj() @ p_closed
                          @ mf.mo_coeff[:, i]))
        if w < 0.5:
            return False
    for i in opens:
        w = float(np.real(mf.mo_coeff[:, i].conj() @ p_open
                          @ mf.mo_coeff[:, i]))
        if w < 0.5:
            return False
    return True


def arm0_mechanism(key="Sc-unc"):
    """One cycle, instrumented: F[dm_trusted] -> MOs -> the hook's
    ranking vs the columns' Euclidean norms."""
    spec = CAMPAIGN[key]
    head = spec.series[0]
    mol = build_mol(spec.atom_config())
    _, mf = run_rohf(spec.atom_config(), f"elements/{key}/rungs/rohf.chk")
    s = mol.intor_symmetric("int1e_ovlp")
    dm = mf.make_rdm1()
    f = mf.get_fock(dm=dm)
    from pyscf.scf import hf as _hf
    e, c = _hf.eig(f, s)
    # the honest projector weights, the hook's own algebra
    from atomterm.windows import select_span, span_projector
    p_open = span_projector(
        s, sum((select_span(mol, p) for p in head.open), []))
    p_closed = span_projector(
        s, sum((select_span(mol, p) for p in head.closed), []))
    w_open = np.einsum("ai,ab,bi->i", c.conj(), p_open, c).real
    w_closed = np.einsum("ai,ab,bi->i", c.conj(), p_closed, c).real
    norms = np.linalg.norm(c, axis=0)
    n_closed = (mol.nelectron - mol.spin) // 2
    picked = np.argsort(w_closed)[::-1][:n_closed]
    print(f"arm 0  the CLOSED ranking (what the hook fills first):")
    for i in picked:
        print(f"    w_closed = {w_closed[i]:.6f}  ||c||_2 = "
              f"{norms[i]:8.2f}  eps = {e[i]:12.2f}")
    bad = [i for i in picked if norms[i] > 10]
    print(f"    -> {len(bad)} of {n_closed} closed picks are "
          f"Euclidean-garbage (norm > 10): the tie at w ~ 1 lands the "
          f"core on near-null s/p combos" if bad else
          f"    -> closed picks all physical this draw")
    print(f"    open-span weight of the norm>10 columns: "
          f"max = {w_open[norms > 10].max():.4f} (the garbage is "
          f"s/p-family -- the CLOSED spans)")
    return len(bad) > 0


def mom_hook(mol, c_seed, occ_seed):
    """The MOM-anchored get_occ: occupations follow the seed's
    occupied columns by max S-overlap (garbage never wins the overlap
    contest -- near-null combos are S-orthogonal to physical
    columns).  The contest state RESETS each call: MO indices
    re-diagonalize every cycle, so a cross-cycle 'taken' set would
    assign cycle 2 to whatever is left (the bug the first prototype
    caught)."""
    s = mol.intor_symmetric("int1e_ovlp")
    order = np.argsort(-np.asarray(occ_seed))     # closed pairs first

    def _occ(mo_energy, mo_coeff=None):
        if mo_coeff is None:
            raise ValueError("MOM occupations need the MOs")
        occ = np.zeros(len(mo_energy))
        ov = np.abs(mo_coeff.conj().T @ s @ c_seed)
        taken = np.zeros(len(mo_energy), dtype=bool)
        for j in order[np.asarray(occ_seed) > 0]:
            free = np.flatnonzero(~taken)
            bi = free[np.argmax(ov[free, j])]
            occ[bi] = occ_seed[j]
            taken[bi] = True
        return occ

    return _occ


def seeded_run(mol, c_seed, occ_seed, accel="diis", max_cycle=60):
    m = scf.ROHF(mol)
    m.get_occ = mom_hook(mol, c_seed, occ_seed)
    m.conv_tol = 1e-12
    m.max_cycle = max_cycle
    m.verbose = 0
    if accel == "newton":
        m = m.newton()
    dm0 = (c_seed * occ_seed) @ c_seed.conj().T
    m.run(dm0=dm0.real.astype(np.float64))
    return m


def armAB(key, contracted):
    spec = CAMPAIGN[key]
    head = spec.series[0]
    mol_u = build_mol(spec.atom_config())
    _, mf_c = run_rohf(CAMPAIGN[contracted].atom_config(),
                       f"elements/{contracted}/rungs/rohf.chk")
    occ = np.asarray(mf_c.mo_occ)
    c_b, rep = project_orbitals(mol_u,
                                build_mol(CAMPAIGN[contracted]
                                          .atom_config()),
                                mf_c.mo_coeff[:, occ > 0], thresh=1e-8)
    mf = seeded_run(mol_u, c_b, occ[occ > 0])
    post = honest_post(mol_u, mf, head.closed, head.open)
    hit = abs(mf.e_tot - head.e_ref)
    return mf.converged and post and hit < 1e-3, (
        f"{key} MOM-seeded: E = {mf.e_tot:.8f} conv={mf.converged} "
        f"post={post} | d_anchor = {hit*1e3:.3f} mEh")


def p1_mom(key, contracted, tag, source_basis="6-31g"):
    spec_u = CAMPAIGN[key]
    member = next(s for s in CAMPAIGN[contracted].series if s.tag == tag)
    cfg_s = replace(CAMPAIGN[contracted].atom_config(),
                    basis=source_basis, spin=member.spin)
    chk = pathlib.Path(f"sandbox/scratch/{key.lower()}_"
                       f"{source_basis}_{tag}.chk")
    chk.unlink(missing_ok=True)
    mol_s, mf_s = run_rohf(cfg_s, str(chk),
                           ref=ReferenceSpec(occ_closed=member.closed,
                                             occ_open=member.open))
    occ = np.asarray(mf_s.mo_occ)
    mol_u = build_mol(replace(spec_u.atom_config(), spin=member.spin))
    c_b, rep = project_orbitals(mol_u, mol_s, mf_s.mo_coeff[:, occ > 0],
                                thresh=1e-8)
    best = None
    for accel in ("diis", "newton"):
        mf = seeded_run(mol_u, c_b, occ[occ > 0], accel=accel)
        post = honest_post(mol_u, mf, member.closed, member.open)
        de = (mf.e_tot - member.e_ref) * 1e3
        if best is None or (mf.converged and post
                            and not best[0]):
            best = (mf.converged and post, accel, mf.e_tot, de, post)
        if mf.converged and post:
            return True, (f"{key} {tag} <- {source_basis} MOM: "
                          f"E = {mf.e_tot:.8f} [{accel}] "
                          f"conv={mf.converged} post={post}; vs "
                          f"contracted anchor {de:+.1f} mEh; "
                          f"capture_min={rep['capture'].min():.4f}")
        print(f"    [{accel}] E = {mf.e_tot:.6f} conv={mf.converged} "
              f"post={post} (d_contracted {de:+.1f} mEh)")
    return False, (f"{key} {tag}: MOM-seeded fails "
                   f"(best: {best[1]} E={best[2]:.6f})")


if __name__ == "__main__":
    print("=== arm 0: the mechanism (one instrumented cycle) ===")
    arm0_mechanism()
    print("=== arms A'/B': the MOM-seeded head transfer ===")
    r = []
    for key, contracted in (("Sc-unc", "Sc"), ("Ti-unc", "Ti")):
        ok, msg = armAB(key, contracted)
        r.append(ok)
        print(("PASS  " if ok else "FAIL  ") + msg)
    print("=== P1' the decision runs (6-31G seeds, MOM) ===")
    for key, contracted, tag in (("Sc-unc", "Sc", "d3"),
                                 ("Sc-unc", "Sc", "d2s1"),
                                 ("Ti-unc", "Ti", "d4"),
                                 ("Ti-unc", "Ti", "d3s1")):
        ok, msg = p1_mom(key, contracted, tag)
        r.append(ok)
        print(("PASS  " if ok else "FAIL  ") + msg)
    print("VERDICT:", "P1' (MOM-seeded) WINS" if all(r)
          else ("PARTIAL" if any(r) else "FAILS"))
