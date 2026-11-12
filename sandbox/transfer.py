"""STAGE 1b (ticket 37), part 2: the known-answer transfer (battery
anchor 5) + the P1 decision runs on the twins -- the SANDBOX form of
reference.projected_source and the dm0 seed (the src seam is Stage 1c).

ANCHOR 5 -- the known-answer transfer: the contracted heads' occupied
MOs projected into the unc bases, re-relaxed under the heads' own
forced-occupation hooks.  Certificate: convergence + landing a
DECLARED basin exactly (the trusted unc anchor, or the contracted
basin's unc image -- the map's 2026-11-07/08 record: the seeded Ti
transfer reproduced the contracted basin exactly; the deepest basin
came from the census).

P1 DECISION RUNS -- the disease legs: the series' non-head
configurations, seeded from the ON-DEMAND 6-31G source (forced SCF
there is tame), projected, re-relaxed on the unc basis.  Certificate:
convergence + occupation_postcheck + a sane energy vs the contracted
member's anchor.  The d3 leg on Sc-unc is the configuration whose
census diverged to +1.3e7 on every arm -- if THIS converges, P1 wins.

Sandbox protocol: in-repo, imports src, imported by nothing; the
seed mechanics replicate what Stage 1c formalizes (run_rohf's dm0
thread) without touching src.
"""
import pathlib
import sys

import numpy as np
from dataclasses import replace

sys.path.insert(0, "sandbox")
from project import project_orbitals  # noqa: E402

from atomterm.cases import CAMPAIGN, ReferenceSpec  # noqa: E402
from atomterm.reference import (_forced_rohf_occ, occupation_postcheck,
                                run_rohf)
from atomterm.system import build_mol  # noqa: E402
from pyscf import scf  # noqa: E402

SCRATCH = pathlib.Path("sandbox/scratch")


def seeded_rohf(mol, occ_closed, occ_open, dm0, accel="diis"):
    """One forced-occupation ROHF on `mol` seeded with dm0 -- the
    sandbox replica of the Stage-1c seam (run_rohf(dm0=))."""
    m = scf.ROHF(mol)
    m.get_occ = _forced_rohf_occ(mol, list(occ_closed), list(occ_open))
    m.conv_tol = 1e-12
    m.verbose = 0
    if accel == "newton":
        m = m.newton()
    m.run(dm0=dm0)
    return m


def projected_dm0(mol_src, mf_src, mol_tgt, thresh=1e-8):
    """The projected occupied density: the source's occupied MOs
    projected into the target basis (the sandbox form of
    reference.projected_source's product)."""
    occ = np.asarray(mf_src.mo_occ)
    cols = occ > 0
    c_b, rep = project_orbitals(mol_tgt, mol_src,
                                mf_src.mo_coeff[:, cols], thresh=thresh)
    dm0 = (c_b * occ[cols]) @ c_b.conj().T
    return dm0, rep


def anchor5_transfer(key, contracted):
    """The head->head transfer: contracted head -> unc basis -> the
    unc head's own forced SCF, seeded."""
    spec_u = CAMPAIGN[key]
    head_u = spec_u.series[0]
    mol_u = build_mol(spec_u.atom_config())
    mol_c, mf_c = run_rohf(CAMPAIGN[contracted].atom_config(),
                           f"elements/{contracted}/rungs/rohf.chk")
    dm0, rep = projected_dm0(mol_c, mf_c, mol_u)
    mf = seeded_rohf(mol_u, head_u.closed, head_u.open, dm0.real)
    basins = {"unc anchor": head_u.e_ref,
              "contracted head": CAMPAIGN[contracted].series[0].e_ref}
    hits = {n: abs(mf.e_tot - e) for n, e in basins.items()}
    best = min(hits, key=hits.get)
    post = occupation_postcheck(mol_u, mf, head_u.closed, head_u.open)
    ok = mf.converged and hits[best] < 1e-4 and post
    return ok, (f"{key}: E = {mf.e_tot:.8f} conv={mf.converged} "
                f"post={post} -> {best} (d = {hits[best]*1e3:.3f} mEh); "
                f"capture_min = {rep['capture'].min():.4f}")


def p1_decision(key, contracted, tag, source_basis="6-31g"):
    """One disease leg: the non-head configuration, seeded from the
    on-demand 6-31G source."""
    spec_u = CAMPAIGN[key]
    member = next(s for s in CAMPAIGN[contracted].series
                  if s.tag == tag)
    # the on-demand source: forced SCF on the small basis, fresh
    cfg_s = replace(CAMPAIGN[contracted].atom_config(),
                    basis=source_basis, spin=member.spin)
    chk = SCRATCH / f"{key.lower()}_{source_basis}_{tag}.chk"
    chk.parent.mkdir(parents=True, exist_ok=True)
    chk.unlink(missing_ok=True)
    mol_s, mf_s = run_rohf(cfg_s, str(chk),
                           ref=ReferenceSpec(occ_closed=member.closed,
                                             occ_open=member.open))
    assert occupation_postcheck(mol_s, mf_s, member.closed, member.open)
    # project + re-relax on the unc basis
    mol_u = build_mol(replace(spec_u.atom_config(), spin=member.spin))
    dm0, rep = projected_dm0(mol_s, mf_s, mol_u)
    for accel in ("diis", "newton"):
        mf = seeded_rohf(mol_u, member.closed, member.open,
                         dm0.real, accel=accel)
        post = occupation_postcheck(mol_u, mf, member.closed,
                                    member.open)
        if mf.converged and post:
            de = (mf.e_tot - member.e_ref) * 1e3
            sane = mf.e_tot < member.e_ref + 0.05
            return (mf.converged and post and sane), (
                f"{key} {tag} <- {source_basis}: E = {mf.e_tot:.8f} "
                f"[{accel}] conv={mf.converged} post={post}; "
                f"vs contracted anchor {de:+.1f} mEh; "
                f"capture_min = {rep['capture'].min():.4f}")
    return False, (f"{key} {tag} <- {source_basis}: NO arm converged "
                   f"occupation-safe (the P1 branch fails)")


if __name__ == "__main__":
    print("=== anchor 5: the known-answer transfer (heads) ===")
    a5 = []
    for key, contracted in (("Sc-unc", "Sc"), ("Ti-unc", "Ti")):
        ok, msg = anchor5_transfer(key, contracted)
        a5.append(ok)
        print(("PASS  " if ok else "FAIL  ") + msg)
    print("=== the P1 decision runs (the disease legs) ===")
    p1 = []
    for key, contracted, tag in (("Sc-unc", "Sc", "d3"),
                                 ("Sc-unc", "Sc", "d2s1"),
                                 ("Ti-unc", "Ti", "d4"),
                                 ("Ti-unc", "Ti", "d3s1")):
        ok, msg = p1_decision(key, contracted, tag)
        p1.append(ok)
        print(("PASS  " if ok else "FAIL  ") + msg)
    print(f"anchor 5: {'ALL PASS' if all(a5) else 'FAILURES'}; "
          f"P1: {'WINS' if all(p1) else ('PARTIAL' if any(p1) else 'FAILS')}")
