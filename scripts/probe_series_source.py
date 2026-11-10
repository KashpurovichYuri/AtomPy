"""The series-source probe: the config-ladder averaged distillation vs
the 1:1 authors' contractions -- the universal-construction evidence.

For each atom (the contracted transition metals, where the 1:1
identity realization is the incumbent): run the low-configuration
forced-occupation ROHF ladder, average the VALENCE density's
member-summed shell blocks over the ladder, and solve the s_r-metric
eigenproblem -- the occupation-averaged radial prior, with its
spectrum read directly (single configs are rank-1; the average
carries the config-spread as a second radial).

The priors compared on the SAME CASCI path (build_operators eigh,
anchored to the stored 1:1 totals through the prior-independent core
constant):

- ``1:1``    the authors' bare contractions (the incumbent records);
- ``neutral`` the reference configuration's own block top (rank-1);
- ``series1`` the ladder average's mean radial (k=1);
- ``series2`` mean + config-spread (k=2, where the spectrum certifies);
- ``occ1``   the 4s occupied-only variant (the empty-s configs dropped
  from the s average -- the one-probe question).

Usage:
    python scripts/probe_series_source.py Sc > elements/Sc/logs/series_source_probe.log
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pyscf import scf  # noqa: E402

from atomterm.operators import build_operators, window_closure  # noqa: E402
from atomterm.reference import ReferenceSpec, run_rohf  # noqa: E402
from atomterm.system import AtomConfig  # noqa: E402
from atomterm.windows import (build_full_mo, occupied_core_idx,  # noqa: E402
                              shell_slots, span_projector)

# the declaration per atom: the config ladder (closed shells, open
# shells, spin), the valence realizations (the measured greedy lists),
# the comparison cells (address, sector), the anchor column
PLAN = {
    "Sc": {
        "basis": "ccpvdz",
        "ar": ("1s", "2s", "2p", "3s", "3p"),
        "ladder": (("d1s2", ("1s", "2s", "2p", "3s", "3p", "4s"), ("3d",), 1),
                   ("d2s1", ("1s", "2s", "2p", "3s", "3p"), ("3d", "4s"), 3),
                   ("d3",   ("1s", "2s", "2p", "3s", "3p"), ("3d",), 3)),
        "bands": {"3d": ("3d", "4d"), "4s": ("4s", "6s")},
        "cells": (("sc_3d_2-1_casci_ao_-_pol_avg", ("3d",), (2, 1)),
                  ("sc_3d4s_2-1_casci_ao_-_pol_avg", ("3d", "4s"), (2, 1))),
    },
    "Ti": {
        "basis": "ccpvdz",
        "ar": ("1s", "2s", "2p", "3s", "3p"),
        "ladder": (("d2s2", ("1s", "2s", "2p", "3s", "3p", "4s"), ("3d",), 2),
                   ("d3s1", ("1s", "2s", "2p", "3s", "3p"), ("3d", "4s"), 4),
                   ("d4",   ("1s", "2s", "2p", "3s", "3p"), ("3d",), 4)),
        "bands": {"3d": ("3d", "4d"), "4s": ("4s", "6s")},
        "cells": (("ti_3d_2-0_casci_ao_-_pol_avg", ("3d",), (2, 0)),
                  ("ti_3d4s_3-1_casci_ao_-_pol_avg", ("3d", "4s"), (3, 1))),
    },
}
NCOL = 9                                   # the fold boundary (valence)


def shell_geometry(mol, labels, shell, bands):
    """The member chains and per-band per-member widths (the shell
    walker -- never label lists)."""
    ell = "spdfgh".index(shell[-1])
    n_members = 2 * ell + 1
    bm = {b: {} for b in bands}
    for _ib, sl, _c, rows in shell_slots(mol):
        if sl != ell:
            continue
        nm = labels[rows[0]]
        nm = nm[: nm.index(shell[-1]) + 1]
        if nm in bm:
            for m, r in enumerate(rows):
                bm[nm].setdefault(m, []).append(int(r))
    absent = [b for b in bands if not bm[b]]
    if absent:
        raise ValueError(f"bands {absent} absent")
    chain = [[r for b in bands for r in bm[b][m]]
             for m in range(n_members)]
    widths = [len(bm[b][0]) for b in bands]
    return chain, widths


def valence_block(mf, chain):
    """The member-summed shell block of the VALENCE density (the core
    columns never enter -- the measured core-dominance rejection)."""
    cv = mf.mo_coeff[:, NCOL:]
    wv = np.asarray(mf.mo_occ)[NCOL:]
    d = (cv * wv) @ cv.conj().T
    b = sum(d[np.ix_(chain[m], chain[m])] for m in range(len(chain)))
    return 0.5 * (b + b.T).real


def metric_radials(block, s_r, k):
    """The s_r-metric eigenproblem's top-k radials: eigh of the
    transported operator, vectors back to s_r-orthonormal chain
    space, and the normalized spectrum (the k certificate)."""
    l = np.linalg.cholesky(s_r)
    li = np.linalg.inv(l)
    t = li @ block @ li.T
    t = 0.5 * (t + t.T)
    lam, q = np.linalg.eigh(t)
    lam, q = lam[::-1].real, q[:, ::-1].real
    v = li.T @ q[:, :k]                    # s_r-orthonormal columns
    err = np.abs(v.T @ s_r @ v - np.eye(k)).max()
    assert err < 1e-9, f"s_r-orthonormality {err:.2e}"
    return v, lam


def assemble(mol_nao, radial, chain):
    """Radial columns x the full m-set, radial-major (k=1 reads like a
    label-selected m-set)."""
    n_members = len(chain)
    kv = radial.shape[1]
    blk = np.zeros((mol_nao, kv * n_members))
    col = 0
    for j in range(kv):
        for m in range(n_members):
            blk[chain[m], col] = radial[:, j]
            col += 1
    return blk


def authors_vector(s_r, widths, bands, band):
    """The bare contraction (the incumbent 1:1 prior)."""
    v = np.zeros(len(s_r))
    i = 0
    for b, w in zip(bands, widths):
        if b == band:
            v[i:i + w] = 1.0
        i += w
    return (v / np.sqrt(v @ s_r @ v))[:, None]


def capture(c_win, mf, s, chain, thresh=0.05):
    """Each significant reference column's S-projection on the window
    span (the honest quality readout)."""
    rows = sorted(r for m in range(len(chain)) for r in chain[m])
    p = span_projector(s, rows)
    c = mf.mo_coeff[:, np.asarray(mf.mo_occ) > 0]
    norm = np.einsum("ai,ab,bi->i", c, s, c).real
    w = np.einsum("ai,ab,bi->i", c, p, c).real / norm
    sig = np.flatnonzero(w >= thresh)
    cr = c_win.T @ (s @ c[:, sig])
    cap = np.einsum("ij,ij->j", cr, cr).real / norm[sig]
    return cap


def main(atom):
    plan = PLAN[atom]
    cfg = AtomConfig(atom, plan["basis"])
    mol, rec = scf.chkfile.load_scf(f"elements/{atom}/rungs/rohf.chk")
    mf = scf.RHF(mol)
    mf.__dict__.update(rec)
    s = mf.get_ovlp()
    labels = [line.split()[-1] for line in mol.ao_labels()]
    nao = mol.nao_nr()
    print(f"# {atom} {plan['basis']} -- the series-source probe "
          f"(reference E = {mf.e_tot:.8f})")

    # the configuration ladder (forced occupations; the basins of the
    # degenerate fillings wander -- the energies recorded per run)
    series = []
    scratch = Path("probe_tmp")
    scratch.mkdir(exist_ok=True)
    for tag, closed, opens, spin in plan["ladder"]:
        cfg_x = AtomConfig(atom, plan["basis"], spin=spin)
        ref = ReferenceSpec(occ_closed=closed, occ_open=opens,
                            occ_e_ref=None)
        _m, mf_x = run_rohf(cfg_x, str(scratch / f"{atom}-{tag}.chk"),
                            ref=ref)
        series.append((tag, mf_x))
        print(f"  config {tag}: E = {mf_x.e_tot:.8f}")

    # per-shell geometry + the priors
    geo = {sh: shell_geometry(mol, labels, sh, bands)
           for sh, bands in plan["bands"].items()}
    neutral_tag = plan["ladder"][0][0]
    priors = {}                            # (shell, name) -> (radial, k)
    spectra = {}
    for shell, bands in plan["bands"].items():
        chain, widths = geo[shell]
        s_r = s[np.ix_(chain[0], chain[0])]
        blocks = {t: valence_block(m, chain) for t, m in series}
        b_avg = sum(blocks.values()) / len(blocks)
        v1, lam = metric_radials(b_avg, s_r, 1)
        v2, lam2 = metric_radials(b_avg, s_r, 2)
        spectra[shell] = lam2 / lam2[0]
        priors[(shell, "1:1")] = authors_vector(s_r, widths, bands,
                                                bands[0]), 1
        priors[(shell, "neutral")] = metric_radials(
            blocks[neutral_tag], s_r, 1)[0], 1
        priors[(shell, "series1")] = v1, 1
        priors[(shell, "series2")] = v2, 2
        if shell.endswith("s") and len(series) > 2:
            b_occ = (blocks[series[0][0]] + blocks[series[1][0]]) / 2
            priors[(shell, "occ1")] = metric_radials(b_occ, s_r, 1)[0], 1
        print(f"  shell {shell} {bands}: averaged spectrum "
              f"{np.round(spectra[shell], 5)[:3]}")

    # the CASCI legs: each prior combination per cell, anchored to the
    # stored 1:1 total through the prior-independent core constant
    import json
    store = json.loads(Path(f"elements/{atom}/cells.json").read_text())
    for addr, shells, (na, nb) in plan["cells"]:
        stored = store["cells"][addr]
        stored_e, l2_target = stored["e_ground"], stored["ground_l2"]
        print(f"\n  == {addr} (stored 1:1 total {stored_e:.9f}, "
              f"term <L^2> {l2_target:.1f})")
        names = ["1:1", "neutral", "series1", "series2"]
        if all((sh, "occ1") in priors for sh in shells):
            names.append("occ1")
        for name in names:
            blocks = [assemble(nao, priors[(sh, name)][0], geo[sh][0])
                      for sh in shells]
            c_win = np.hstack(blocks)
            err0 = np.abs(c_win.T @ s @ c_win
                          - np.eye(c_win.shape[1])).max()
            # the PRODUCTION treatment: project against the core span
            # and renormalize (a raw window's core overlap contaminates
            # the active-space H at Eh scale -- measured)
            mo_full = build_full_mo(mf, occupied_core_idx(mf, NCOL),
                                    c_win)
            c_act = mo_full[:, NCOL:NCOL + c_win.shape[1]]
            err = np.abs(c_act.T @ s @ c_act
                         - np.eye(c_act.shape[1])).max()
            ops = build_operators(mol, mf, c_act, NCOL, na, nb)
            ev, vec = np.linalg.eigh(ops["H"].real)
            l2m = ops["L2"].real
            # the TARGETED TERM's root: the lowest eigenstate of the
            # stored cell's <L^2> (a bigger window can host lower
            # states of other terms -- never the comparison axis)
            sel = next(r for r in range(len(ev))
                       if abs(float(vec[:, r] @ l2m @ vec[:, r])
                              - l2_target) < 0.1)
            e_act = float(ev[sel])
            l2_0 = float(vec[:, sel] @ l2m @ vec[:, sel])
            if name == "1:1":
                const = stored_e - e_act       # the core bookkeeping
            total = e_act + const
            cert = window_closure(mol, c_act)
            print(f"    {name:>8}: total {total:.9f}  "
                  f"dE vs 1:1 {(total - stored_e) * 1e3:+8.3f} mEh  "
                  f"L2 {l2_0:7.4f}  ortho {err:.1e}  "
                  f"closure {cert['lz']:.1e}/{cert['lperp']:.1e}")
        # the capture certificates on the k=2 span
        for sh in shells:
            v2, k2 = priors[(sh, "series2")]
            c2 = assemble(nao, v2, geo[sh][0])
            cap = capture(c2, mf, s, geo[sh][0])
            print(f"    series2 {sh} capture over the reference's "
                  f"significant columns: {cap.min():.4f}-{cap.max():.4f}")
    for f in scratch.glob(f"{atom}-*.chk"):
        f.unlink()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Sc")
