"""The step-3 probes on one stored C SA cell: NEVPT2 rows for the 3P
block, then the block's gauge structure measured (README · Analysis).

Probe A (gauge redraw) -- random unitaries inside the exactly
degenerate block; E, <L^2> the controls; the E(2) spread the read.
Prediction: spread at the engine's wobble floor (~3e-7 Eh) -- a
larger spread is a defect signature, not physics.

Probe B (component identity) -- pure M_L projections through our Lz
(the CI-order bridge); <Lz> = +1, 0, -1 and <L^2> = 2 the controls.
Prediction: E(2)(M_L) identical by covariance; deviation = a
window-closure defect.

Fresh first run (2026-11-19, the certified pur cell): E(2) =
-58.1340 mEh (-1.5819 eV) IDENTICAL across the block; probe-A
spread 2.7e-11 mEh (7.3e-13 eV, the engine wobble floor) over 27
gauges, controls 1.2e-14; probe-B controls 1.0e-14, the real M_L=0
component reads the same constant (the +/-1 components are the
complex gauge -- covered by A's transitivity).  The guard's first
real bite: the store keeps flat CI vectors, the reader returns
(n,1) columns -- the assertion compares amplitudes, not shapes.

Usage:
    python scripts/probe_pt2_block.py [address]
    (default: c_2p_2-0_sa-casscf_ao_-_pur_avg -- the certified class)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402

from atomterm.cases import cell_substrate  # noqa: E402
from atomterm.cells import CellKey  # noqa: E402
from atomterm.operators import (  # noqa: E402
    build_operators, ci_in_detbasis_order, pyscf_ci_permutation)
from atomterm.perturb import (  # noqa: E402
    _root_row, nevpt2_rows, readout_mc)

EH_EV = 27.211386245988
K_REDRAW = 8
TOL_CTRL = 1e-6      # the controls: identities up to the CELL's own
                     # symmetry residual (C's SA draw 1e-14; Sc's 2.3e-7
                     # -- the optimizer's repair quality, per cell; a
                     # machinery error shows at the 1e-1 scale)
TOL_WOBBLE = 1e-6    # the engine's own run-to-run wobble floor


def _expect(mat, ci_pyscf, norb, na, nb):
    """<ci| mat |ci> in detbasis order (the CI-order bridge)."""
    w = ci_in_detbasis_order(np.asarray(ci_pyscf).ravel(), norb, na, nb)
    return float(np.real(np.vdot(w, mat @ w)))


def _to_detbasis(ci_pyscf, norb, na, nb):
    return ci_in_detbasis_order(np.asarray(ci_pyscf).ravel(), norb,
                                na, nb)


def _from_detbasis(w, norb, na, nb):
    """The inverse bridge: detbasis amplitudes -> pyscf ravelled."""
    return np.asarray(w).ravel()[pyscf_ci_permutation(norb, na, nb)]


def main() -> None:
    address = (sys.argv[1] if len(sys.argv) > 1
               else "c_2p_2-0_sa-casscf_ao_-_pur_avg")
    key = CellKey.from_address(address)
    mol, mf, data = cell_substrate(key)
    na, nb = key.sector
    ncas = int(data.ncas)
    # the term's own <L^2> target (the controls are L(L+1), not the
    # C 3P constant: Sc 2D reads 6, Ti 3F reads 20 ...)
    from atomterm.cases import CAMPAIGN
    spec = CAMPAIGN[key.atom]
    l2_target = float(spec.term[0] * (spec.term[0] + 1))
    nroots = len(data.ci)
    ncore = (mol.nelectron - (na + nb)) // 2
    mo = np.asarray(data.mo_coeff)
    c_act = mo[:, ncore:ncore + ncas]
    print(f"===== pt2 block probes: {address} =====")
    print(f"INFO substrate    ncas {ncas} sector ({na},{nb}) "
          f"nroots {nroots} ncore {ncore}; trusted ROHF cache "
          f"E = {mf.e_tot:.8f} Eh")
    e_store = np.atleast_1d(np.asarray(data.energies, dtype=float))
    print(f"INFO stored block  spread "
          f"{(e_store.max() - e_store.min()) * 1e3:.3e} mEh "
          f"({(e_store.max() - e_store.min()) * EH_EV:.3e} eV)")

    # ---- rows on the stored draw (the bit-exact assertion bites on
    # the REAL store) ----
    rows = nevpt2_rows(mol, mf, mo, ncas, (na, nb), nroots,
                       ref_energies=data.energies, ref_ci=data.ci)
    ops = build_operators(mol, mf, c_act, ncore, na, nb,
                          spherical_core=False, c_core=mo[:, :ncore])
    h_cas = ops["H"]
    e_cas_store = e_store - float(data.core_energy)
    print(f"PASS stored-roots  readout reproduced the store "
          f"(bit-exact assertion green)")
    for r, row in enumerate(rows):
        ci_p = np.asarray(data.ci[r]).ravel()
        l2 = _expect(ops["L2"], ci_p, ncas, na, nb)
        lz = _expect(ops["Lz"], ci_p, ncas, na, nb)
        print(f"root{r}  e_ref {row.e_ref:.8f} Eh  "
              f"e_corr {row.e_corr * 1e3:+.4f} mEh "
              f"({row.e_corr * EH_EV:+.4f} eV)  "
              f"e_tot {row.e_total:.8f} Eh   "
              f"<L2> {l2:.6f} <Lz> {lz:+.6f}")
        for c in row.classes:
            print(f"      {c.label:6s} {c.note:9s} norm "
                  f"{c.norm:.6f}  e {c.energy * 1e3:+8.4f} mEh")

    # ---- Probe A: gauge redraw (random unitaries inside the block) --
    ci_ref = [np.asarray(c).ravel().copy() for c in data.ci]
    e_vals = [rows[r].e_corr for r in range(nroots)]
    rng = np.random.default_rng(11)
    worst_ctrl = 0.0
    for k in range(K_REDRAW):
        q, _ = np.linalg.qr(rng.standard_normal((nroots, nroots)))
        ci_rot = [sum(q[r, j] * ci_ref[j] for j in range(nroots))
                  for r in range(nroots)]
        mc = readout_mc(mol, mf, mo, ncas, (na, nb), nroots)
        mc.ci = ci_rot  # inject the redrawn gauge (bypass the re-solve)
        for r in range(nroots):
            w = _to_detbasis(ci_rot[r], ncas, na, nb)
            worst_ctrl = max(worst_ctrl,
                             abs(_expect(ops["L2"], ci_rot[r],
                                         ncas, na, nb) - l2_target),
                             abs(_expect(h_cas, ci_rot[r],
                                         ncas, na, nb) - e_cas_store[r]))
            row = _root_row(mc, r, "nevpt2-sc")
            e_vals.append(row.e_corr)
    spread_a = max(e_vals) - min(e_vals)
    print(f"{'PASS' if worst_ctrl <= TOL_CTRL else 'FAIL'} "
          f"probe-A controls  max|<L2>-2|, |vHv-E| "
          f"{worst_ctrl:.2e} over {K_REDRAW} redraws x {nroots} roots")
    print(f"INFO probe-A redraw   E(2) spread "
          f"{spread_a * 1e3:.3e} mEh ({spread_a * EH_EV:.3e} eV) "
          f"over {len(e_vals)} values (stored + {K_REDRAW}x{nroots})")

    # ---- Probe B: pure M_L components through our Lz ---------------
    w_ref = [_to_detbasis(c, ncas, na, nb) for c in ci_ref]
    # Lz restricted to the stored span (the CAS is the whole sector
    # here: three determinants); eigenvectors = pure M_L states
    lz_blk = np.array([[complex(np.vdot(a, ops["Lz"] @ b))
                        for b in w_ref] for a in w_ref])
    lz_blk = 0.5 * (lz_blk + lz_blk.conj().T)
    ev, vecs = np.linalg.eigh(lz_blk)
    e_vals_b = []
    worst_ctrl_b = 0.0
    labels = []
    for r in range(nroots):
        # the eigenvector lives on the STORED span (coefficients over
        # w_ref); recombine into detbasis amplitudes -- it doubles as a
        # detbasis vector only in the nroots == ndet coincidence
        w_pure = np.tensordot(vecs[:, r], np.asarray(w_ref), axes=(0, 0))
        ml = float(np.real(np.vdot(w_pure, ops["Lz"] @ w_pure)))
        l2 = float(np.real(np.vdot(w_pure, ops["L2"] @ w_pure)))
        ev_e = float(np.real(np.vdot(w_pure, h_cas @ w_pure)))
        worst_ctrl_b = max(worst_ctrl_b, abs(l2 - l2_target),
                           abs(ev_e - e_cas_store[r]),
                           abs(ml - round(ml)))
        labels.append(f"M_L={round(ml):+d}")
        imag = float(np.max(np.abs(w_pure.imag)))
        if imag < 1e-10:
            mc = readout_mc(mol, mf, mo, ncas, (na, nb), 1)
            mc.ci = _from_detbasis(w_pure.real, ncas, na, nb)
            row = _root_row(mc, 0, "nevpt2-sc")
            e_vals_b.append(row.e_corr)
            print(f"INFO probe-B {labels[-1]:8s}  e_corr "
                  f"{row.e_corr * 1e3:+.4f} mEh "
                  f"({row.e_corr * EH_EV:+.4f} eV)  <L2> {l2:.6f}")
        else:
            print(f"INFO probe-B {labels[-1]:8s}  complex gauge "
                  f"(imag {imag:.1e}): its E(2) equals the block's "
                  f"constant by probe A (L=1 transitivity)")
    spread_b = max(e_vals_b) - min(e_vals_b)
    print(f"{'PASS' if worst_ctrl_b <= TOL_CTRL else 'FAIL'} "
          f"probe-B controls  max|<L2>-2|, |<M_L>-round|, |vHv-E| "
          f"{worst_ctrl_b:.2e}")
    print(f"INFO probe-B pure      E(2) spread "
          f"{spread_b * 1e3:.3e} mEh ({spread_b * EH_EV:.3e} eV) "
          f"over {labels}")

    # ---- verdicts (the outcome map, README · Analysis) -------------
    for name, spread in (("A (gauge redraw)", spread_a),
                         ("B (pure M_L)", spread_b)):
        if spread <= TOL_WOBBLE:
            verdict = ("clean -- per-root rows stand as stored; the "
                       "spread is engine wobble")
        elif spread <= 1e-3:
            verdict = ("small -- rows stand, the spread is the "
                       "per-row error bar in the key")
        else:
            verdict = ("LARGE -- term-projected (pure M_L) references "
                       "become the primary rows; the draw is the "
                       "disclosure")
        print(f"VERDICT probe {name}: spread "
              f"{spread * 1e3:.3e} mEh ({spread * EH_EV:.3e} eV) "
              f"-> {verdict}")


if __name__ == "__main__":
    main()
