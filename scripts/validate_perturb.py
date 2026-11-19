"""Validate the perturbation analyzation (`perturb.py`) at campaign
sizes: the three-layer certificate chain printed as certificate lines
on the stored C/Si SA roots.  Its green output is the GATE the
PERTURBATION_TABLE render cites (spec: wayfinder/specs/
pt2-perturbation.md; tickets 43 -> 45).  A run instrument (the
``run.py`` precedent), not a pytest tier -- the campaign-size
correctness evidence is priced at doors.

Scope (the declared class, the ``regen_cells.py`` grammar): the
``grid_keys`` + ``planned`` cells of the given atoms, ``--method``
(default ``sa-casscf``), ``--start`` (default ``ao`` -- the
m-complete kind; the MO/MOC/NO controls join by declaration), free
route only (core-tagged lane cells run on the D2h substrate -- not
the PT2 reference class).

The chain, per cell:

*  (a) twins at campaign size -- the production coupling column vs
   PySCF's own C-contraction kernels (:func:`operators.pyscf_apply`,
   an independent construction through a different code path), and
   the EN E(2) production vs the dense Slater-Condon twin;
*  (b) the external anchor for the REFERENCE energy -- the
   recomputed ``e0`` vs the store's PySCF CI energy, on the run's
   OWN core field (the raw density of the run's core prefix: the
   stored energies are the kernel's values on exactly that field).
   The E(2) external anchor is ABSENT from the env (pyscf 2.13 ships
   no NEVPT2/CASPT2); the note prints with the summary -- the
   pyscf-forge install stays the open user-gated option;
*  (c) exact identities -- the whole-space anchor (CAS = valence
   space => E(2) = 0 exactly), the FOIS completeness certificate
   (the production coupled pattern == the brute substitution count,
   on a synthetic full-support CI -- the stored roots' own support is
   gauge, disclosed per root), and the stored gauges' containment
   (no coupling outside the Slater-Condon reach).

Usage:
    python scripts/validate_perturb.py [atoms] [--method sa-casscf]
        [--start ao]
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
from pyscf import ao2mo  # noqa: E402

from atomterm.cases import CAMPAIGN, campaign_spec, scenarios_for  # noqa: E402
from atomterm.cells import (  # noqa: E402
    _key_core, grid_keys, load_cell, planned, void_reason)
from atomterm.detbasis import (  # noqa: E402
    determinant_basis, operator_apply, slater_condon_matrix)
from atomterm.operators import (  # noqa: E402
    active_space_integrals, ci_in_detbasis_order, pyscf_apply,
    spin_block_one_body, spin_block_two_body)
from atomterm.perturb import (  # noqa: E402
    epstein_nesbet, epstein_nesbet_dense, dyall_e2, dyall_e2_dense,
    dyall_h0, fois_indices, rdm_cumulants, rdm_reconstruct, spin_rdm)
from atomterm.reference import run_rohf  # noqa: E402

# the certificate tolerances (the observed fresh values sit two to six
# orders below each; a graze is a finding, not a band to widen)
TOL_E0 = 1e-6        # observed fresh with the run's raw fold: 1e-15
                     # (the averaged-field fold sits at 5e-05 -- the
                     # ruler the anchor exists to catch)
TOL_COUPLE = 1e-10   # observed fresh: 4.4e-16
TOL_DENSE = 1e-10    # observed fresh: 0.0 (same-arithmetic agreement)
DENSE_BUDGET = 20000  # the dense twin is quadratic in the sector


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="the perturbation analyzation's campaign-size "
                    "certificate chain (the render gate)")
    p.add_argument("atoms", nargs="*", default=["C", "Si"],
                   help="campaign atoms (default: C Si -- the wave-1 "
                        "declarations)")
    p.add_argument("--method", default="sa-casscf",
                   help="the cell method axis (default sa-casscf)")
    p.add_argument("--start", default="ao",
                   help="the cell start axis (default ao)")
    return p.parse_args(argv)


def scoped_keys(atom, spec, args):
    """The declared validation scope: (key, status) pairs -- status is
    ``None`` (runnable), a void reason string, or ``"no chkfile"``."""
    scens = scenarios_for(spec)
    out = []
    seen = set()
    for k in grid_keys(spec) + planned(spec):
        if (k.method != args.method or k.start != args.start
                or k.core_tag is not None or k.address in seen):
            continue
        seen.add(k.address)
        if Path(f"elements/{atom}/rungs/{k.address}.chk").exists():
            out.append((k, None))
        elif k.source_tag is not None and k.mask == "radial":
            # the declared construction refusal (run_cell raises it):
            # the member-labeled prior substrate has no radial meaning
            out.append((k, "construction-refused (the prior-window "
                          "radial mask is undeclared)"))
        else:
            v = void_reason(k, scens)
            out.append((k, f"void ({v})" if v else "no chkfile"))
    return out


def cell_tensors(key, mol, mf, spec, data):
    """The valence-space substrate of a stored cell: the folded-core
    one-body + raw physicist tensor over the run's OWN columns, read
    the way the run produced them -- ``c_core`` is the run's own core
    prefix and the field is the RAW core density
    (``spherical_core=False``): the stored CI energies are PySCF's
    kernel values on exactly that field, and the anchor must read the
    run as stored (fact 14: optimized stages fold their OWN inactive
    columns; the optimized core's gauge can carry a p-admixture whose
    anisotropy moves ``e0`` at the 1e-4 Eh scale).  Spin-blocked into
    the detbasis convention; the spatial pair returned for
    ``pyscf_apply``."""
    na, nb = key.sector
    ncore, _, _ = _key_core(key, mol, mf, spec.shell_bands)
    ncas = data.ncas
    mo = np.asarray(data.mo_coeff)
    c_val = mo[:, ncore:]
    n_val = c_val.shape[1]
    h1_s, g_s = active_space_integrals(mol, mf, c_val, ncore,
                                       c_core=mo[:, :ncore],
                                       spherical_core=False)
    return (ncore, ncas, n_val,
            spin_block_one_body(h1_s), spin_block_two_body(g_s),
            h1_s, ao2mo.kernel(mol, c_val))


class Certificates:
    """The run's certificate counters: one home for the PASS/FAIL
    tally and the whole-space anchor's once-per-atom placement."""

    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.whole_space_done = False

    def line(self, label, ok):
        """One certificate line + its tally."""
        print(f"{'PASS' if ok else 'FAIL'} {label}", flush=True)
        if ok:
            self.passed += 1
        else:
            self.failed += 1

    @property
    def total(self):
        return self.passed + self.failed


def validate_cell(key, mol, mf, spec, data, certs):
    """One cell's certificates (the lines + the tally live in
    ``certs``)."""
    na, nb = key.sector
    ncore, ncas, n_val, h1v, gv, h1_s, eri_pyscf = cell_tensors(
        key, mol, mf, spec, data)
    dets = determinant_basis(n_val, na, nb)
    ndet = len(dets)
    cas_idx = {k for k, d in enumerate(dets)
               if all(o < ncas or n_val <= o < n_val + ncas for o in d)}
    brute = set(fois_indices(n_val, na, nb, ncas).tolist())
    energies = np.atleast_1d(data.energies)

    print(f"---- {key.address}  n_val={n_val} ncas={ncas} ndet={ndet} "
          f"nroots={len(data.ci)}", flush=True)

    # (c) FOIS completeness: production vs dense ON A SYNTHETIC
    # FULL-SUPPORT CAS CI -- two independent constructions must agree
    # exactly (n_fois AND e2).  The brute substitution count is the
    # structural REACH (an upper bound): symmetry-adapted tensors
    # (the pur substrates) carry selection-rule zeros inside it -- a
    # declared class, disclosed on the INFO line, never a failure.
    rng = np.random.default_rng(0)
    ci_full = rng.standard_normal(len(cas_idx))
    ci_full /= np.linalg.norm(ci_full)
    res_synth = epstein_nesbet(h1v, gv, ncas, na, nb, ci_full)
    den_synth = epstein_nesbet_dense(h1v, gv, ncas, na, nb, ci_full)
    d_syn = float(np.abs(res_synth.e2 - den_synth.e2))
    certs.line(
        f"completeness production n_fois {res_synth.n_fois} == dense "
        f"{den_synth.n_fois}, e2 diff {d_syn:.1e}",
        res_synth.n_fois == den_synth.n_fois and d_syn <= TOL_DENSE)
    v_synth = np.zeros(ndet)
    v_synth[list(cas_idx)] = ci_full
    hv_synth = operator_apply(dets, one_body=h1v, two_body=gv)(v_synth)
    coupled_synth = {k for k in np.nonzero(hv_synth != 0)[0].tolist()
                     if k not in cas_idx}
    certs.line(
        f"reach-containment synthetic coupled {len(coupled_synth)}, "
        f"all in the brute reach {len(brute)}",
        coupled_synth <= brute)
    print(f"INFO reach          coupled {res_synth.n_fois} of the "
          f"brute reach {len(brute)} (selection-rule zeros on "
          f"symmetry-adapted tensors are declared class)", flush=True)

    # (c) whole-space anchor, once per atom on the first cell's
    # tensors: CAS = valence => the FOIS is empty => E(2) = 0 exact;
    # then the density layer on the same cell: the RDM stack's
    # conventions against the Hamiltonian (the energy identity), the
    # in-repo gamma anchor, the contraction chain, and the cumulant
    # reconstruction
    if not certs.whole_space_done:
        ci_ws = np.random.default_rng(1).standard_normal(ndet)
        ci_ws /= np.linalg.norm(ci_ws)
        res_ws = epstein_nesbet(h1v, gv, n_val, na, nb, ci_ws)
        certs.line(
            f"whole-space  E(2) = {float(np.real(res_ws.e2)):.1e} exact "
            f"(CAS = the valence space)",
            res_ws.e2 == 0.0 and res_ws.n_fois == 0)
        certs.whole_space_done = True

        ci0 = ci_in_detbasis_order(np.asarray(data.ci[0]), ncas, na,
                                   nb)
        g1 = spin_rdm(ci0, ncas, na, nb, 1)
        G2 = spin_rdm(ci0, ncas, na, nb, 2)
        emb_sp = np.arange(ncas)
        emb_so = np.concatenate([np.arange(ncas),
                                 n_val + np.arange(ncas)]).astype(int)
        h1c = spin_block_one_body(h1_s[np.ix_(emb_sp, emb_sp)])
        gc = gv[np.ix_(emb_so, emb_so, emb_so, emb_so)]
        H_cas = slater_condon_matrix(determinant_basis(ncas, na, nb),
                                     one_body=h1c, two_body=gc)
        e0_cas = float(np.real(np.vdot(
            ci0, H_cas @ ci0)))
        e_id = float(np.real(np.sum(h1c * g1) + 0.5 * np.sum(gc * G2)))
        certs.line(
            f"density       energy identity diff "
            f"{abs(e_id - e0_cas):.1e} (E0 = h1:g1 + G:G2/2)",
            abs(e_id - e0_cas) <= 1e-9)
        from atomterm.natural import ci_rdm1
        d_nat = float(np.max(np.abs(
            g1[:ncas, :ncas] + g1[ncas:, ncas:]
            - ci_rdm1(determinant_basis(ncas, na, nb), ci0, ncas))))
        certs.line(
            f"density       gamma anchor   maxdiff {d_nat:.1e} (vs "
            f"natural.ci_rdm1)", d_nat <= 1e-9)
        d_c = float(np.max(np.abs(np.trace(G2, axis1=1, axis2=2)
                                  + g1)))
        certs.line(
            f"density       contraction    maxdiff {d_c:.1e} "
            f"(G2 middle-pair trace == -gamma)", d_c <= 1e-9)
        try:
            G3 = spin_rdm(ci0, ncas, na, nb, 3)
            lam = rdm_cumulants({1: g1, 2: G2, 3: G3}, 3, 2 * ncas)
            d_r = float(np.max(np.abs(
                rdm_reconstruct(lam, 3, 2 * ncas) - G3)))
            certs.line(
                f"density       reconstruction maxdiff {d_r:.1e} "
                f"(Gamma3 = lambda3 + wedges)", d_r <= 1e-9)
        except ValueError as exc:
            certs.line(f"density       cumulant stack REFUSED: {exc}",
                       False)

    apply_pyscf = pyscf_apply(h1_s, eri_pyscf, n_val, na, nb)
    v = np.zeros(ndet)

    dense_ok_size = ndet <= DENSE_BUDGET
    for r, ci_pyscf in enumerate(np.asarray(data.ci)):
        ci = ci_in_detbasis_order(np.asarray(ci_pyscf), ncas, na, nb)
        res = epstein_nesbet(h1v, gv, ncas, na, nb, ci)
        v[:] = 0.0
        v[list(cas_idx)] = ci
        hv_ours = operator_apply(dets, one_body=h1v, two_body=gv)(v)
        hv_pys = apply_pyscf(v)
        coupled = {k for k in np.nonzero(hv_ours != 0)[0].tolist()
                   if k not in cas_idx}

        # (b) external anchor for the reference energy
        e_ref = float(energies[r]) - float(data.core_energy)
        d_e0 = abs(float(np.real(res.e0)) - e_ref)
        certs.line(
            f"root{r} e0-anchor    diff {d_e0:.1e} (tol {TOL_E0:.0e}; "
            f"store PySCF CI energy)", d_e0 <= TOL_E0)

        # (a) coupling twin at campaign size
        d_c = float(np.max(np.abs(hv_ours - hv_pys)))
        certs.line(
            f"root{r} coupling      maxdiff {d_c:.1e} "
            f"(tol {TOL_COUPLE:.0e}; vs pyscf_apply)",
            d_c <= TOL_COUPLE)

        # (a) the dense Slater-Condon twin
        den = epstein_nesbet_dense(h1v, gv, ncas, na, nb, ci)
        d_d = abs(float(np.real(res.e2 - den.e2)))
        certs.line(
            f"root{r} dense-twin    diff {d_d:.1e} "
            f"(tol {TOL_DENSE:.0e}; EN production vs dense)",
            d_d <= TOL_DENSE)

        # containment: no stored-gauge coupling outside the reach
        certs.line(
            f"root{r} containment   {len(coupled)} coupled, all in the "
            f"brute FOIS ({len(brute)})", coupled <= brute)

        support = int(np.count_nonzero(ci))
        print(f"INFO root{r} gauge        support {support}/"
              f"{len(cas_idx)} dets, max|c| "
              f"{np.max(np.abs(ci)):.3f} (degenerate-block gauge is "
              f"chaotic -- the meters read the stored draw)",
          flush=True)
        print(f"INFO root{r} EN          e2 {float(np.real(res.e2)) * 1e3:+.3f} "
              f"mEh  floor {res.floor:.4f}  w_ref {res.w_ref:.6f}  "
              f"n_fois {res.n_fois}", flush=True)

        # rung 2: the Dyall partition + correction on the same root
        try:
            h0 = dyall_h0(h1v, gv, ncas, na, nb, ci)
            res_d = dyall_e2(h0, ncas, na, nb, ci)
            res_dd = dyall_e2_dense(h0, ncas, na, nb, ci)
            d_dd = float(np.abs(res_d.e2 - res_dd.e2))
            certs.line(
                f"root{r} dyall         e0 anchor "
                f"{abs(float(np.real(res_d.e0)) - e_ref):.1e} "
                f"(tol {TOL_E0:.0e}); twin diff {d_dd:.1e} "
                f"(tol {TOL_DENSE:.0e})",
                abs(float(np.real(res_d.e0)) - e_ref) <= TOL_E0
                and d_dd <= TOL_DENSE)
            print(f"INFO root{r} Dyall       e2 "
                  f"{float(np.real(res_d.e2)) * 1e3:+.3f} mEh  "
                  f"floor {res_d.floor:.4f}  w_ref "
                  f"{res_d.w_ref:.6f}  n_fois {res_d.n_fois}",
                  flush=True)
        except ValueError as exc:
            certs.line(f"root{r} dyall         REFUSED: {exc}", False)
    return


def main() -> None:
    args = parse_args(sys.argv[1:])
    unknown = [a for a in args.atoms if a not in CAMPAIGN]
    if unknown:
        sys.exit(f"not declared: {', '.join(unknown)}; "
                 f"declared: {', '.join(CAMPAIGN)}")
    total_fail = 0
    for atom in args.atoms:
        spec = campaign_spec(atom)
        cfg = spec.atom_config()
        mol, mf = run_rohf(cfg, f"elements/{atom}/rungs/rohf.chk",
                           ref=spec.ref)
        if not mf.converged:
            sys.exit(f"{atom}: the reference ROHF did not converge")
        print(f"\n===== perturb validation: {atom} =====", flush=True)
        print(f"PASS substrate    ROHF E = {mf.e_tot:.8f} Eh "
              f"(trusted cache elements/{atom}/rungs/rohf.chk)",
              flush=True)
        keys = scoped_keys(atom, spec, args)
        runnable = [(k, s) for k, s in keys if s is None]
        print(f"NOTE scope        method={args.method} start={args.start} "
              f"free-route: {len(runnable)} runnable of {len(keys)} "
              f"declared; lane-class (core-tagged) cells excluded (the "
              f"D2h substrate is not the PT2 reference class)",
          flush=True)
        stats = Certificates()
        for key, status in keys:
            if status is not None:
                print(f"SKIP {key.address}  {status}", flush=True)
                continue
            data = load_cell(key)
            if data.ncas is None:
                print(f"SKIP {key.address}  no ncas in the store",
                      flush=True)
                continue
            validate_cell(key, mol, mf, spec, data, stats)
        print(f"NOTE external-anchor absent from the env for E(2): "
              f"pyscf 2.13 ships no NEVPT2/CASPT2; the pyscf-forge "
              f"install stays the open user-gated option (the spec "
              f"records it)", flush=True)
        print(f"===== {atom}: {stats.passed} PASS, {stats.failed} "
              f"FAIL =====", flush=True)
        total_fail += stats.failed
    if total_fail:
        sys.exit(f"\n{total_fail} certificate(s) FAILED -- the render "
                 "gate is RED (ticket 45 must not run)")


if __name__ == "__main__":
    main()
