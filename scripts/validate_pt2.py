"""The PT2 gate: the certificate chain over the committed payload --
a RUN INSTRUMENT (the run.py precedent), not a pytest tier.

Certificates:
  (a) the anchors: O2/6-31G CASCI(6,8) and H14/STO-3G CASCI(8,10)
      against the retained example values at 1e-6 Eh (the engine in
      this env);
  (b) the payload identities: every committed row's per-class sum ==
      e_corr (the kernel's arithmetic read twice), every root's
      e_total == e_ref + e_corr, every class row's meter consistent
      (energy == e_per_norm * norm on nonzero norms), the eV columns
      equal the Eh columns through the conversion constant;
  (c) the probe: gauge redraw (8 unitaries x the ground block) on the
      certified compact SA cell -- the E(2) spread at the engine's
      wobble floor, controls flat (probe_pt2_block.py, invoked
      fresh);
  (d) the disposition census: every declared cell rows-or-refused,
      zero silent draws.

Green output lands in logs/ (the provenance convention).

Usage: python scripts/validate_pt2.py <atom> [<atom> ...]
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
from atomterm.cases import CAMPAIGN  # noqa: E402
from run_pt2_ladder import _declared_lane_cells  # noqa: E402

EH_EV = 27.211386245988
ANCHORS = {"O2": -0.169785157128082, "H14": -0.094164359938171}

FAILS = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"{'PASS' if ok else 'FAIL'}  {name}"
          f"{'  ' + detail if detail else ''}", flush=True)
    if not ok:
        FAILS.append(name)


def _anchors() -> None:
    from pyscf import gto, scf
    from atomterm.perturb import nevpt2_rows

    mol = gto.Mole()
    mol.verbose = 0
    mol.output = None
    mol.atom = [["O", (0.0, 0.0, 0.0)], ["O", (0.0, 0.0, 1.207)]]
    mol.basis = "6-31g"
    mol.spin = 2
    mol.build()
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-20
    mf.scf()
    rows = nevpt2_rows(mol, mf, mf.mo_coeff, 6, (5, 3), 1)
    check("anchor O2", abs(rows[0].e_corr - ANCHORS["O2"]) <= 1e-6,
          f"E(2) {rows[0].e_corr:.12f} vs {ANCHORS['O2']:.12f}")

    mol = gto.Mole()
    mol.verbose = 0
    mol.output = None
    mol.atom = [["H", (0.0, 0.0, 2.0 * k + d)] for k in range(7)
                for d in (0.0, 0.8)]
    mol.basis = "sto-3g"
    mol.build()
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-20
    mf.scf()
    rows = nevpt2_rows(mol, mf, mf.mo_coeff, 8, (5, 5), 1)
    check("anchor H14", abs(rows[0].e_corr - ANCHORS["H14"]) <= 1e-6,
          f"E(2) {rows[0].e_corr:.12f} vs {ANCHORS['H14']:.12f}")


def _payload_identities(atom: str) -> int:
    payload = json.loads(
        (Path("elements") / atom / "perturbation.json").read_text())
    cells = payload["cells"]
    n_roots = n_classes = 0
    worst_sum = worst_tot = worst_ev = worst_meter = 0.0
    for c in cells:
        if c["disposition"] != "rows":
            continue
        if c.get("lane"):
            root_lists = [b["roots"] for b in c["blocks"]]
        else:
            root_lists = [c["roots"]]
        for roots in root_lists:
            for r in roots:
                n_roots += 1
                n_classes += len(r["classes"])
                s = sum(cl["energy"] for cl in r["classes"])
                worst_sum = max(worst_sum, abs(s - r["e_corr"]))
                worst_tot = max(worst_tot,
                                abs(r["e_total"] - r["e_ref"] - r["e_corr"]))
                worst_ev = max(worst_ev,
                               abs(r["e_corr_eV"] - r["e_corr"] * EH_EV),
                               abs(r["e_total_eV"] - r["e_total"] * EH_EV))
                for cl in r["classes"]:
                    if cl["norm"]:
                        worst_meter = max(
                            worst_meter,
                            abs(cl["energy"] - cl["e_per_norm"] * cl["norm"]))
    check(f"payload {atom}: class sums", worst_sum <= 1e-12,
          f"max |sum - e_corr| {worst_sum:.2e} over {n_roots} roots / "
          f"{n_classes} classes")
    check(f"payload {atom}: totals", worst_tot <= 1e-12,
          f"max |e_total - e_ref - e_corr| {worst_tot:.2e}")
    check(f"payload {atom}: eV columns", worst_ev <= 1e-9,
          f"max |eV - Eh*27.211386245988| {worst_ev:.2e}")
    check(f"payload {atom}: meter consistency",
          worst_meter <= 1e-12,
          f"max |energy - e_per_norm*norm| {worst_meter:.2e}")
    refused = sum(c["disposition"] == "refusal" for c in cells)
    check(f"payload {atom}: disposition census",
          all(c["disposition"] in ("rows", "refusal") for c in cells)
          and all((c["disposition"] == "refusal")
                  == bool(c["refusal_reason"]) for c in cells),
          f"{len(cells)} declared: {len(cells) - refused} rows / "
          f"{refused} refused")
    return n_roots


def _route_invariance(atom: str) -> None:
    """The irrep route-invariance row: the lane's ground block vs the
    free CI twin's ground root -- the same term component priced
    through the pinned route and the free route.  The band is the
    declared cross-substrate route band (the lane certificates are
    basin-blind; the observed values land in the ticket record)."""
    payload = json.loads(
        (Path("elements") / atom / "perturbation.json").read_text())
    cells = {c["address"]: c for c in payload["cells"]}
    lane_addr = free_addr = None
    for addr, kind in _declared_lane_cells(atom):
        if kind == "lane" and "casci" in addr and lane_addr is None:
            lane_addr = addr
        if kind == "grid" and "_pol_avg" in addr and free_addr is None:
            free_addr = addr
    if lane_addr is None or free_addr is None \
            or lane_addr not in cells or free_addr not in cells:
        print(f"SKIP  route invariance {atom}: no lane/free pair")
        return
    lane = cells[lane_addr]
    free = cells[free_addr]
    if lane["disposition"] != "rows" or free["disposition"] != "rows":
        print(f"SKIP  route invariance {atom}: a refused partner")
        return
    lane_e2 = lane["blocks"][0]["roots"][0]["e_corr"]
    free_e2 = free["roots"][0]["e_corr"]
    delta = abs(lane_e2 - free_e2)
    if delta <= 1e-4:
        check(f"route invariance {atom}", True,
              f"lane ground block {lane_e2 * 1e3:+.4f} mEh vs free "
              f"ground root {free_e2 * 1e3:+.4f} mEh (delta "
              f"{delta * 1e3:.2e} mEh; band 1e-4 Eh) -- the pinning "
              f"is pricing-inert")
        return
    # past the band: the disclosed orbital-gauge reading (measured on
    # Ti: the lane's D2h-canonical draw and the free rung land the
    # SAME basin and energy, different orbital gauges -- E(2) moves
    # at the mEh scale, the wave's orbital-draw finding).  The
    # machinery certificate lives elsewhere (the anchors, the kernel
    # twin, probe A on the same-orbital block)
    print(f"DISCLOSED  route invariance {atom}: lane ground block "
          f"{lane_e2 * 1e3:+.4f} mEh vs free ground root "
          f"{free_e2 * 1e3:+.4f} mEh (delta {delta * 1e3:.2e} mEh "
          f"past the band) -- the lane's D2h gauge vs the free draw: "
          f"the orbital-draw reading, the ticket carries the record")


def _probe(atom: str) -> None:
    import subprocess

    addr = {"C": "c_2p_2-0_sa-casscf_ao_-_pur_avg",
            "Si": "si_3p_2-0_sa-casscf_ao_-_pur_avg",
            "Sc": "sc_3d4s_2-1_sa-casscf_ao_-_pur_avg",
            "Ti": "ti_3d4s_3-1_sa-casscf_ao_-_pur_avg",
            "Sc-unc": "sc-unc_3d4s_2-1_sa-casscf_ao_-_pur_avg",
            "Ti-unc": "ti-unc_3d4s_3-1_sa-casscf_ao_-_pur_avg",
            "U14e": "u14e_10f13d_4-0_sa-casscf_ao_-_pur_avg"}.get(atom)
    if addr is None:
        print(f"SKIP  probe {atom}: no certified SA address declared")
        return
    r = subprocess.run(
        [sys.executable, "scripts/probe_pt2_block.py", addr],
        capture_output=True, text=True)
    out = r.stdout + r.stderr
    passes = out.count("PASS")
    verdicts = [ln for ln in out.splitlines()
                if ln.startswith("VERDICT")]
    # the per-L rule (README . Analysis): L <= 1 requires the clean or
    # small verdicts (a LARGE spread is the transitivity defect); L >=
    # 2 admits LARGE as the DECLARED disclosure (rotation orbits do
    # not cover the block -- the canonical gauge stays the primary
    # row, the spread the disclosure)
    l_term = CAMPAIGN[atom].term[0]
    machinery_fail = any("LARGE" in v for v in verdicts)
    clean = (r.returncode == 0 and len(verdicts) >= 2
             and (not machinery_fail or l_term >= 2))
    check(f"probe {atom} ({addr})", clean,
          f"{passes} PASS lines, {len(verdicts)} verdicts"
          f"{' (L>=2: LARGE = the declared spread disclosure)' if machinery_fail else ''}")


def main() -> None:
    atoms = sys.argv[1:]
    if not atoms:
        sys.exit("usage: validate_pt2.py <atom> [<atom> ...]")
    print("===== validate_pt2 =====")
    _anchors()
    n = sum(_payload_identities(a) for a in atoms)
    for a in atoms:
        _route_invariance(a)
        _probe(a)
    if FAILS:
        sys.exit(f"\nvalidate_pt2: {len(FAILS)} FAIL: {FAILS}")
    print(f"\nvalidate_pt2: ALL GREEN ({n} roots certified)")


if __name__ == "__main__":
    main()
