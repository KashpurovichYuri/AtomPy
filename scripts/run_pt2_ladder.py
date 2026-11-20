"""The PT2 ladder collector: NEVPT2 rows over the DECLARED cell set,
written to ``elements/<Atom>/perturbation.json`` -- the renderer's
only source (the render contract: the table and its rows cannot
diverge).

THE DECLARED SET (the intent, one manifest per atom -- stale cells
cannot leak in; a missing address refuses before any run):
  the span ladder     CASCI {compact, @r1, @s1[, @s2]} x {ao, mo, moc}
  the NO twins        CASCI per declared window tail
  the reference axis  {casscf, sa-casscf} x declared windows x
                      {ao, mo}
  the frozen pair     the core-relaxation disclosure rows
  the CAS-size axis   U14e's promoted trio (10f13d / 10f13d7s /
                      10f13d14p) -- the flagship R-test (Ti's same-l
                      axis is VOID-declared: the distilled same-l
                      windows wait on window_distilled's future seam)
  the lane            the core-partition axis: the pinned CI lane +
                      its pinned SA twin + the free CI twin (the
                      route-invariance partner) per declared leg;
                      priced through the lane store's own readout
                      twin (cells.readout_lane) -- the lane payload
                      class carries no ci_stack/core_energy
Every stored root is a row (the declared-roots policy -- no silent
draws).  Per-root meters (L2 / Lz / S2) are measured LIVE through the
CI-order bridge; cells.json is never a source.  A stored-roots
refusal does not abort the ladder: the cell records its refusal
(reason string) and the ladder moves on -- the recovery loop is
``scripts/regen_cells.py`` scoped to the refused address, then a
collector re-run fills the row.

Usage:
    python scripts/run_pt2_ladder.py <atom> [<atom> ...]
"""

import json
import resource
import sys
import time
from datetime import date, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
from pyscf import __version__ as pyscf_version  # noqa: E402

from atomterm.cases import campaign_spec, cell_substrate, lane_substrate  # noqa: E402
from atomterm.cells import CellKey, cell_chk  # noqa: E402
from atomterm.detbasis import determinant_basis, operator_apply  # noqa: E402
from atomterm.operators import (ci_in_detbasis_order, l2_operator_parts,  # noqa: E402
                                lz_operator, s2_operator_parts)
from atomterm.perturb import nevpt2_rows  # noqa: E402

EH_EV = 27.211386245988


def _declared(atom: str) -> dict:
    """The declared set per atom (the manifest -- edit here).

    ``windows``      the span ladder's windows (compact first);
    ``no_tail``      the NO twins' suffix (``no1``/``no2``/...);
    ``no_windows``   explicit NO-twin windows (U14e's own pattern);
    ``size``         the CAS-size axis (window -> its NO-tail suffix),
                     the flagship R-test;
    ``dist``         the distilled prior family (casci ao);
    ``window_sector``  per-window sector overrides (else ``sector``);
    ``sector``       the ladder's primary sector (the maximal-Ms class
                     the compact cells carry);
    ``ref_windows``  the windows the reference axis declares (the unc
                     twins carry CI-only @s tails);
    ``frozen``       the frozen pair's windows;
    ``lane_sectors`` the lane legs -> their own core-derived sector;
    ``lane_sa``      the legs carrying a pinned SA lane cell.
    """
    return {
        "C": dict(
            windows=("c_2p", "c_2p@r1", "c_2p@s1"), no_tail="no1",
            sector="2-0", frozen=("c_2p", "c_2p@r1"),
            lane_sectors={"c_2p@1s2s": "2-0"},
            lane_sa=("c_2p@1s2s",)),
        "Si": dict(
            windows=("si_3p", "si_3p@r1", "si_3p@s1", "si_3p@s2"),
            no_tail="no1", sector="2-0",
            frozen=("si_3p", "si_3p@r1"),
            lane_sectors={"si_3p@1s2s2p3s": "2-0"},
            lane_sa=("si_3p@1s2s2p3s",)),
        "Sc": dict(
            windows=("sc_3d4s", "sc_3d4s@r1", "sc_3d4s@s1", "sc_3d4s@s2"),
            no_tail="no2", sector="2-1",
            frozen=("sc_3d4s", "sc_3d4s@r1"),
            lane_sectors={"sc_3d4s@ar": "2-1", "sc_3d@ar": "2-1",
                          "sc_3d3p@semi": "4-3", "sc_3d3p4d@semi": "4-3",
                          "sc_3d3p4s@1s2s2p3s": "5-4",
                          "sc_3s3p3d@1s2s2p4s": "5-4"},
            lane_sa=("sc_3d4s@ar", "sc_3d@ar", "sc_3d3p@semi",
                     "sc_3d3p4d@semi", "sc_3d3p4s@1s2s2p3s",
                     "sc_3s3p3d@1s2s2p4s")),
        "Ti": dict(
            windows=("ti_3d4s", "ti_3d4s@r1", "ti_3d4s@s1", "ti_3d4s@s2"),
            no_tail="no2", sector="3-1",
            # the CAS-size axis (3d4d / 3d4s4p4d) is VOID-declared at
            # pur: the same-l distilled windows wait on the
            # member-preserving cross-shell orthogonalization
            # (window_distilled's future seam; since the distillation
            # wave 2026-11-10) -- the CAS-size axis is U14e's trio
            frozen=("ti_3d4s", "ti_3d4s@r1"),
            # four lane legs OUT (the fresh construction refuses on
            # today's windows.py; their stores are schema-stale, no
            # win_lz): the same-l pairs 3d4d@ar / 3d4s4p4d@ar (the
            # distilled window's S-orthonormality -- the declared
            # future seam) and the semi legs 3d3p@semi / 3d3p4d@semi
            # (the Loewdin needle in the raw shell span) -- the rows
            # wait on a windows.py resolution, recorded in ticket 45
            lane_sectors={"ti_3d4s@ar": "3-1", "ti_3d@ar": "3-1",
                          "ti_3s3p3d@1s2s2p4s": "6-4",
                          "ti_3d3p4s@1s2s2p3s": "6-4"},
            lane_sa=("ti_3d4s@ar", "ti_3d@ar",
                     "ti_3s3p3d@1s2s2p4s", "ti_3d3p4s@1s2s2p3s")),
        "Sc-unc": dict(
            windows=("sc-unc_3d4s", "sc-unc_3d4s@r1", "sc-unc_3d4s@s1",
                     "sc-unc_3d4s@s2"),
            no_tail="no2", sector="2-1",
            ref_windows=("sc-unc_3d4s", "sc-unc_3d4s@r1"),
            frozen=("sc-unc_3d4s", "sc-unc_3d4s@r1"),
            # the lane legs OUT: the lane substrate's basin steer
            # cannot build today (the plain D2h draw lands a wrong
            # basin and the C1-seeded retry diverges) -- recorded in
            # ticket 45; the rows wait on a lane_substrate resolution
            lane_sectors={}, lane_sa=()),
        "Ti-unc": dict(
            windows=("ti-unc_3d4s", "ti-unc_3d4s@r1", "ti-unc_3d4s@s1",
                     "ti-unc_3d4s@s2"),
            no_tail="no2", sector="3-1",
            ref_windows=("ti-unc_3d4s", "ti-unc_3d4s@r1"),
            frozen=("ti-unc_3d4s", "ti-unc_3d4s@r1"),
            lane_sectors={}, lane_sa=()),  # the substrate steer: see Sc-unc
        "U14e": dict(
            windows=("u14e_10f13d", "u14e_10f13d7s", "u14e_10f13d14p"),
            size_group=("u14e_10f13d", "u14e_10f13d7s", "u14e_10f13d14p"),
            sector="4-0",
            window_sector={"u14e_10f13d7s": "5-1"},
            no_windows=("u14e_10f13dno2", "u14e_10f13d14pno3"),
            dist=("u14e_5f6d@r1", "u14e_5f6d@s1"),
            frozen=("u14e_10f13d",),
            lane_sectors={"u14e_5f6d@6s6p7s": "4-0",
                          "u14e_5f6d7s@6s6p": "5-1"},
            lane_sa=("u14e_5f6d@6s6p7s", "u14e_5f6d7s@6s6p")),
    }[atom]


def _declared_addresses(atom: str) -> list[str]:
    """The declared cell set per atom (the grid rows; the lane legs
    live in :func:`_declared_lane_cells`)."""
    d = _declared(atom)
    windows = d["windows"]

    def sect_of(x: str) -> str:
        return d.get("window_sector", {}).get(x, d["sector"])

    addrs = []
    for x in windows:
        addrs += [f"{x}_{sect_of(x)}_casci_{s}_-_pur_avg"
                  for s in ("ao", "mo", "moc")]
    if d.get("no_tail"):
        head = windows[0]
        addrs.append(f"{head}{d['no_tail']}_{sect_of(head)}_casci_no_-_pur_avg")
        addrs += [f"{head}@{w.split('@', 1)[1]}{d['no_tail']}"
                  f"_{sect_of(head)}_casci_no_-_pur_avg" for w in windows[1:]]
    if d.get("no_windows"):
        addrs += [f"{x}_{sect_of(x)}_casci_no_-_pur_avg"
                  for x in d["no_windows"]]
    if d.get("dist"):
        for x in d["dist"]:
            addrs.append(f"{x}_{sect_of(x)}_casci_ao_-_pur_avg")
    ref_windows = d.get("ref_windows", windows)
    for x in ref_windows:
        addrs += [f"{x}_{sect_of(x)}_{m}_{s}_-_pur_avg"
                  for m in ("casscf", "sa-casscf") for s in ("ao", "mo")]
    for x in d["frozen"]:
        addrs.append(f"{x}_{sect_of(x)}_sa-casscf_ao_frozen_pur_avg")
    return addrs


def _declared_lane_cells(atom: str) -> list[tuple[str, str]]:
    """The declared lane legs as ``(address, kind)``: the pinned CI
    lane (``pur``), its pinned SA twin where the store carries one,
    and the free CI twin (``pol`` -- the route-invariance partner,
    priced through the grid seam)."""
    d = _declared(atom)
    out = []
    for leg, sect in d["lane_sectors"].items():
        out.append((f"{leg}_{sect}_casci_ao_-_pur_avg", "lane"))
        if leg in d.get("lane_sa", ()):
            out.append((f"{leg}_{sect}_sa-casscf_ao_-_pur_avg", "lane"))
        out.append((f"{leg}_{sect}_casci_ao_-_pol_avg", "grid"))
    return out


def _meters(mol, mf, c_act, ncore, na, nb):
    """The meter machinery: the operators as their ACTION
    (``detbasis.operator_apply`` -- the matrix-free twin; the
    ``len(dets)^2`` container never forms, any sector size).  The
    many-body meters read per root through the CI-order bridge."""
    dets = determinant_basis(c_act.shape[1], na, nb)
    apply_l2 = operator_apply(dets, *l2_operator_parts(mol, c_act))
    apply_s2 = operator_apply(dets, *s2_operator_parts(c_act.shape[1]))
    apply_lz = operator_apply(dets, one_body=lz_operator(mol, c_act))

    def expect(kind, ci_pyscf):
        w = ci_in_detbasis_order(np.asarray(ci_pyscf).ravel(),
                                 c_act.shape[1], na, nb)
        if kind == "l2":
            return float(np.real(w @ apply_l2(w)))
        if kind == "s2":
            return float(np.real(w @ apply_s2(w)))
        return float(np.real(w @ apply_lz(w)))

    return expect


def _collect_lane(key, lane_sub, spec, workdir, t0) -> dict:
    """One lane cell's entry: the stored-roots twin's re-solve
    (asserted) + the per-block PT2 rows.  The meters ride the stored
    payload (the twin asserted them); the block spread is the term's
    across-component E(2) ptp."""
    from atomterm.cells import readout_lane
    from atomterm.perturb import nevpt2_rows_on

    mol, mf = lane_sub
    entry = {
        "address": key.address, "lane": True,
        "sector": list(key.sector),
        "term": None, "blocks": [],
        "disposition": "rows", "refusal_reason": None,
    }
    try:
        payload, blocks = readout_lane(
            key, mol, mf, shell_bands=spec.shell_bands,
            shell_radials=spec.shell_radials, workdir=workdir)
    except ValueError as exc:
        entry["disposition"] = "refusal"
        entry["refusal_reason"] = str(exc)
        print(f"  [refused] {key.address}: {exc}", flush=True)
        return entry
    entry["term"] = list(payload["term"])
    entry["core_energy"] = float(payload["total"]) - float(
        np.mean(payload["energies"]))
    e2s = []
    i0 = 0  # the blocks' global meter offset (the payload's l2/s2
    # lists run in branch order across ALL components -- the same
    # accumulation readout_lane asserted against)
    for irrep, mc in blocks.items():
        mult = int(mc.fcisolver.nroots)
        blk = {"irrep": irrep, "mult": mult, "roots": []}
        for r in nevpt2_rows_on(mc):
            i = i0 + r.root
            e2s.append(r.e_corr)
            blk["roots"].append({
                "root": r.root,
                "e_ref": r.e_ref,
                "e_corr": r.e_corr,
                "e_total": r.e_total,
                "e_corr_eV": r.e_corr * EH_EV,
                "e_total_eV": r.e_total * EH_EV,
                "l2": float(payload["l2"][i]),
                "s2": float(payload["s2"][i]),
                "classes": [{
                    "label": c.label, "note": c.note,
                    "norm": c.norm, "energy": c.energy,
                    "energy_meh": c.energy * 1e3,
                    "e_per_norm": c.e_per_norm,
                } for c in r.classes],
            })
        entry["blocks"].append(blk)
        i0 += mult
    entry["term_e2_mean"] = float(np.mean(e2s))
    entry["term_e2_ptp"] = float(np.ptp(e2s))
    print(f"  [lane] {key.address}: {len(blocks)} blocks, "
          f"E(2) mean {entry['term_e2_mean'] * 1e3:+.4f} mEh "
          f"({entry['term_e2_mean'] * EH_EV:+.4f} eV), ptp "
          f"{entry['term_e2_ptp'] * 1e3:.2e} mEh "
          f"({time.time() - t0:.1f} s)", flush=True)
    return entry


def collect_atom(atom: str, workdir: str = "elements") -> dict:
    """One atom's declared ladder: rows + meters + provenance."""
    declared = _declared_addresses(atom) + \
        [a for a, _ in _declared_lane_cells(atom)]
    missing = [a for a in declared
               if not cell_chk(CellKey.from_address(a), workdir).exists()]
    if missing:
        sys.exit(f"{atom}: {len(missing)} declared addresses are not in "
                 f"the store: {missing[:3]}{'...' if len(missing) > 3 else ''}")
    spec = campaign_spec(atom)
    lane_sub = None
    entries = []
    for addr in declared:
        key = CellKey.from_address(addr)
        na, nb = key.sector
        t0 = time.time()
        is_lane = dict(_declared_lane_cells(atom)).get(addr) == "lane"
        if is_lane:
            if lane_sub is None:
                lane_sub = lane_substrate(spec)
            entries.append(_collect_lane(key, lane_sub, spec, workdir,
                                         t0))
            continue
        mol, mf, data = cell_substrate(key, workdir)
        ncas = int(data.ncas)
        nroots = len(data.ci)
        ncore = (mol.nelectron - (na + nb)) // 2
        mo = np.asarray(data.mo_coeff)
        c_act = mo[:, ncore:ncore + ncas]
        expect = _meters(mol, mf, c_act, ncore, na, nb)
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
        print(f"  [rss] {addr} at {rss:.1f} GB", flush=True)
        entry = {
            "address": addr,
            "ncas": ncas, "nroots": nroots,
            "sector": [na, nb],
            "core_energy": float(data.core_energy),
            "kernel_energies": np.asarray(
                data.energies, dtype=float).ravel().tolist(),
            "converged": data.converged,
            "cycles": data.cycles,
            "roots": [], "disposition": "rows", "refusal_reason": None,
        }
        try:
            rows = nevpt2_rows(mol, mf, mo, ncas, (na, nb), nroots,
                               ref_energies=data.energies, ref_ci=data.ci)
        except ValueError as exc:
            entry["disposition"] = "refusal"
            entry["refusal_reason"] = str(exc)
            print(f"  [refused] {addr}: {exc}", flush=True)
            entries.append(entry)
            continue
        for r, row in enumerate(rows):
            ci_p = np.asarray(data.ci[r]).ravel()
            entry["roots"].append({
                "root": r,
                "e_ref": row.e_ref,
                "e_corr": row.e_corr,
                "e_total": row.e_total,
                "e_corr_eV": row.e_corr * EH_EV,
                "e_total_eV": row.e_total * EH_EV,
                "l2": expect("l2", ci_p),
                "lz": expect("lz", ci_p),
                "s2": expect("s2", ci_p),
                "classes": [{
                    "label": c.label, "note": c.note,
                    "norm": c.norm, "energy": c.energy,
                    "energy_meh": c.energy * 1e3,
                    "e_per_norm": c.e_per_norm,
                } for c in row.classes],
            })
        print(f"  [rows] {addr}: nroots {nroots}, "
              f"E(2) root0 {rows[0].e_corr * 1e3:+.4f} mEh "
              f"({rows[0].e_corr * EH_EV:+.4f} eV), "
              f"<L2> {entry['roots'][0]['l2']:.6f} "
              f"({time.time() - t0:.1f} s)", flush=True)
        entries.append(entry)
    return {
        "atom": atom,
        "provenance": {
            "engine": "nevpt2-sc (pyscf.mrpt.NEVPT, strongly contracted)",
            "pyscf_version": pyscf_version,
            "collected": datetime.now().isoformat(timespec="seconds"),
            "date": date.today().isoformat(),
            "declared_roots_policy": "every stored root is a row; no "
                                     "silent draws",
            "e_ref_partition": "kernel totals only (the readout's own "
                               "root totals); e_total = e_ref + e_corr",
            "wobble_note": "SC-NEVPT2 wobbles ~3e-7 Eh run-to-run "
                           "(thread-parallel reductions); the committed "
                           "payload is one draw, the store's chk is the "
                           "reproducibility unit",
            "stored_roots_bands": "energies 1e-11 absolute; CI "
                                  "amplitudes 1e-8 per-root; the "
                                  "degenerate block's joint projector "
                                  "(p_tol 1e-4, the solver-tail "
                                  "scale); the re-solve angle "
                                  "(davidson noise) 1e-3; the "
                                  "same-state energy noise 1e-6",
            "meters": "L2/Lz/S2 measured live through the CI-order "
                      "bridge on the stored draw",
            "lane_note": "the lane rows are the D2h-pinned components "
                         "(the blocks' own roots; the meters carry "
                         "their state class) -- the route-invariance "
                         "anchor is the lane's ground block against "
                         "the free CI twin's ground root",
        },
        "cells": entries,
    }


def main() -> None:
    atoms = sys.argv[1:]
    if not atoms:
        sys.exit("usage: run_pt2_ladder.py <atom> [<atom> ...]")
    for atom in atoms:
        print(f"===== pt2 ladder: {atom} =====", flush=True)
        payload = collect_atom(atom)
        out = Path("elements") / atom / "perturbation.json"
        out.write_text(json.dumps(payload, indent=1) + "\n")
        n_rows = sum(c["disposition"] == "rows" for c in payload["cells"])
        n_ref = len(payload["cells"]) - n_rows
        print(f"wrote {out}: {n_rows} rows / {n_ref} refused "
              f"of {len(payload['cells'])} declared cells", flush=True)


if __name__ == "__main__":
    main()
