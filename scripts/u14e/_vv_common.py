"""The u14e V&V family's shared core: metadata, the trusted reference,
the committed anchors, the analysis meters, the verdict ledger, the
row store and the one rendered table.  The DERIVATIONS live in the
legs -- each leg is readable top-to-bottom, PySCF-direct: the SCF
loaded from the trusted rung, the window written out, the CAS/SA
kernels called directly, the engine's own NEVPT per root.

This module owns only what the legs share (the metadata-from-a-common-
place allowance): the atom identity (read from the trusted rung and
cross-checked against the campaign declaration), the committed payload
anchors, the analysis meters (atomterm's angular operators -- the
project's double-certified construction; re-implementing L^2 per leg
would be its own error surface), the ledger and the table.

Run the legs from the repo root:
    python scripts/u14e/u14e_vv_4_12.py [--pinned] [--ground-only]
    python scripts/u14e/u14e_vv_4_15.py [--ground-only]
    python scripts/u14e/u14e_vv_6_13.py
    python scripts/u14e/u14e_vv_uhf.py
or all of them through the driver:
    python scripts/u14e/u14e_vv_orchestration.py
"""

import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import numpy as np  # noqa: E402
from pyscf import scf  # noqa: E402

from atomterm.cases import campaign_spec  # noqa: E402
from atomterm.detbasis import determinant_basis, operator_apply  # noqa: E402
from atomterm.operators import (  # noqa: E402
    ci_in_detbasis_order, l2_operator_parts, lz_operator,
    s2_operator_parts)

EH_EV = 27.211386245988
L_L2 = 72.0  # the 5L parent term: L = 8 -> L(L+1)
S_S2 = 6.0   # spin 4 -> S = 2 -> S(S+1)

ATOM_DIR = Path("elements") / "U14e"
ROWS = ATOM_DIR / "vv.json"
TABLE = ATOM_DIR / "VV.md"
PAYLOAD = ATOM_DIR / "perturbation.json"

# the assertion bands (declared constants; each cites the first green
# run's observed values -- logs/u14e_vv.log, 2026-11-20):
BAND_L2 = 1e-5     # the Casimir meters are gauge-invariant; observed:
                   # free 8.9e-15, pinned 1.9e-07 (17 roots each)
BAND_E2_MEH = 2.0  # a fresh SA draw is a new orbital gauge and E(2)
                   # is state-specific: the gauge-movement scale is
                   # the route record's Ti disclosure (1.76 mEh); the
                   # first green run's observed delta was 0.0000 mEh
BAND_ROUTE = 1e-6  # Eh; the comparison surface is the ENSEMBLE
                   # energy (variational, gauge-invariant across the
                   # routes), not E(2) -- E(2) moves at the mEh scale
                   # across gauges and is printed, never asserted;
                   # observed 5.3e-09 Eh


def load_rung(chk: str) -> tuple:
    """One rung, LOADED (never re-run): the molecule from the chk, an
    SCF object rebuilt around it, the converged record spliced in
    (``converged`` set True -- the trusted cache IS the converged
    record).  Detects the spin class from the record's mo_coeff."""
    from pyscf.lib import chkfile as _chk
    mol = _chk.load_mol(chk)
    rec = _chk.load(chk, "scf")
    if not isinstance(rec, dict) or "mo_coeff" not in rec:
        sys.exit(f"{chk}: no converged SCF record")
    uhf = np.asarray(rec["mo_coeff"]).ndim == 3
    mf = scf.UHF(mol) if uhf else scf.ROHF(mol)
    mf.__dict__.update(rec)
    mf.converged = True
    return mol, mf


def load_reference(chk: str | None = None):
    """The trusted rung, LOADED -- never re-run (the reproducibility
    unit).  The molecule comes from the chk itself (the campaign's own
    basis/ECP/spin); one cross-check asserts the chk's identity
    scalars still match a freshly built declared molecule, so the two
    cannot drift apart silently."""
    chk = chk or str(ATOM_DIR / "rungs" / "rohf.chk")
    mol, mf = load_rung(chk)
    from atomterm.system import build_mol
    declared = build_mol(campaign_spec("U14e").atom_config())
    if (mol.nelectron != declared.nelectron or mol.spin != declared.spin
            or mol.nao != declared.nao):
        sys.exit(f"{chk}: the rung's identity (nelectron "
                 f"{mol.nelectron}, spin {mol.spin}, nao {mol.nao}) "
                 f"drifted from the declaration (nelectron "
                 f"{declared.nelectron}, spin {declared.spin}, nao "
                 f"{declared.nao})")
    print(f"  ROHF (trusted rung, loaded): E = {mf.e_tot:.8f}, "
          f"nao {mol.nao}, nelectron {mol.nelectron}, spin "
          f"{mol.spin}", flush=True)
    return mol, mf


def anchor(addr: str) -> dict:
    """The committed payload row an assertion reads -- fail fast
    unless it prices roots (the declared-roots discipline)."""
    payload = json.loads(PAYLOAD.read_text())
    for c in payload["cells"]:
        if c["address"] == addr:
            if c.get("disposition") != "rows" or not c.get("roots"):
                sys.exit(f"anchor {addr} carries no priced roots")
            return c
    sys.exit(f"anchor {addr} not in {PAYLOAD}")


def meters(mol, c_act, na: int, nb: int):
    """The analysis meters as their ACTION (atomterm's angular
    operators through the matrix-free application; the many-body
    meters read per root through the CI-order bridge).  This is the
    one piece of atomterm analysis the legs keep -- the construction
    is double-certified in angular.py, and the meters are the terms'
    identity, not part of the derivation under test."""
    norb = c_act.shape[1]
    dets = determinant_basis(norb, na, nb)
    apply_l2 = operator_apply(dets, *l2_operator_parts(mol, c_act))
    apply_s2 = operator_apply(dets, *s2_operator_parts(norb))
    apply_lz = operator_apply(dets, one_body=lz_operator(mol, c_act))
    apply = {"l2": apply_l2, "s2": apply_s2, "lz": apply_lz}

    def expect(kind: str, ci_pyscf) -> float:
        w = ci_in_detbasis_order(np.asarray(ci_pyscf).ravel(),
                                 norb, na, nb)
        return float(np.real(w @ apply[kind](w)))

    return expect


# the verdict ledger --------------------------------------------------

class Verdict:
    """Every check lands one row; the exit code is the verdict.
    ``note`` rows are measured and recorded, never gating (the
    new-physics class)."""

    def __init__(self):
        self.rows = []

    def check(self, name, ok, detail):
        ok = bool(ok)
        self.rows.append((name, ok, "check"))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}",
              flush=True)
        return ok

    def note(self, name, detail):
        self.rows.append((name, True, "note"))
        print(f"  [note] {name}: {detail}", flush=True)

    def report(self) -> bool:
        print("\n===== verdict =====", flush=True)
        for name, ok, kind in self.rows:
            tag = ("PASS" if ok else "FAIL") if kind == "check" \
                else "REPORT"
            print(f"  {tag}  {name}", flush=True)
        checks = [(n, ok) for n, ok, k in self.rows if k == "check"]
        good = sum(ok for _, ok in checks)
        print(f"{good}/{len(checks)} assertions hold", flush=True)
        return good == len(checks)


# the row store + the one table ---------------------------------------

def drop_row(row: dict) -> None:
    """One leg's result row, merged into the family's store
    (``leg`` + ``route`` is the key)."""
    rows = json.loads(ROWS.read_text()) if ROWS.exists() else {}
    rows[f"{row['leg']}|{row['route']}"] = row
    ROWS.write_text(json.dumps(rows, indent=1) + "\n")
    render_table()


_ROW_ORDER = (("4_12", "free"), ("4_12", "pinned"), ("4_15", "free"),
              ("6_13", "free"), ("ref", "uhf"))


def render_table() -> None:
    """The small committed table (elements/U14e/VV.md): one row per
    run leg, beside its committed-payload anchor."""
    rows = json.loads(ROWS.read_text()) if ROWS.exists() else {}
    order = [k for k in _ROW_ORDER if k in rows]
    order += [k for k in rows if k not in order]
    out = [
        "# U14e -- the V&V family",
        "",
        "One row per leg of `scripts/u14e/u14e_vv_*.py` (the "
        "standalone, PySCF-direct derivations), beside the committed "
        "payload anchor (`perturbation.json`) it asserts against. "
        "Re-rendered by every leg after it drops its row; the full "
        "record lives in `logs/u14e_vv_*.log`.",
        "",
        "| leg | route | CAS | E_ref (Eh) | E(2) (mEh) | total (Eh) "
        "| anchor E(2) (mEh) | delta (mEh) | &lt;L²&gt; | &lt;S²&gt; "
        "| date |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for k in order:
        r = rows[k]
        a = ("—" if r.get("anchor_e2_meh") is None
             else f"{r['anchor_e2_meh']:+.4f}")
        out.append(
            f"| {r['leg']} | {r['route']} | {r['cas']} "
            f"| {r['e_ref']:.8f} | {r['e2_meh']:+.4f} "
            f"| {r['e_total']:.8f} | {a} "
            f"| {r['delta_meh']:+.4f} | {r['l2']} | {r['s2']} "
            f"| {r['date']} |")
    out += ["", "Atomic units; E(2) = SC-NEVPT2 (pyscf.mrpt, strongly "
            "contracted); the meters are atomterm's exact angular "
            "operators.", ""]
    TABLE.write_text("\n".join(out))
