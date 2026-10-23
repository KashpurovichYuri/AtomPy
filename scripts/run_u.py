"""U campaign (u_14e): Stuttgart RLC ECP78 + cut-Dyall basis.

Serial stages, each inspectable: gated Newton reference (cached) ->
window selection BY MEASUREMENT (open-MO Mulliken weights per primitive,
printed at start; exponent bands are the CONTIGUOUS weight sets around
the measured peaks, never hand-picked intervals) -> CASCI
per-configuration term tables (term x parity x <N_f>/<N_d> x dominant
configuration x E, per-block max|dE| vs PySCF) -> closure + commutator
fingerprints. The pur core substrate is mandatory on U (raw [H, L^2] =
2.2; pur 5.1e-13 -- the fact-17 defect at U scale).

Window bands (measured open-MO weights; f: 1.44->1.12, 0.674->0.91,
2.92->0.53, 0.294->0.41, 0.113->0.06; d: 0.172->0.41, 0.372->0.30,
0.076->0.22):
    1f    f (1.0, 2.0)                 7 + 5 AOs
    2f    f (0.5, 2.0)                14 + 5
    3f    f (0.2, 2.0)                21 + 5   (the converged window)
    4f    f (0.2, 4.0)                28 + 5   (tight-side stability)
    2f2d  f (0.5, 2.0), d (0.1, 0.5)  14 + 10  (d-side bias)
    3f2d  f (0.2, 2.0), d (0.1, 0.5)  21 + 10
Conclusions require neighboring-band stability (one more / one fewer
primitive per side -> same ordering).

Sectors >= LEAN_MIN determinants run on the TRUNCATED-roots path:
lowest eigenpairs only (scipy subset eigh), every expectation contracted
against those columns -- analyze()'s full-spectrum eigh + eigvals is the
O(n^3) trap at these sizes (measured: 370 s eigh at 14,950 dets; hours
at 40,920). Small sectors go through analyze() (the library label
contract) verbatim.

7s^2 is FROZEN in the core (ncore = 5: 6s 6p x3 7s): 7s^2 is 1S-closed,
so the ground term of 5f^3 6d^1 7s^2 IS the term of 5f^3 6d^1, and the
f^4 manifold lives in the same window.

Usage: python scripts/run_u.py [tier ...]      (default: all tiers)
Threads: OPENBLAS_NUM_THREADS=128 OMP_NUM_THREADS=128 (255 livelocks).
"""

import math
import os
import sys
from pathlib import Path

sys.path.insert(0, "src")
import numpy as np
from pyscf import gto
from scipy.sparse.linalg import eigsh

from atomterm.system import AtomConfig
from atomterm import drivers
from atomterm.analysis import analyze, term_symbol
from atomterm.detbasis import commutator, frobenius_norm, slater_condon_matrix
from atomterm.drivers import ReferenceSpec
from atomterm.operators import (
    build_operator_parts,
    build_operators,
    parity_operator,
    population_operator,
    window_closure,
)


def low_roots(h, nroots):
    """Lowest ``nroots`` eigenpairs by Lanczos (ARPACK eigsh).

    The measured wall: dense LAPACK eig phases are memory-bound at
    ~40-80 GFLOPs effective regardless of threads (370 s at 14,950 dets,
    part 1; >45 min without completing at 18,411, part 2) -- a subset
    dsyevr still pays the full tridiagonalization. ARPACK pays only
    matvecs (one gemm each). ``k = nroots + 16`` resolves the exactly
    degenerate multiplets (13-17 members) comfortably; ``v0`` fixed for
    determinism. Validated against the both-sides-exact 2f2d spectrum
    (7e-15 cross-check)."""
    n = h.shape[0]
    k = min(nroots + 16, n - 1)
    v0 = np.ones(n) / np.sqrt(n)
    evals, evecs = eigsh(h, k=k, which="SA", v0=v0,
                         ncv=min(n, max(6 * nroots, 256)), tol=0)
    order = np.argsort(evals)[:nroots]
    return evals[order], evecs[:, order]

HERE = Path(__file__).resolve().parent
ECP = (HERE / ".." / "data" / "u_stuttgart_rlc_ecp78.txt").read_text()
BAS = {"U": gto.parse((HERE / ".." / "data" / "u_dyall_cut_basis.txt")
                      .read_text())}

CFG = AtomConfig(symbol="U", basis=BAS, ecp=ECP, spin=4,
                 tag="U_stuttgart78_dyall")
REF = ReferenceSpec(occ_closed=("l0", "l1"), occ_open=("l2", "l3"),
                    occ_e_ref=-50.38185233, conv_accelerator="newton")
OUT = Path("results/u_14e")
OUT.mkdir(parents=True, exist_ok=True)
CHK = str(OUT / "U_stuttgart78_dyall_rohf.chk")
NELEC_ACT, NCORE = 4, 5  # 5f^3 6d^1 active; 6s 6p^3 7s core
LEAN_MIN = 6000           # dets; above: truncated-roots path
EXACT_MAX = 12000         # dets; PySCF side exact (pspace) up to here
# (measured: 10,626-det sector exact in ~6 min; 18,411 spent >45 min in
# dense LAPACK eig phases -- above this, PySCF runs Davidson at 1e-12
# and the dE column is labeled solver-class for those tiers)

CFG_LABEL = {0: "d4", 1: "f1d3", 2: "f2d2", 3: "f3d1", 4: "f4"}

TIERS = [  # science-critical first: the stability question, then the rest
    ("3f",   {"f": (0.2, 2.0), "d": (0.1, 0.25)}, [(4, 0, 48)]),
    ("4f",   {"f": (0.2, 4.0), "d": (0.1, 0.25)}, [(4, 0, 48)]),
    ("2f2d", {"f": (0.5, 2.0), "d": (0.1, 0.5)},  [(4, 0, 48)]),
    ("3f2d", {"f": (0.2, 2.0), "d": (0.1, 0.5)},  [(4, 0, 48)]),
    ("4f2d", {"f": (0.2, 4.0), "d": (0.1, 0.5)},  [(4, 0, 48)]),
    ("1f",   {"f": (1.0, 2.0), "d": (0.1, 0.25)},
     [(4, 0, 30), (3, 1, 40), (2, 2, 40)]),
    ("2f",   {"f": (0.5, 2.0), "d": (0.1, 0.25)},
     [(4, 0, 48), (3, 1, 40), (2, 2, 40)]),
    ("1f3d", {"f": (1.0, 2.0), "d": (0.05, 0.5)}, [(4, 0, 30)]),
]

# set in main(): the OCCUPIED closed block (canonical order is
# meaningless on near-degenerate uncontracted bases) and the pur core
CLOSED_IDX = None
C_PUR = None


def prim_weights(mol, mf, l):
    """Open-MO Mulliken S-weight per primitive of angular momentum l,
    summed over the open (occ == 1) reference MOs: the measured quantity
    the exponent bands are contiguous sets of."""
    chain = drivers.primitive_chain(mol, l)
    s = mol.intor_symmetric("int1e_ovlp")
    w = {}
    for exp, rows in chain:
        tot = 0.0
        for k in np.flatnonzero(mf.mo_occ == 1.0):
            c = mf.mo_coeff[:, k]
            tot += float(np.real(c[rows] @ (s @ c)[rows]))
        w[exp] = tot
    return w


def build_window(mol, bands):
    """(f_rows, d_rows, win): AO rows of the {l: exponent_band} window,
    f rows first -- the column layout the configuration tables index."""
    f_rows = drivers.select_primitives(mol, 3, exponent_band=bands["f"])
    d_rows = drivers.select_primitives(mol, 2, exponent_band=bands["d"])
    return f_rows, d_rows, f_rows + d_rows


def root_expectations(mol, c_act, f_rows, na, nb, evecs, nroots, s2, l2):
    """Per-root labels, truncated to the roots read (the O(n^3) trap is
    any contraction over the full eigenvector set). Returns dict with
    parity, <N_f>, <S^2>, <L^2> arrays and per-root configuration
    weights over integer N_f (window columns are l-pure, so every
    determinant has an EXACT integer f-occupation)."""
    from atomterm.detbasis import determinant_basis

    norb = c_act.shape[1]
    dets = determinant_basis(norb, na, nb)
    nf_cols = len(f_rows)
    nf_per_det = np.array(
        [sum(1 for so in det if so % norb < nf_cols) for det in dets])

    vec = evecs[:, :nroots]
    par = np.diag(parity_operator(mol, c_act, na, nb))
    nf_op = population_operator(mol, c_act, f_rows, na, nb)
    out = {
        "par": np.einsum("ik,i,ik->k", vec.conj(), par, vec).real,
        "nf": np.einsum("ik,ij,jk->k", vec.conj(), nf_op, vec).real,
        "s2": np.einsum("ik,ij,jk->k", vec.conj(), s2, vec).real,
        "l2": np.einsum("ik,ij,jk->k", vec.conj(), l2, vec).real,
        "w": [np.bincount(nf_per_det, weights=np.abs(vec[:, k]) ** 2,
                          minlength=NELEC_ACT + 1)
              for k in range(nroots)],
    }
    return out


def config_label(w):
    """Dominant-configuration string with runner-up (weight >= 0.10)."""
    nf = int(np.argmax(w))
    txt = f"{CFG_LABEL[nf]} ({w[nf]:.2f})"
    order = np.argsort(-w)
    for j in order:
        if j != nf and w[j] >= 0.10:
            txt += f" | {CFG_LABEL[j]} {w[j]:.2f}"
            break
    return txt


def render_table(tag, ndet, nroots, e_tot, exp, cas_energies,
                 comm_lz, comm_l2, h_norm):
    """The STANDARD report: per-root term x parity x <N_f>/<N_d> x
    dominant configuration x E; per-block max|dE| vs PySCF by
    sorted-position matching (both spectra exact: ours Lanczos, theirs
    pspace-full up to EXACT_MAX -- the residual measures the fold,
    machine class on the m-scalar pur core; above EXACT_MAX PySCF runs
    Davidson and the column is solver-class, noted in the header)."""
    lines = [f"[{tag}] ndet={ndet}  roots={nroots}"]
    lines.append(f"||[H,Lz]|| = {comm_lz:.2e} (rel {comm_lz / h_norm:.2e})"
                 f"  ||[H,L2]|| = {comm_l2:.2e} "
                 f"(rel {comm_l2 / h_norm:.2e})")
    d_e = float(np.max(np.abs(e_tot - cas_energies)))
    note = ("pur core: exactly m-scalar, the two fold conventions "
            "coincide" if d_e < 1e-8 else
            "PySCF Davidson class (exact only below EXACT_MAX)")
    lines.append(f"max|dE| vs PySCF = {d_e:.2e}  ({note})")
    lines.append(f"{'#':>3} {'term':>5} {'par':>4} {'<Nf>':>6} {'<Nd>':>6}"
                 f"  {'config (w)':<26} {'E (Eh)':>12} {'rel mEh':>9}")
    blocks: dict[str, list[float]] = {}
    for k in range(nroots):
        blk = CFG_LABEL[int(np.argmax(exp["w"][k]))]
        blocks.setdefault(blk, []).append(abs(e_tot[k] - cas_energies[k]))
        term = term_symbol(exp["s2"][k], exp["l2"][k])
        if exp["par"][k] < 0:
            term += "\u00b0"  # odd-parity marker (superscript circle)
        lines.append(
            f"{k:3d} {term:>5} {'odd' if exp['par'][k] < 0 else 'even':>4}"
            f" {exp['nf'][k]:6.3f} {NELEC_ACT - exp['nf'][k]:6.3f}"
            f"  {config_label(exp['w'][k]):<26} {e_tot[k]:12.6f}"
            f" {(e_tot[k] - e_tot[0]) * 1e3:9.2f}")
    lines.append("")
    lines.append(f"{'block':>6} {'n':>3} {'E_ground (Eh)':>13}"
                 f" {'max|dE|':>9}")
    grounds = {b: min(e_tot[k] for k in range(nroots)
                      if CFG_LABEL[int(np.argmax(exp["w"][k]))] == b)
               for b in blocks}
    for blk in sorted(grounds, key=lambda b: grounds[b]):
        lines.append(f"{blk:>6} {len(blocks[blk]):3d} {grounds[blk]:13.6f}"
                     f" {max(blocks[blk]):9.2e}")
    deg = int(np.sum(np.abs(e_tot - e_tot[0]) < 1e-6))
    lines.append(f"ground multiplet: x{deg} within 1e-6 Eh")
    return "\n".join(lines)


def run_sector(mol, mf, tier, bands, na, nb, nroots):
    """One sector of one tier: fold on the pur substrate, CASCI vs PySCF,
    term table via analyze() (small) or the truncated-roots path (big)."""
    f_rows, d_rows, win = build_window(mol, bands)
    ncas, tag = len(win), f"tier{tier}_ms{na - nb}"
    ndet = math.comb(ncas, na) * math.comb(ncas, nb)
    nroots = min(nroots, ndet)

    mo = drivers.build_full_mo(mf, CLOSED_IDX, win, c_core=C_PUR)
    # analyze the ORTHONORMALIZED window (build_full_mo Loewdin-rotates
    # overlapping same-l primitives): span- and closure-preserving, and
    # the basis the CASCI actually solved in; f/d blocks stay separate
    # (exactly S-orthogonal between different l)
    c_act = mo[:, NCORE:NCORE + ncas]
    # PySCF side EXACT (pspace = full space) where affordable: Davidson
    # straggles inside the degenerate multiplets (measured: 38 mEh in
    # the 2f2d clusters); above EXACT_MAX only a tight tolerance
    fci_opts = ({"pspace_size": ndet} if ndet <= EXACT_MAX
                else {"conv_tol": 1e-12})
    cas = drivers.run_casci(mf, ncas, (na, nb), mo, CHK, nroots=nroots,
                            name=tag, fci_opts=fci_opts)

    if ndet <= LEAN_MIN:
        ops = build_operators(mol, mf, c_act, NCORE, na, nb,
                              c_core=mo[:, :NCORE])
        res = analyze(ops["H"], ops["Lz"], ops["L2"], ops["S2"], tag)
        ev, vec = np.linalg.eigh(ops["H"].real)
        exp = root_expectations(mol, c_act, f_rows, na, nb, vec, nroots,
                                ops["S2"], ops["L2"])
        comm = (res.comm_norm_h_lz, res.comm_norm_h_l2, res.h_norm)
    else:
        parts = build_operator_parts(mol, mf, c_act, NCORE, na, nb,
                                     c_core=mo[:, :NCORE])
        dets = parts.pop("dets")
        parts.pop("h1")
        parts.pop("g")
        mats = {k: slater_condon_matrix(dets, one_body=ob, two_body=tb)
                for k, (ob, tb) in parts.items() if k != "N"}
        h = mats.pop("H")
        ev, vec = low_roots(h, nroots)
        comm_lz = frobenius_norm(commutator(h, mats["Lz"]))
        comm_l2 = frobenius_norm(commutator(h, mats["L2"]))
        comm = (comm_lz, comm_l2, max(frobenius_norm(h), 1e-300))
        exp = root_expectations(mol, c_act, f_rows, na, nb, vec, nroots,
                                mats["S2"], mats["L2"])
        del mats, h

    e_tot = ev[:nroots] + cas.core_energy
    print(render_table(tag, ndet, nroots, e_tot, exp, cas.energies, *comm))
    return e_tot, exp


def main():
    wanted = sys.argv[1:]
    mol, mf = drivers.run_rohf(CFG, CHK, ref=REF)
    print(f"reference: E = {mf.e_tot:.8f} Eh  (5f^3 6d^1 7s^2, "
          f"Newton-gated)  threads: "
          f"{os.environ.get('OPENBLAS_NUM_THREADS', 'default')}")
    global CLOSED_IDX, C_PUR
    closed_idx = list(np.flatnonzero(mf.mo_occ == 2.0))
    assert len(closed_idx) == NCORE, (len(closed_idx), NCORE)
    CLOSED_IDX = closed_idx
    # pur substrate: the canonical core's foreign-l tails break window
    # closure at the 1e-2 level on U (measured); pur closes at 5e-13
    C_PUR, contam = drivers.purified_core_orbitals(mol, mf, NCORE)
    print(f"pur core: max contamination {max(contam):.2e}")

    # the measurement the bands are contiguous sets of (never asserted
    # from labels: uncontracted n-numbers are fake, fact-1 amendment)
    for l, name in [(3, "f"), (2, "d")]:
        w = prim_weights(mol, mf, l)
        print(f"{name} chain open-MO weights:",
              " ".join(f"{e:.4f}->{v:.3f}" for e, v in w.items()))

    for tier, bands, sectors in TIERS:
        if wanted and tier not in wanted:
            continue
        _, _, win = build_window(mol, bands)
        cl = window_closure(mol, np.eye(mol.nao_nr())[:, win])
        print(f"\n{'=' * 68}\ntier {tier}: {len(win)} AOs  raw-window "
              f"closure lz={cl['lz']:.1e} lperp={cl['lperp']:.1e}\n"
              f"{'=' * 68}")
        for na, nb, nroots in sectors:
            run_sector(mol, mf, tier, bands, na, nb, nroots)


if __name__ == "__main__":
    main()
