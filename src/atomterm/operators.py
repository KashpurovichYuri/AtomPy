"""Assembly of physical operators as full many-electron matrices.

Bridges PySCF (AO/MO integrals, angular.py) to detbasis (Slater-Condon rules).
Everything is plain first quantization: an operator enters the determinant
space as (one-body matrix, two-body integral tensor) and leaves as an
``(ndet, ndet)`` matrix.

Spin-orbital conventions follow atomterm.detbasis: alpha block first, beta
block second, two-body integrals in chemist notation over spin-orbitals.

Operators provided
------------------
* Hamiltonian in the active space (core folded into the one-body part,
  the same convention as PySCF's mcscf ``h1eff``; eigenvalues are directly
  comparable to ``mc.e_cas``).
* Lz (one-body); L^2 = Sum_i l^2(i) + 2 Sum_{i<j} l_j.l_i
  (one-body ``Sum_a l_a.l_a`` + two-body ``2 Sum_a l_a (x) l_a``);
* S^2 = Sum_i 3/4 + 2 Sum_{i<j} s_j.s_i;
* number operator N.

The matrix products ``l_a.l_a`` are exact when the one-particle space is
rotation-closed (true for a shell of atomic AOs, and true up to the orbital
optimization residual for CASSCF MOs -- which is exactly what the
``[H, L]`` commutator norms measure).
"""

from __future__ import annotations

import numpy as np

from atomterm.detbasis import determinant_basis, slater_condon_matrix
from atomterm.angular import l_matrices

__all__ = [
    "spin_block_one_body",
    "spin_block_two_body",
    "active_space_integrals",
    "full_space_integrals",
    "lz_operator",
    "l2_operator_parts",
    "s2_operator_parts",
    "number_operator",
    "build_operator_parts",
    "build_operators",
]


def spin_block_one_body(a_spatial: np.ndarray) -> np.ndarray:
    """Spatial one-body matrix -> spin-orbital (alpha | beta) block-diagonal."""
    z = np.zeros_like(a_spatial)
    return np.block([[a_spatial, z], [z, a_spatial]]).astype(complex)


def spin_block_two_body(g_spatial: np.ndarray) -> np.ndarray:
    """Spatial chemist integrals (pq|rs) -> spin-orbital chemist integrals.

    Spin of electron 1 must match between bra and ket (p with r) and likewise
    for electron 2 (q with s): blocks (aa;aa), (ab;ab), (ba;ba), (bb;bb).
    """
    n = g_spatial.shape[0]
    gso = np.zeros((2 * n, 2 * n, 2 * n, 2 * n), dtype=g_spatial.dtype)
    ia, ib = slice(0, n), slice(n, 2 * n)
    gso[ia, ia, ia, ia] = g_spatial
    gso[ia, ib, ia, ib] = g_spatial
    gso[ib, ia, ib, ia] = g_spatial
    gso[ib, ib, ib, ib] = g_spatial
    return gso


def _mo_transform(a_ao: np.ndarray, c_act: np.ndarray) -> np.ndarray:
    """C^dag A C."""
    return c_act.conj().T @ a_ao @ c_act


def _spherical_average_density(mol, p: np.ndarray) -> np.ndarray:
    """l=0 (totally symmetric) component of an AO density matrix.

    Within every pair of same-l shells the density block is replaced by its
    average over the shared m index: P[(c1,m),(c2,m')] -> delta_mm' * mean_m
    P[(c1,m),(c2,m)]. Cross-l blocks vanish by exact symmetry and are zeroed.

    For an atom the core density is spherical by construction; this removes
    nothing but numerical noise from the SCF solution and restores the exact
    rotational invariance of the analysis Hamiltonian.
    """
    ao_loc = mol.ao_loc_nr()
    shells = [
        (mol.bas_angular(k), mol.bas_nctr(k), 2 * mol.bas_angular(k) + 1, ao_loc[k])
        for k in range(mol.nbas)
    ]
    ps = np.zeros_like(p)
    for l1, nc1, nk1, o1 in shells:
        for l2, nc2, nk2, o2 in shells:
            if l1 != l2:
                continue  # zero by exact symmetry
            for c1 in range(nc1):
                for c2 in range(nc2):
                    blk = np.zeros((nk1, nk2), dtype=p.dtype)
                    for m in range(nk1):
                        blk[m, m] = np.mean(
                            [p[o1 + c1 * nk1 + j, o2 + c2 * nk2 + j] for j in range(nk1)]
                        )
                    ps[o1 + c1 * nk1 : o1 + (c1 + 1) * nk1,
                       o2 + c2 * nk2 : o2 + (c2 + 1) * nk2] = blk
    return ps


def active_space_integrals(mol, mf, c_act: np.ndarray, ncore: int):
    """(h1, g) of the active space, core contributions folded into h1.

    mcscf convention: ``h1 = C^dag (h_core + V_eff(core)) C`` with the
    spherically-averaged closed-shell core density (the l=0 component: the
    physical density of closed s shells in a central field) and
    ``V_eff[m,n] = sum_ls P[l,s] (mn|ls) - 0.5 (ml|ns)``.

    ``g`` follows the atomterm.detbasis convention (physicist arrangement):
    ``g[p,q,r,s] = <pq|rs> = int p(1) q(2) r(1) s(2)`` = chemist ``(pr|qs)``.
    CI eigenvalues of (h1, g) equal ``mc.e_cas``.
    """
    c_core = mf.mo_coeff[:, :ncore]
    d_core = _spherical_average_density(mol, 2.0 * c_core @ c_core.conj().T)

    eri = mol.intor("int2e", aosym="s1")  # chemist (mu nu | la si)
    veff_core = (
        np.einsum("ls,mnls->mn", d_core, eri, optimize=True)
        - 0.5 * np.einsum("ls,mlns->mn", d_core, eri, optimize=True)
    )
    h1_ao = mf.get_hcore() + veff_core
    h1 = _mo_transform(h1_ao, c_act)

    # detbasis convention: g[p,q,r,s] = <pq|rs> = int p(1) q(2) r(1) s(2)
    # (physicist arrangement) = chemist (pr|qs), hence the transpose.
    g = np.einsum(
        "up,vq,wr,xs,uvwx->pqrs",
        c_act.conj(),
        c_act.conj(),
        c_act,
        c_act,
        eri,
        optimize=True,
    ).transpose(0, 2, 1, 3)
    return h1, g


def _composite_two_body(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Chemist integrals of the two-electron operator A(1)B(2):
    ``(pq|rs) = A[p, r] B[q, s]``."""
    return np.einsum("pr,qs->pqrs", a, b)


def lz_operator(mol, c_act: np.ndarray) -> np.ndarray:
    """Lz as a spin-orbital one-body matrix (spatially diagonal)."""
    return spin_block_one_body(_mo_transform(l_matrices(mol)["z"], c_act))


def l2_operator_parts(mol, c_act: np.ndarray):
    """(one_body, two_body) parts of the total L^2."""
    n = c_act.shape[1]
    one = np.zeros((n, n), dtype=complex)
    two = np.zeros((n, n, n, n), dtype=complex)
    for a in "xyz":
        la = _mo_transform(l_matrices(mol)[a], c_act)
        one += la @ la
        two += 2.0 * _composite_two_body(la, la)
    return spin_block_one_body(one), spin_block_two_body(two)


def s2_operator_parts(norb: int):
    """(one_body, two_body) parts of the total S^2.

    S^2 = Sz^2 + (S+ S- + S- S+)/2 built from first-quantized composites
    ``Sum_{i<j} A(i) B(j)``; equivalently ``Sum_i 3/4 + 2 Sum_{i<j} s_j.s_i``.
    """
    d = 2 * norb
    sz_diag = np.zeros(d)
    sz_diag[:norb], sz_diag[norb:] = 0.5, -0.5

    s_plus = np.zeros((d, d), dtype=complex)
    rows = np.arange(norb)
    s_plus[rows, rows + norb] = 1.0  # <p alpha| s+ |p beta> = 1
    s_minus = s_plus.conj().T

    one = 0.75 * np.eye(d, dtype=complex)
    two = 2.0 * (
        _composite_two_body(np.diag(sz_diag), np.diag(sz_diag))
        + 0.5 * _composite_two_body(s_plus, s_minus)
        + 0.5 * _composite_two_body(s_minus, s_plus)
    )
    return one, two


def number_operator(norb: int) -> np.ndarray:
    """Electron number N (one-body identity in spin-orbital space)."""
    return spin_block_one_body(np.eye(norb))


def full_space_integrals(mol, mf, w: np.ndarray):
    """(h1, g) over the given orbital columns with NO core folding.

    Raw one-electron h_core and raw two-electron ERIs (detbasis physicist
    convention) in the orbital basis ``w`` -- for core-explicit (all-electron
    active space) second-quantized Hamiltonians.
    """
    h1 = _mo_transform(mf.get_hcore(), w)
    eri = mol.intor("int2e", aosym="s1")
    g = np.einsum(
        "up,vq,wr,xs,uvwx->pqrs",
        w.conj(),
        w.conj(),
        w,
        w,
        eri,
        optimize=True,
    ).transpose(0, 2, 1, 3)
    return h1, g


def build_operator_parts(mol, mf, c_act: np.ndarray, ncore: int, na: int, nb: int) -> dict:
    """Operator content of build_operators as raw parts, evaluation-agnostic:
    ``{name: (one_body, two_body)}`` plus the determinant basis and integrals.
    """
    h1, g = active_space_integrals(mol, mf, c_act, ncore)
    dets = determinant_basis(c_act.shape[1], na, nb)
    s2_one, s2_two = s2_operator_parts(c_act.shape[1])
    l2_one, l2_two = l2_operator_parts(mol, c_act)
    return {
        "dets": dets,
        "h1": h1,
        "g": g,
        "H": (spin_block_one_body(h1), spin_block_two_body(g)),
        "Lz": (lz_operator(mol, c_act), None),
        "L2": (l2_one, l2_two),
        "S2": (s2_one, s2_two),
        "N": (number_operator(c_act.shape[1]), None),
    }


def build_operators(
    mol, mf, c_act: np.ndarray, ncore: int, na: int, nb: int, method: str = "slater_condon"
) -> dict:
    """All analysis operators as full matrices in the determinant basis.

    method: "slater_condon" (first quantization, explicit excitation rules)
    or "secondquant" (fermionic a+/a strings applied to occupation vectors).
    Both consume identical operator parts and produce elementwise identical
    matrices (pinned by tests).
    """
    parts = build_operator_parts(mol, mf, c_act, ncore, na, nb)
    dets = parts.pop("dets")
    integrals = (parts.pop("h1"), parts.pop("g"))

    if method == "slater_condon":
        from atomterm.detbasis import slater_condon_matrix

        def build(ob, tb):
            return slater_condon_matrix(dets, one_body=ob, two_body=tb)
    elif method == "secondquant":
        from atomterm.secondquant import secondquant_matrix

        def build(ob, tb):
            return secondquant_matrix(dets, one_body=ob, two_body=tb)
    else:
        raise ValueError(f"unknown method: {method!r}")

    out = {"dets": dets, "h1": integrals[0], "g": integrals[1]}
    for name, (ob, tb) in parts.items():
        out[name] = build(ob, tb)
    return out
