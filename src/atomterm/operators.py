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

from atomterm.detbasis import determinant_basis, frobenius_norm, slater_condon_matrix
from atomterm.angular import l_matrices

__all__ = [
    "spin_block_one_body",
    "spin_block_two_body",
    "active_space_integrals",
    "full_space_integrals",
    "lz_operator",
    "parity_one_body",
    "parity_operator",
    "l2_operator_parts",
    "window_closure",
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


def active_space_integrals(
    mol, mf, c_act: np.ndarray, ncore: int, spherical_core: bool = True,
    c_core: np.ndarray = None,
):
    """(h1, g) of the active space, core contributions folded into h1.

    mcscf convention: ``h1 = C^dag (h_core + V_eff(core)) C`` with the
    spherically-averaged closed-shell core density (the l=0 component: the
    physical density of closed s shells in a central field) and
    ``V_eff[m,n] = sum_ls P[l,s] (mn|ls) - 0.5 (ml|ns)``.

    ``c_core`` overrides the folded core columns (default ``mf.mo_coeff[:,
    :ncore]``). Orbital-optimized stages (CASSCF/SA) must pass their OWN
    converged inactive columns: the optimizer may legally reorganize the
    core span (core-active swaps), and folding a different core than the
    one the energies were computed with breaks the identity
    ``E_tot = E_core + e_cas`` (observed as O(Eh) nonsense, not small drift).

    ``spherical_core=False`` skips the density averaging: this reproduces
    the raw SCF core field, which carries the reference's anisotropy (the
    core is polarized by ONE member of the degenerate open-shell manifold,
    e.g. the |px py> mean field). Formally ``H_eff = H0 + eps*W`` with
    H0 rotationally invariant and W the anisotropic core field of weight
    eps ~ 1e-5: commutators [H, Lz], [H, L^2] = O(eps), and the 3P
    multiplet splits as {M_L=+-1}/{0} (a k=2, m=0 field couples dM_L = 0,
    +-2; the dipole k=1 is excluded exactly by inversion symmetry). The
    averaged density is the ensemble average over the degenerate manifold
    -- the scalar projection -- so it restores exact invariance.

    ``g`` follows the atomterm.detbasis convention (physicist arrangement):
    ``g[p,q,r,s] = <pq|rs> = int p(1) q(2) r(1) s(2)`` = chemist ``(pr|qs)``.
    CI eigenvalues of (h1, g) equal ``mc.e_cas``.
    """
    c_core = mf.mo_coeff[:, :ncore] if c_core is None else c_core
    d_raw = 2.0 * c_core @ c_core.conj().T
    d_core = _spherical_average_density(mol, d_raw) if spherical_core else d_raw

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


def parity_one_body(mol, c_act: np.ndarray) -> np.ndarray:
    """Per-window-orbital parity signs ``(-1)^l`` (diagonal matrix).

    Requires RAW AO columns (each column one basis function of definite
    l) -- asserted: optimized windows mix l within a column and per-
    column parity is undefined there. Parity labels therefore live on
    the raw-window census stages.
    """
    from atomterm.drivers import lm_sectors

    sectors = lm_sectors(mol)
    row_l = np.zeros(mol.nao_nr(), dtype=int)
    for (l, m), rows in sectors.items():
        for r in rows:
            row_l[r] = l
    unit = np.eye(mol.nao_nr())
    signs = np.zeros(c_act.shape[1])
    for j in range(c_act.shape[1]):
        r = int(np.argmax(np.abs(c_act[:, j])))
        if np.abs(c_act[:, j] - unit[:, r]).max() > 1e-10:
            raise ValueError(
                "parity_one_body needs raw AO columns (definite l per "
                "column); optimized windows break per-column parity")
        signs[j] = (-1.0) ** row_l[r]
    return np.diag(signs)


def parity_operator(mol, c_act: np.ndarray, na: int, nb: int) -> np.ndarray:
    """Parity as a diagonal matrix in the determinant basis: each
    determinant's eigenvalue is the PRODUCT ``prod (-1)^{l_i}`` over its
    occupied window orbitals.

    Structurally outside the (one_body, two_body) parts protocol: that
    machinery SUMS one-body contributions (N, Lz, H1), while parity is
    multiplicative -- so this is a direct construction over the
    determinant occupations, not a Slater-Condon evaluation.
    """
    from atomterm.detbasis import determinant_basis

    signs = np.diag(parity_one_body(mol, c_act))
    norb = c_act.shape[1]
    dets = determinant_basis(norb, na, nb)
    out = np.ones(len(dets))
    for k, det in enumerate(dets):
        for so in det:
            out[k] *= signs[so % norb]
    return np.diag(out)


def window_closure(mol, c_win: np.ndarray) -> dict[str, float]:
    """One-particle rotational-closure defects of a window SPAN.

    The span of the S-orthonormal columns ``c_win`` is rotation-closed
    iff its projector commutes with all three rotation generators. In
    the symmetrically orthonormalized AO basis (fact 4 metric dance) the
    projector is ``P = Q Q^dag`` with ``Q = S^(1/2) c_win`` (QR-cleaned)
    and the generators act as ``S^(-1/2) l_a S^(-1/2)``. Reported as
    Frobenius norms of the commutators:

    ``lz``    = ||[P, Lz]||  -- exactly the part the many-body ``[H, Lz]``
    can see: Lz mixes only within +-m pairs, so m=0-member tilts and
    EQUAL tilts of a full +-m pair are invisible here (fact 18c blind
    spots);
    ``lperp`` = sqrt(||[P, Lx]||^2 + ||[P, Ly]||^2) -- the defect only
    L+- (hence L^2, trivial in maximal-Ms single-irrep sectors) sees.

    ``lz = lperp = 0`` (to roundoff) is the COMPLETE closure certificate
    of the one-particle space: no blind spots, no core fold, no CI --
    symmetry breaking is attributed to the span itself. The many-body
    commutator norms respond to these defects at O(defect) while the
    energy responds at O(defect^2) at stationary points (facts 16/18).
    """
    s = mol.intor("int1e_ovlp")
    ev, u = np.linalg.eigh(s)
    s_half = (u * np.sqrt(ev)) @ u.T
    s_ihalf = (u / np.sqrt(ev)) @ u.T
    q, _ = np.linalg.qr(s_half @ c_win)
    p = q @ q.conj().T
    defect = {}
    for axis, l_ao in l_matrices(mol).items():
        a = s_ihalf @ l_ao @ s_ihalf
        defect[axis] = frobenius_norm(p @ a - a @ p)
    return {
        "lz": float(defect["z"]),
        "lperp": float(np.hypot(defect["x"], defect["y"])),
    }


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


def core_anisotropy(mol, mf, ncore: int, c_core: np.ndarray = None) -> float:
    """Relative Frobenius norm of the discarded core-density anisotropy
    ``||rho - rho_sph||_F / ||rho_sph||_F`` -- the ``eps`` weight of the
    internal-field term in ``H_eff = H0 + eps*W`` (see
    :func:`active_space_integrals`)."""
    c_core = mf.mo_coeff[:, :ncore] if c_core is None else c_core
    d = 2.0 * c_core @ c_core.conj().T
    d_sph = _spherical_average_density(mol, d)
    return float(frobenius_norm(d - d_sph) / frobenius_norm(d_sph))


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


def build_operator_parts(
    mol, mf, c_act: np.ndarray, ncore: int, na: int, nb: int,
    spherical_core: bool = True, c_core: np.ndarray = None,
) -> dict:
    """Operator content of build_operators as raw parts, evaluation-agnostic:
    ``{name: (one_body, two_body)}`` plus the determinant basis and integrals.
    """
    h1, g = active_space_integrals(mol, mf, c_act, ncore, spherical_core, c_core)
    dets = determinant_basis(c_act.shape[1], na, nb)
    s2_one, s2_two = s2_operator_parts(c_act.shape[1])
    l2_one, l2_two = l2_operator_parts(mol, c_act)
    parts = {
        "dets": dets,
        "h1": h1,
        "g": g,
        "H": (spin_block_one_body(h1), spin_block_two_body(g)),
        "Lz": (lz_operator(mol, c_act), None),
        "L2": (l2_one, l2_two),
        "S2": (s2_one, s2_two),
        "N": (number_operator(c_act.shape[1]), None),
    }
    return parts


def build_operators(
    mol, mf, c_act: np.ndarray, ncore: int, na: int, nb: int,
    method: str = "slater_condon", spherical_core: bool = True,
    c_core: np.ndarray = None,
) -> dict:
    """All analysis operators as full matrices in the determinant basis.

    method: "slater_condon" (first quantization, explicit excitation rules)
    or "secondquant" (fermionic a+/a strings applied to occupation vectors).
    Both consume identical operator parts and produce elementwise identical
    matrices (pinned by tests). ``spherical_core``/``c_core``: see
    :func:`active_space_integrals`.
    """
    parts = build_operator_parts(mol, mf, c_act, ncore, na, nb,
                                 spherical_core, c_core)
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
