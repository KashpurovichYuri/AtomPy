"""One-electron angular momentum operators Lx, Ly, Lz, L^2 in the AO
basis -- the angular-momentum layer every term classification and
symmetry certificate rests on.

PURPOSE: exact one-electron L matrices for a single atom. GOAL: TWO
independent constructions that certify each other elementwise, so
the operator layer is never the weak link: the libcint integral path
(primary) and the exact tesseral algebra via overlaps (cross-check).

The primary construction (libcint, via PySCF)
---------------------------------------------
PySCF exposes the libcint integral ``int1e_irp`` with 9 components;
decoding (verified empirically and pinned by tests) gives

    T[3*a + b] = Integral mu | r_a d_b nu      (a, b = 0..2 for x, y, z)

The angular momentum operator ``L = -i r x grad`` is the antisymmetric
part:

    Lz = -i (x dy - y dx) = -i (T[1] - T[3])
    Lx = -i (y dz - z dy) = -i (T[5] - T[7])
    Ly = -i (z dx - x dz) = -i (T[6] - T[2])

Exact for any pure-spherical AO basis on a single center at the
origin.

The analytic cross-check
------------------------
For a single center at the origin every AO is a radial function
times a real spherical harmonic, so the generators act only within
each (l, slot) block (slots walked by ``windows.shell_slots`` -- fused
shells included). The real-to-complex harmonic transform ``U^(l)``
follows from exact tesseral algebra in the Condon-Shortley convention
(with PySCF's p-shell ordering),

    gen^(l) = U^dag gen_c^(l) U,   gen_c = generators in the |l,m> basis

``Lz|mu> = Sum_k blk[k, pos_mu] |chi_k>`` for k in the same slot, hence
``<mu|Lz|nu> = Sum_k blk[k, pos_nu] S[mu, k]``.  For ``L^2`` see
:func:`l2_matrix`.

Layering: `l_matrices` (libcint) is primary; `l_ladder_matrices` and
`l2_matrix` derive from it; `real_to_complex` +
`_complex_generators` + the shell walker build the independent
analytic path. Valid for a single atom at the origin with a spherical
basis (enforced here and by ``system.build_mol``).
"""

from __future__ import annotations

import numpy as np

from atomterm.windows import s_metric_roots, shell_slots

__all__ = [
    "l_matrices",
    "l_ladder_matrices",
    "l2_matrix",
    "l_matrices_analytic",
    "real_to_complex",
]


# ------------------------------------------------- the libcint construction

def _check_single_center(mol) -> None:
    """The L operators are defined about r = 0 on spherical functions
    of ONE center (atoms built through system.build_mol satisfy this;
    direct mol consumers are checked here)."""
    if mol.cart:
        raise ValueError("spherical (pure) basis required")
    if mol.natm != 1 or not np.allclose(mol.atom_coords()[0], 0.0, atol=1e-10):
        raise ValueError("operators are defined for a single atom at the origin")


def l_matrices(mol) -> dict[str, np.ndarray]:
    """AO matrices of Lx, Ly, Lz via the libcint ``int1e_irp``
    integrals (the decoding: the module docstring).

    Returns
    -------
    ``{"x", "y", "z"}``: complex Hermitian matrices, purely imaginary
    on a real AO basis.

    Raises
    ------
    ValueError
        For a cartesian basis, several atoms, or an off-origin atom.
    """
    _check_single_center(mol)
    t = mol.intor("int1e_irp", comp=9)
    img = -1j
    return {
        "x": img * (t[5] - t[7]),
        "y": img * (t[6] - t[2]),
        "z": img * (t[1] - t[3]),
    }


def l_ladder_matrices(mol) -> dict[str, np.ndarray]:
    """AO matrices of the orbital ladders ``L± = Lx ± i Ly`` (built
    from :func:`l_matrices`; the analytic path assembles its Lx/Ly
    FROM these ladders in the |l,m> basis -- :func:`_complex_generators`).
    """
    l = l_matrices(mol)
    return {"+": l["x"] + 1j * l["y"], "-": l["x"] - 1j * l["y"]}


# ------------------------------------------------------------------- L^2

def l2_matrix(mol) -> np.ndarray:
    """AO matrix of the one-electron ``L^2 = Lx.Lx + Ly.Ly + Lz.Lz``
    (exact).

    Notes
    -----
    Matrix products are exact only in an orthonormal basis, so the
    product is taken in the symmetrically orthonormalized AO basis
    and transported back (:func:`windows.s_metric_roots`):

    ``L2_ao = S^(1/2) [ sum_a (S^-1/2 l_a S^-1/2)^2 ] S^(1/2)``.
    """
    _check_single_center(mol)
    s = mol.intor("int1e_ovlp")
    s_half, s_ihalf = s_metric_roots(s)
    l_ops = l_matrices(mol)
    l2_orth = sum(
        (s_ihalf @ l @ s_ihalf) @ (s_ihalf @ l @ s_ihalf) for l in l_ops.values()
    )
    return s_half @ l2_orth @ s_half


# ------------------------------------------- the analytic cross-check

def real_to_complex(l: int) -> np.ndarray:
    """``U^(l)``: real harmonics (PySCF AO order) -> complex
    ``|l,m>``, ``m = -l..l`` -- the per-shell transform every M_L
    bookkeeping of real AOs needs (real shell AOs are ``+-|m|``
    combinations, not m eigenstates).

    Exact tesseral algebra (Condon-Shortley complex harmonics, the
    same convention as scipy/PySCF):

        sine slot (column l-m): (-1)^(m+1) * sqrt(2) * Im Y_l^m
        cos  slot (column l+m): (-1)^m     * sqrt(2) * Re Y_l^m

    with PySCF's p-shell reordering (px, py, pz) applied for ``l = 1``.

    Returns
    -------
    Unitary ``(2l+1, 2l+1)``: ``real_k = sum_m U[m, k] Y_l^m``.
    """
    k = 2 * l + 1
    u = np.zeros((k, k), dtype=complex)
    u[l, l] = 1.0
    sq2 = np.sqrt(2.0)
    for m in range(1, l + 1):
        rp, rm = l + m, l - m  # rows of Y^{+m} and Y^{-m}
        # sine-type slot at column l-m
        u[rp, l - m] = (-1) ** (m + 1) * 1j / sq2
        u[rm, l - m] = (-1) ** (m + 1) * (-1j * (-1) ** m) / sq2
        # cosine-type slot at column l+m
        u[rp, l + m] = (-1) ** m / sq2
        u[rm, l + m] = 1.0 / sq2
    if l == 1:  # PySCF AO order for p shells: (px, py, pz)
        u = u[:, [2, 0, 1]]
    return u


def _complex_generators(l: int) -> dict[str, np.ndarray]:
    """Lx, Ly, Lz in the complex |l,m> basis, m = -l..l:
    Lz diagonal, the ladder L+ strictly upper with the
    ``sqrt((l-m)(l+m+1))`` Clebsch factors, Lx/Ly its symmetric
    combinations."""
    k = 2 * l + 1
    m_vals = np.arange(-l, l + 1)
    lz = np.diag(m_vals).astype(complex)
    lp = np.zeros((k, k), dtype=complex)
    for row, m in enumerate(m_vals[:-1]):
        lp[row + 1, row] = np.sqrt((l - m) * (l + m + 1))  # m -> m + 1
    lm = lp.conj().T
    return {
        "x": (lp + lm) / 2.0,
        "y": (lp - lm) / (2.0j),
        "z": lz,
    }


def l_matrices_analytic(mol) -> dict[str, np.ndarray]:
    """AO matrices of Lx, Ly, Lz built from overlaps and exact
    solid-harmonic blocks -- the INDEPENDENT second construction that
    certifies :func:`l_matrices` elementwise (pinned by tests).

    Notes
    -----
    Works per contraction slot (``windows.shell_slots``): fused shells
    (``nctr > 1``) are handled exactly like segmented ones -- each
    slot is one radial times its m-set, so the generators act within
    the slot and the overlap does the transport.
    """
    _check_single_center(mol)
    nao = mol.nao_nr()
    s = mol.intor("int1e_ovlp")
    l_ops = {a: np.zeros((nao, nao), dtype=complex) for a in "xyz"}

    gen_cache: dict[int, dict[str, np.ndarray]] = {}
    for _ib, l, _c, rows in shell_slots(mol):
        if l == 0:
            continue
        if l not in gen_cache:
            u = real_to_complex(l)
            gen_cache[l] = {
                a: u.conj().T @ g @ u for a, g in _complex_generators(l).items()
            }
        for a in "xyz":
            blk = gen_cache[l][a]
            # L|nu> = sum_k blk[k, pos_nu] |chi_k>  (same slot => same radial)
            l_ops[a][:, rows] += s[:, rows] @ blk

    return l_ops
