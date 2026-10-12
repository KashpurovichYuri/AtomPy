"""One-electron angular momentum operators Lx, Ly, Lz, L^2 in the AO basis.

Primary construction (libcint, via PySCF)
-----------------------------------------
PySCF exposes the libcint integral ``int1e_irp`` with 9 components; decoding
(verified empirically and pinned by tests) gives

    T[3*a + b] = Integral mu | r_a d_b nu      (a, b = 0..2 for x, y, z)

The angular momentum operator ``L = -i r x grad`` is the antisymmetric part:

    Lz = -i (x dy - y dx) = -i (T[1] - T[3])
    Lx = -i (y dz - z dy) = -i (T[5] - T[7])
    Ly = -i (z dx - x dz) = -i (T[6] - T[2])

Exact for any pure-spherical AO basis on a single center at the origin.

Analytic cross-check construction
---------------------------------
For a single center at the origin every AO is a radial function times a real
spherical harmonic, so the generators act only within each (l, shell) block.
The real-to-complex harmonic transform ``U^(l)`` follows from exact tesseral
algebra in the Condon-Shortley convention (with PySCF's p-shell ordering),

    gen^(l) = U^dag gen_c^(l) U,   gen_c = generators in the |l,m> basis

``Lz|mu> = Sum_k blk[k, pos_mu] |chi_k>`` for k in the same shell, hence
``<mu|Lz|nu> = Sum_k blk[k, pos_nu] S[mu, k]``.  For ``L^2`` simply
``L^2 = Lx.Lx + Ly.Ly + Lz.Lz`` (matrix product in the orthonormalized AO
basis), which is exact because the single-center space is closed under L.

Limitation: valid for a single atom at the origin with spherical basis
(enforced by atomterm.system.build_mol).
"""

from __future__ import annotations

import numpy as np

__all__ = ["l_matrices", "l2_matrix", "l_matrices_analytic"]


def _check_single_center(mol) -> None:
    if mol.cart:
        raise ValueError("spherical (pure) basis required")
    if mol.natm != 1 or not np.allclose(mol.atom_coords()[0], 0.0, atol=1e-10):
        raise ValueError("operators are defined for a single atom at the origin")


def l_matrices(mol) -> dict[str, np.ndarray]:
    """AO matrices of Lx, Ly, Lz (complex, Hermitian, imaginary) via libcint."""
    _check_single_center(mol)
    t = mol.intor("int1e_irp", comp=9)
    img = -1j
    return {
        "x": img * (t[5] - t[7]),
        "y": img * (t[6] - t[2]),
        "z": img * (t[1] - t[3]),
    }


def l2_matrix(mol) -> np.ndarray:
    """AO matrix of the one-electron L^2 = Lx.Lx + Ly.Ly + Lz.Lz (exact).

    Matrix products are exact only in an orthonormal basis, so the product is
    taken in the symmetrically orthonormalized AO basis and transformed back:
    L2_ao = S^(1/2) [ sum_a (S^-1/2 l_a S^-1/2)^2 ] S^(1/2).
    """
    _check_single_center(mol)
    s = mol.intor("int1e_ovlp")
    ev, u = np.linalg.eigh(s)
    s_half = (u * np.sqrt(ev)) @ u.T
    s_ihalf = (u / np.sqrt(ev)) @ u.T
    l_ops = l_matrices(mol)
    l2_orth = sum(
        (s_ihalf @ l @ s_ihalf) @ (s_ihalf @ l @ s_ihalf) for l in l_ops.values()
    )
    return s_half @ l2_orth @ s_half


def _shell_blocks(mol):
    """Yield (l, first_ao_index, n_components) for every shell."""
    ao_loc = mol.ao_loc_nr()
    for k in range(mol.nbas):
        l = mol.bas_angular(k)
        yield l, ao_loc[k], ao_loc[k + 1] - ao_loc[k]


def _real_to_complex(l: int) -> np.ndarray:
    """U^(l): real-harmonic (PySCF AO order) -> complex |l,m>, m=-l..l.

    Exact tesseral algebra (Condon-Shortley complex harmonics, the same
    convention as scipy/PySCF):

        sine slot (column l-m): (-1)^(m+1) * sqrt(2) * Im Y_l^m
        cos  slot (column l+m): (-1)^m     * sqrt(2) * Re Y_l^m

    with PySCF's p-shell reordering (px, py, pz) applied for l = 1.
    Columns are the coefficients of each real harmonic in the complex |l,m>
    basis, i.e. real_k = sum_m U[m, k] Y_l^m.
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
    """Lx, Ly, Lz in the complex |l,m> basis, m = -l..l."""
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
    """AO matrices of Lx, Ly, Lz built from overlaps + solid-harmonic blocks.

    Independent of the libcint path; used to cross-validate ``l_matrices``.
    """
    _check_single_center(mol)
    nao = mol.nao_nr()
    s = mol.intor("int1e_ovlp")
    l_ops = {a: np.zeros((nao, nao), dtype=complex) for a in "xyz"}

    gen_cache: dict[int, dict[str, np.ndarray]] = {}
    for l, ao0, ncomp in _shell_blocks(mol):
        if l == 0:
            continue
        if l not in gen_cache:
            u = _real_to_complex(l)
            gen_cache[l] = {
                a: u.conj().T @ g @ u for a, g in _complex_generators(l).items()
            }
        for a in "xyz":
            blk = gen_cache[l][a]
            # Lz |nu> = sum_k blk[k, pos_nu] |chi_k>  (same shell => same radial)
            l_ops[a][:, ao0 : ao0 + ncomp] += s[:, ao0 : ao0 + ncomp] @ blk

    return l_ops
