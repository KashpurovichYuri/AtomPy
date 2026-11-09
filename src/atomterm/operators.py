"""Operators: the first-quantized assembly layer.

PURPOSE
-------
bridge PySCF (AO/MO integrals) and angular (the exact L
matrices) to detbasis' determinant space, producing every operator
the term analysis reads as a plain ``(ndet, ndet)`` matrix -- H, Lz,
Lz^2, L^2, S^2, Sz, N, the parity and population labels, and the
closure certificate. An operator enters the determinant space as a
(one-body matrix, two-body tensor) pair; either twin evaluator turns
the pair into the matrix -- detbasis' Slater-Condon rules or
secondquant's a+/a strings, elementwise identical behind
``build_operators(method=...)`` (pinned by tests).

Symbols
-------
Every operator is the pair (o, g): ``O1 = Σ_pq o_pq a+_p a_q`` and
``O2 = ½ Σ_pqrs g_pqrs a+_p a+_q a_s a_r`` (spin-orbitals ordered
(alpha; beta); ``L^a = C† L^a_AO C`` -- the AO angular-momentum
matrices transformed into the window)::

    N   : o = 1
    Lz  : o = L^z
    Lz^2: o = (L^z)^2,  g = 2 (L^z)_pr (L^z)_qs
    L^2 : o = Σ_a (L^a)^2,  g = 2 Σ_a (L^a)_pr (L^a)_qs
    Sz  : o = diag(½·1, −½·1)
    S^2 : o = ¾·1,  g = 2[(s_z)_pr (s_z)_qs + ½(s+)_pr (s−)_qs
          + ½(s−)_pr (s+)_qs]
    L±  : o = L^x ± i L^y          (orbital ladders; sector operators)
    S±  : o = s±  (s+ = |α><β|)    (inter-sector maps: M_S → M_S ± 1)
    H   : o = h_eff (the folded core field),  g = (pq|rs)

The ladder algebra is pinned elementwise: L^2 = Lz^2 + ½(L+L− + L−L+),
S^2 = Sz^2 + ½(S−S+ + S+S−), [Lz, L±] = ±L±. Sz certifies nothing
(the field is spin-free: [H, Sz] ≡ 0 identically) -- it labels the
sector; its matrix is exactly M_S·1.

Conventions
-----------
* Spin orbitals follow detbasis: alpha ``0..norb-1``, beta
  ``norb..2norb-1``.
* Two-body tensors are detbasis' physicist tensor
  ``g[p,q,r,s] = <pq|rs> = int p(1) q(2) r(1) s(2)``
  = chemist ``(pr|qs)``.
* PySCF is chemist everywhere (``int2e``, ``ao2mo``,
  ``fci.direct_spin1``). The transpose happens exactly ONCE, in
  ``_physicist_eri``; passing ``g`` to PySCF FCI needs
  ``g.transpose(0, 2, 1, 3)``.
* PySCF FCI vectors order determinants by (alpha-string address,
  beta-string address) in integer-ascending bit order, while detbasis
  enumerates lexicographic combinations -- the two coincide only in
  trivial sectors (norb <= 3, or <= 1 electron per spin block) and
  differ everywhere else (first at norb=4, na=2).
  :func:`pyscf_ci_permutation` is the ONE bridge: every consumer that
  indexes or contracts a PySCF CI vector against detbasis-order
  objects maps through it first.
"""

from __future__ import annotations

import numpy as np

from atomterm.angular import l_ladder_matrices, l_matrices
from atomterm.detbasis import (
    determinant_basis,
    frobenius_norm,
    slater_condon_matrix,
)
from atomterm.secondquant import annihilate, create, secondquant_matrix
from atomterm.windows import (
    l_per_row,
    occupied_core_idx,
    s_metric_roots,
    shell_slots,
)

__all__ = [
    "spin_block_one_body",
    "spin_block_two_body",
    "active_space_integrals",
    "full_space_integrals",
    "lz_operator",
    "l_ladder_operator",
    "parity_one_body",
    "parity_operator",
    "population_operator",
    "l2_operator_parts",
    "window_closure",
    "core_anisotropy",
    "ground_block",
    "pyscf_ci_permutation",
    "ci_in_detbasis_order",
    "pyscf_apply",
    "s2_operator_parts",
    "spin_ladders",
    "spin_ladder_matrix",
    "sz_operator",
    "number_operator",
    "build_operator_parts",
    "build_operators",
]


# -------------------------------------------------------- the PySCF boundary


def spin_block_one_body(a_spatial: np.ndarray) -> np.ndarray:
    """Spatial one-body matrix -> spin-orbital block-diagonal form.

    Returns
    -------
    (2n, 2n) complex ``diag(a, a)``: the same operator lifted to the
    spin-orbital space (alpha block first, beta second).

    Examples
    --------
    >>> float(spin_block_one_body(np.eye(2))[1, 1].real)
    1.0
    """
    z = np.zeros_like(a_spatial)
    return np.block([[a_spatial, z], [z, a_spatial]]).astype(complex)


def spin_block_two_body(g_spatial: np.ndarray) -> np.ndarray:
    """Spatial physicist tensor -> spin-orbital physicist tensor.

    In ``g[p,q,r,s] = <pq|rs>`` the index pairs (p, r) and (q, s)
    ride the same electron, so spin conservation keeps only the
    blocks (aa;aa), (ab;ab), (ba;ba), (bb;bb).
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
    """AO-basis operator expressed in the window basis (``C^dag A C``)."""
    return c_act.conj().T @ a_ao @ c_act


def _physicist_eri(eri: np.ndarray, c: np.ndarray) -> np.ndarray:
    """AO chemist ERIs -> physicist tensor in the orbital basis ``c``.

    The ONE transpose at the PySCF convention boundary:
    ``(uv|wx)  ->  g[p,q,r,s] = <pq|rs> = chemist (pr|qs)``.
    """
    return np.einsum(
        "up,vq,wr,xs,uvwx->pqrs",
        c.conj(),
        c.conj(),
        c,
        c,
        eri,
        optimize=True,
    ).transpose(0, 2, 1, 3)


def _spherical_average_density(mol, p: np.ndarray) -> np.ndarray:
    """l=0 (totally symmetric) component of an AO density matrix.

    Within every pair of same-l slots (:func:`windows.shell_slots`)
    the block is replaced by its mean over the shared m index,
    ``P[(c1,m),(c2,m')] -> delta_mm' mean_m P[(c1,m),(c2,m)]``;
    cross-l blocks vanish by exact symmetry and are zeroed.

    Notes
    -----
    An atomic core density is spherical by construction: this
    removes nothing but SCF numerical noise and restores exact
    rotational invariance of the analysis Hamiltonian.
    """
    ps = np.zeros_like(p)
    slots = [(l, rows) for _, l, _, rows in shell_slots(mol)]
    for l1, r1 in slots:
        for l2, r2 in slots:
            if l1 != l2:
                continue  # zero by exact symmetry
            nk = len(r1)
            mean_m = np.mean([p[r1[j], r2[j]] for j in range(nk)])
            for m in range(nk):
                ps[r1[m], r2[m]] = mean_m
    return ps


def active_space_integrals(
    mol, mf, c_act: np.ndarray, ncore: int, spherical_core: bool = True,
    c_core: np.ndarray = None,
):
    """(h1, g) of the active space with the core folded into h1.

    The mcscf convention: ``h1 = C^dag (h_core + V_eff(core)) C``,
    ``V_eff[m,n] = sum_ls P[l,s] (mn|ls) - 0.5 (ml|ns)``. CI
    eigenvalues of (h1, g) equal ``mc.e_cas``.

    Parameters
    ----------
    c_act : (nao, ncas) active columns (S-orthonormal).
    ncore : number of folded core pairs.
    spherical_core : average the core density to its l=0 component
        (the default; see Notes).
    c_core : folded core columns. Default: the ROHF CLOSED block
        (``occupied_core_idx`` -- canonical order is junk on
        near-degenerate uncontracted bases). Optimized stages MUST
        pass their own converged inactive columns: folding a
        different core than the energies were computed with breaks
        ``E_tot = E_core + e_cas`` at O(Eh), not small drift.

    Returns
    -------
    (h1, g): h1 carries the core field; g is the raw physicist
    tensor of the module conventions (no folding).

    Notes
    -----
    ``spherical_core=False`` reproduces the raw SCF core field, which
    carries the reference's anisotropy: the core is polarized by ONE
    member of the degenerate open shell (formally ``H_eff = H0 +
    eps*W``, eps ~ 1e-5). Then ``[H, L] = O(eps)`` and e.g. the 3P
    multiplet splits as {M_L=+-1}/{0} (a k=2, m=0 field; the k=1
    dipole is excluded exactly by inversion symmetry). The averaged
    density is the ensemble average over the degenerate manifold --
    the scalar projection -- and restores exact invariance.
    """
    c_core = (mf.mo_coeff[:, occupied_core_idx(mf, ncore)]
              if c_core is None else c_core)
    d_raw = 2.0 * c_core @ c_core.conj().T
    d_core = _spherical_average_density(mol, d_raw) if spherical_core else d_raw

    eri = mol.intor("int2e", aosym="s1")  # chemist (mu nu | la si)
    veff_core = (
        np.einsum("ls,mnls->mn", d_core, eri, optimize=True)
        - 0.5 * np.einsum("ls,mlns->mn", d_core, eri, optimize=True)
    )
    h1_ao = mf.get_hcore() + veff_core
    h1 = _mo_transform(h1_ao, c_act)
    return h1, _physicist_eri(eri, c_act)


def full_space_integrals(mol, mf, w: np.ndarray):
    """(h1, g) over the given orbital columns with NO core folding.

    Raw ``h_core`` and raw ERIs (physicist tensor) in the orbital
    basis ``w`` -- the pieces of the all-electron (core-explicit)
    Hamiltonian.

    Parameters
    ----------
    w : (nao, k) orbital columns; the full AO identity gives the
        AO-basis integrals themselves.
    """
    h1 = _mo_transform(mf.get_hcore(), w)
    g = _physicist_eri(mol.intor("int2e", aosym="s1"), w)
    return h1, g


def pyscf_ci_permutation(norb: int, na: int, nb: int) -> np.ndarray:
    """Map ravelled PySCF CI amplitudes onto the detbasis ordering.

    PySCF's direct_spin1 vector ``v`` (shape ``(nalpha_strings,
    nbeta_strings)``, ravelled alpha-major) indexes determinants by
    string ADDRESSES (integer-ascending bit strings); detbasis
    enumerates lexicographic combinations. These orders differ for
    norb >= 4 once a block holds >= 2 electrons -- coinciding only in
    trivial sectors.

    Returns ``perm`` with ``v[perm[k]] == w_detbasis[k]``-semantics:
    the pyscf amplitude at ravelled index ``k`` belongs to determinant
    ``determinant_basis(norb, na, nb)[perm[k]]``. Scatter with
    ``w[perm] = v`` (or use :func:`ci_in_detbasis_order`).

    Certified by the eigen-equation test: a PySCF-solved vector,
    mapped and contracted against the independently built
    Slater-Condon Hamiltonian, satisfies ``H w = e w`` -- generic
    integrals admit no permutation automorphism, so a wrong mapping
    cannot pass.
    """
    from pyscf.fci import cistring

    dets = determinant_basis(norb, na, nb)
    index = {d: i for i, d in enumerate(dets)}
    a_strs = cistring.make_strings(range(norb), na)  # index == address
    b_strs = cistring.make_strings(range(norb), nb)
    perm = np.empty(len(dets), dtype=np.int64)
    k = 0
    for sa in a_strs:
        occ_a = [i for i in range(norb) if sa >> i & 1]  # orbital i = bit i
        for sb in b_strs:
            occ_b = [i for i in range(norb) if sb >> i & 1]
            det = tuple(sorted(occ_a)) + tuple(o + norb for o in occ_b)
            perm[k] = index[det]
            k += 1
    return perm


def ci_in_detbasis_order(v: np.ndarray, norb: int, na: int,
                         nb: int) -> np.ndarray:
    """One ravelled PySCF CI vector, reordered to detbasis indexing."""
    v = np.asarray(v).ravel()
    w = np.empty_like(v)
    w[pyscf_ci_permutation(norb, na, nb)] = v
    return w


def pyscf_apply(h1e: np.ndarray, eri: np.ndarray, norb: int, na: int,
                nb: int):
    """The C-contraction boundary: the Hamiltonian applied through
    PySCF's own kernels at ANY sector size -- no dense operator, no
    cached couplings, the permutation owned here.

    Parameters
    ----------
    h1e, eri
        The KERNEL FEED's convention (what the lane evaluators pass
        the solvers): ``h1e`` from :func:`active_space_integrals`
        (the folded active one-body), ``eri`` from
        ``ao2mo.kernel(mol, c)`` (PySCF-native chemist).  Real
        tensors -- PySCF's contraction kernels refuse complex input.
    norb, na, nb
        The window width and the sector.

    Returns
    -------
    ``apply`` -- a closure mapping a sector vector in detbasis order
    (real or complex) to ``H @ v`` (same order; a complex vector
    contracts through its real and imaginary halves).  ``h1e`` is
    absorbed into the contracted tensor at the kernel's own ``0.5``
    factor (:func:`pyscf.fci.direct_spin1.absorb_h1e`) -- the exact
    pairing whose roots the lane anchor certifies against our
    Slater-Condon Hamiltonian; the apply pin certifies the
    vector-level action.

    Notes
    -----
    The one-body counterpart does NOT exist here on purpose: the C
    kernel packs the lower triangle, so ``contract_1e`` silently
    symmetrizes -- it serves symmetric matrices only, and the
    angular-momentum operators (antisymmetric-imaginary in a real
    basis) would be corrupted.  Apply those through
    :func:`atomterm.detbasis.operator_apply`, with ``Lz`` rotated to
    its real form ``M = -i Lz`` (norms unchanged) when a real
    intermediate is wanted.
    """
    from pyscf.fci import direct_spin1

    perm = pyscf_ci_permutation(norb, na, nb)
    h2e = direct_spin1.absorb_h1e(
        np.ascontiguousarray(np.real(h1e), dtype=float),
        np.ascontiguousarray(np.real(eri), dtype=float),
        norb, (na, nb), .5)

    def _real(x):
        out = direct_spin1.contract_2e(
            h2e, np.ascontiguousarray(np.real(x)[perm]),
            norb, (na, nb)).reshape(-1)
        back = np.empty(out.size, dtype=out.dtype)
        back[perm] = out
        return back

    def apply(v):
        v = np.asarray(v)
        out = _real(v)
        if np.iscomplexobj(v):
            out = out + 1j * _real(v.imag)
        return out

    return apply


# --------------------------------------------------------- one-body symbols


def lz_operator(mol, c_act: np.ndarray) -> np.ndarray:
    """Lz as a spin-orbital one-body matrix in the window basis
    (spatially diagonal: L never touches spin)."""
    return spin_block_one_body(_mo_transform(l_matrices(mol)["z"], c_act))


def l_ladder_operator(mol, c_act: np.ndarray) -> dict:
    """One-body symbols of the orbital ladders ``L+- = Lx +- i Ly`` in
    the window basis, spin-blocked.

    Many-body L+- are ordinary sector operators (orbital ladders do
    not touch spin); they power the twin construction
    ``L^2 = Lz^2 + (L+ L- + L- L+)/2`` (pinned by test).

    Returns
    -------
    dict ``{"+": ..., "-": ...}``.
    """
    ladders = l_ladder_matrices(mol)
    out = {}
    for key in ("+", "-"):
        out[key] = spin_block_one_body(_mo_transform(ladders[key], c_act))
    return out


def sz_operator(norb: int) -> np.ndarray:
    """One-body S_z over spin orbitals: +1/2 alpha, -1/2 beta.

    Already spin-orbital-resolved (do not spin_block). Inside a fixed
    (na, nb) sector its many-body matrix is exactly ``M_S * identity``
    -- the sector label.

    Examples
    --------
    >>> sz_operator(2).diagonal().tolist()
    [0.5, 0.5, -0.5, -0.5]
    """
    sz = np.zeros((2 * norb, 2 * norb))
    diag = np.arange(2 * norb)
    sz[diag[:norb], diag[:norb]] = 0.5
    sz[diag[norb:], diag[norb:]] = -0.5
    return sz


def spin_ladders(norb: int):
    """One-body ladder symbols ``s+ = |alpha><beta|`` (the beta ->
    alpha flip, ``<p a| s+ |p b> = 1``) and ``s- = s+^dag``.

    Bare ladders change M_S: building blocks of the S^2 pair terms,
    not sector operators.

    Examples
    --------
    >>> sp, _ = spin_ladders(1)
    >>> float(abs(sp[0, 1]))
    1.0
    """
    d = 2 * norb
    s_plus = np.zeros((d, d), dtype=complex)
    rows = np.arange(norb)
    s_plus[rows, rows + norb] = 1.0  # <p alpha| s+ |p beta> = 1
    return s_plus, s_plus.conj().T


def number_operator(norb: int) -> np.ndarray:
    """Electron number N: the identity, spin-blocked."""
    return spin_block_one_body(np.eye(norb))


def spin_ladder_matrix(norb: int, dets_from, dets_to, plus: bool = True) -> np.ndarray:
    """Many-body S+- as an explicit INTER-SECTOR map.

    ``S+ : (na, nb) -> (na+1, nb-1)`` (the beta -> alpha flip),
    ``S- : (na, nb) -> (na-1, nb+1)``;
    ``S+ = sum_i a+(i alpha) a(i beta)``.

    The ONE operator class with no (one_body, two_body) form: an
    inter-sector map has no sector-symbol representation, so a+/a are
    applied directly (the documented convention exception). Built for
    the twin construction ``S^2 = Sz^2 + (S- S+ + S+ S-)/2`` (pinned
    by test).

    Parameters
    ----------
    dets_from, dets_to : determinant lists of the source and target
        sectors (``determinant_basis``).
    plus : True builds S+, False builds S-.

    Returns
    -------
    ``(len(dets_to), len(dets_from))`` dense matrix.
    """
    out = np.zeros((len(dets_to), len(dets_from)))
    index = {d: k for k, d in enumerate(dets_to)}
    for col, occ in enumerate(dets_from):
        for i in range(norb):
            p, q = (i, i + norb) if plus else (i + norb, i)
            if q not in occ or p in occ:
                continue
            sg1, mid = annihilate(occ, q)
            sg2, new = create(mid, p)
            out[index[new], col] = sg1 * sg2
    return out


# --------------------------------------------------------- operator parts


def _composite_two_body(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Physicist tensor of the product operator ``A(1)B(2)``:
    ``g[p,q,r,s] = A[p, r] B[q, s]``."""
    return np.einsum("pr,qs->pqrs", a, b)


def l2_operator_parts(mol, c_act: np.ndarray):
    """(one_body, two_body) parts of the total L^2.

    ``L^2 = sum_i l^2(i) + 2 sum_{i<j} l_j . l_i``: one-body
    ``sum_a l_a . l_a`` + two-body ``2 sum_a l_a (x) l_a``.

    Notes
    -----
    The products are the WINDOW-restricted one-particle operators --
    exact as such in any S-orthonormal window basis. Whether the
    restricted L^2 equals the full one is the window's closure
    (:func:`window_closure`), which the many-body ``[H, L^2]``
    commutators measure.
    """
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

    ``S^2 = Sz^2 + (S+ S- + S- S+)/2`` built from the first-quantized
    composites ``sum_{i<j} A(i) B(j)``; equivalently
    ``sum_i 3/4 + 2 sum_{i<j} s_j . s_i``.
    """
    d = 2 * norb
    sz = sz_operator(norb)
    s_plus, s_minus = spin_ladders(norb)

    one = 0.75 * np.eye(d, dtype=complex)
    two = 2.0 * (
        _composite_two_body(sz, sz)
        + 0.5 * _composite_two_body(s_plus, s_minus)
        + 0.5 * _composite_two_body(s_minus, s_plus)
    )
    return one, two


# ---------------------------------------------------------------- the labels


def parity_one_body(mol, c_act: np.ndarray) -> np.ndarray:
    """The parity operator ``(-1)^l`` in the window basis -- EXACT.

    On the one-particle space parity is ``sum_l (-1)^l P_l`` with
    ``P_l`` the l-subspace projector; l-subspaces are exactly
    S-orthogonal on one center, so in the AO basis

    ``Pi_AO[m, n] = (-1)^{l(m)} S[m, n] delta_{l(m), l(n)}``

    (every AO row has a definite l, :func:`windows.l_per_row`), and
    the window compression is ``Pi_win = C^dag Pi_AO C``.

    Notes
    -----
    l-pure columns give a diagonal +-1 matrix. An l-mixed column is
    DATA, not an error: the eigenvalues leave +-1 (a column half s
    half p pushes the window eigenvalues to ``(1 +- sqrt(3))/2``) --
    the operator is exact, the deviation belongs to the window
    object.
    """
    s = mol.intor("int1e_ovlp")
    ell = l_per_row(mol)
    same_l = ell[:, None] == ell[None, :]
    pi_ao = np.where(same_l, (-1.0) ** ell[:, None] * s, 0.0)
    return c_act.conj().T @ pi_ao @ c_act


def parity_operator(mol, c_act: np.ndarray, na: int, nb: int) -> np.ndarray:
    """Many-body parity: the exterior power ``Lambda^N(Pi)``, exact.

    Matrix elements between determinants are the Cauchy-Binet minors

    ``<D'|Lambda^N(Pi)|D> = det(Pi[D'_a, D_a]) det(Pi[D'_b, D_b])``

    over the occupied spatial column sets (alpha and beta blocks
    factor: Pi is spin-diagonal). On a parity-CLOSED window (Pi maps
    the span into itself -- every l-pure window) this is diagonal
    with the product of occupied signs, squares to I, and has
    eigenvalues exactly ``+-1``. On a mixed window the object is the
    honest ``Lambda^N`` of the compressed Pi: dense, eigenvalues off
    ``+-1`` -- sharp parity is lost, and that loss is the window's
    deviation, not the operator's. Like the spin ladders,
    multiplicative -- outside the additive (one_body, two_body)
    parts protocol, hence built directly over the determinants.
    """
    pi = parity_one_body(mol, c_act)
    norb = c_act.shape[1]
    dets = determinant_basis(norb, na, nb)
    out = np.zeros((len(dets), len(dets)))
    cols = [
        ([d % norb for d in det if d < norb],
         [d % norb for d in det if d >= norb])
        for det in dets
    ]
    for k2, (a2, b2) in enumerate(cols):
        for k1, (a1, b1) in enumerate(cols):
            out[k2, k1] = (np.linalg.det(pi[np.ix_(a2, a1)])
                           * np.linalg.det(pi[np.ix_(b2, b1)]))
    return out


def population_operator(mol, c_act: np.ndarray, rows: list[int], na: int,
                        nb: int) -> np.ndarray:
    """Electron count with S-weight in the given AO ``rows``, as a
    determinant-basis matrix -- the per-root ``<N_f>`` / ``<N_d>`` of
    configuration tables.

    The one-body part is the S-metric projector onto the rows' span,
    ``A = S Sel (Sel^dag S Sel)^{-1} Sel^dag S`` (Hermitian, spectrum
    {0, 1}), evaluated like every other one-body operator. Rows of
    different l are exactly S-orthogonal on one center, so for an
    l-disjoint window (e.g. f rows + d rows) the population operators
    PARTITION the number operator: ``population(f) + population(d)
    == N`` elementwise.
    """
    sel = np.eye(mol.nao_nr())[:, list(rows)]
    s = mol.intor_symmetric("int1e_ovlp")
    a = s @ sel @ np.linalg.solve(sel.T @ s @ sel, sel.T @ s)
    one = spin_block_one_body(_mo_transform(a, c_act))
    dets = determinant_basis(c_act.shape[1], na, nb)
    return slater_condon_matrix(dets, one_body=one, two_body=None)


# --------------------------------------------------------- the certificates


def window_closure(mol, c_win: np.ndarray) -> dict[str, float]:
    """One-particle rotational-closure defects of a window span.

    The span of the S-orthonormal columns ``c_win`` is
    rotation-closed iff its projector commutes with all three
    rotation generators. In the symmetrically orthonormalized AO
    basis (the :func:`windows.s_metric_roots` transport) the
    projector is ``P = Q Q^dag`` with ``Q = S^{1/2} c_win``
    (QR-cleaned) and the generators act as ``S^{-1/2} l_a S^{-1/2}``.

    Returns
    -------
    dict of Frobenius norms of the commutators:

    ``"lz"``
        ``||[P, Lz]||`` -- exactly the part the many-body ``[H, Lz]``
        can see. Lz mixes only +-m pairs: m=0-member tilts and EQUAL
        tilts of a full +-m pair are INVISIBLE here (the fact-18c
        blind spots).
    ``"lperp"``
        ``sqrt(||[P, Lx]||^2 + ||[P, Ly]||^2)`` -- the defect only
        L+- (hence L^2) sees.

    Notes
    -----
    ``lz = lperp = 0`` (to roundoff) is the COMPLETE closure
    certificate: no blind spots, no core fold, no CI -- symmetry
    breaking is attributed to the span itself. The many-body
    commutator norms respond to these defects at O(defect), the
    energy at O(defect^2) at stationary points.
    """
    s = mol.intor("int1e_ovlp")
    s_sqrt, s_isqrt = s_metric_roots(s)
    q, _ = np.linalg.qr(s_sqrt @ c_win)
    p = q @ q.conj().T
    defect = {}
    for axis, l_ao in l_matrices(mol).items():
        a = s_isqrt @ l_ao @ s_isqrt
        defect[axis] = frobenius_norm(p @ a - a @ p)
    return {
        "lz": float(defect["z"]),
        "lperp": float(np.hypot(defect["x"], defect["y"])),
    }


def core_anisotropy(mol, mf, ncore: int, c_core: np.ndarray = None) -> float:
    """Relative weight of the discarded core-field anisotropy,
    ``||rho - rho_sph||_F / ||rho_sph||_F`` -- the eps of
    ``H_eff = H0 + eps*W`` (see :func:`active_space_integrals`).

    Parameters
    ----------
    c_core : folded core columns (default: as in
        :func:`active_space_integrals`).
    """
    c_core = (mf.mo_coeff[:, occupied_core_idx(mf, ncore)]
              if c_core is None else c_core)
    d = 2.0 * c_core @ c_core.conj().T
    d_sph = _spherical_average_density(mol, d)
    return float(frobenius_norm(d - d_sph) / frobenius_norm(d_sph))


def ground_block(ops: dict) -> tuple[float, int, np.ndarray, np.ndarray]:
    """The frozen ground-term block of an operator set, MEASURED: the
    exact eigenspace at ``ev[0]`` sharing its ``<L^2>``.

    Returns ``(ev0, nroots, weights, vec)``: the ground energy, the
    block size, uniform weights over the block, and the eigenvectors
    of the solve (returned so callers that need them -- e.g. the NO
    source density -- never repeat the eigh). One helper for every
    consumer that needs "the ground term's states" rather than "the
    lowest nroots states": the SA-anchoring rule and the NO-source
    averaging rule (never hardcode a term).

    The energy window is SCALE-AWARE: a residual core-field anisotropy
    splits the multiplet's energies AND perturbs the states' ``<L^2>``
    (they are no longer exact eigenstates of a non-scalar H), so the
    term is identified by L^2 within 0.1 -- far above the anisotropy
    perturbation, far below any term spacing (>= 2 between L values)
    -- and the energy window spans one percent of the gap to the
    first state OUTSIDE that L^2 band.  When every near state shares
    the band, the window falls back to roundoff (1e-6).

    Parameters
    ----------
    ops : a :func:`build_operators` dict -- the block is measured on
        that operator set.

    Notes
    -----
    A same-``L^2`` excited term QUASI-DEGENERATE with the ground term
    (within the window) joins the block -- there the "ground block"
    is genuinely ill-defined, whatever the window; the L^2 identity,
    not the energy, carries the term.
    """
    ev, vec = np.linalg.eigh(ops["H"].real)
    nchk = min(60, len(ev))
    l2 = [vec[:, r] @ ops["L2"].real @ vec[:, r] for r in range(nchk)]
    other = [r for r in range(1, nchk) if abs(l2[r] - l2[0]) >= 0.1]
    tol = max(1e-6, 0.01 * (ev[other[0]] - ev[0])) if other else 1e-6
    block = [r for r in range(nchk)
             if abs(l2[r] - l2[0]) < 0.1 and ev[r] - ev[0] < tol]
    w = np.zeros(len(ev))
    w[block] = 1.0 / len(block)
    return float(ev[0]), len(block), w, vec


# -------------------------------------------------------------- the assembly


def build_operator_parts(
    mol, mf, c_act: np.ndarray, ncore: int, na: int, nb: int,
    spherical_core: bool = True, c_core: np.ndarray = None,
) -> dict:
    """Operator content of :func:`build_operators` as raw parts:
    ``{name: (one_body, two_body)}`` plus ``dets`` / ``h1`` / ``g``.

    The evaluation-agnostic seam behind ``method=...``: both twin
    evaluators consume these parts unchanged.
    """
    h1, g = active_space_integrals(mol, mf, c_act, ncore, spherical_core, c_core)
    dets = determinant_basis(c_act.shape[1], na, nb)
    s2_one, s2_two = s2_operator_parts(c_act.shape[1])
    l2_one, l2_two = l2_operator_parts(mol, c_act)
    lz_ao = _mo_transform(l_matrices(mol)["z"], c_act)
    parts = {
        "dets": dets,
        "h1": h1,
        "g": g,
        "H": (spin_block_one_body(h1), spin_block_two_body(g)),
        "Lz": (spin_block_one_body(lz_ao), None),
        # (Lz)^2 as an OPERATOR (the twin-convention square, cf. L2 above):
        # one-body lz@lz PLUS the two-body composite -- the naive
        # one-body-only square is not the operator (intermediate q!=r
        # terms survive the normal ordering)
        "Lz2": (spin_block_one_body(lz_ao @ lz_ao),
                spin_block_two_body(2.0 * _composite_two_body(lz_ao, lz_ao))),
        "L2": (l2_one, l2_two),
        "S2": (s2_one, s2_two),
        "Sz": (sz_operator(c_act.shape[1]), None),
        "N": (number_operator(c_act.shape[1]), None),
    }
    return parts


def build_operators(
    mol, mf, c_act: np.ndarray, ncore: int, na: int, nb: int,
    method: str = "slater_condon", spherical_core: bool = True,
    c_core: np.ndarray = None,
) -> dict:
    """All analysis operators as full matrices in the determinant basis.

    Parameters
    ----------
    method : ``"slater_condon"`` (default; explicit excitation rules)
        or ``"secondquant"`` (a+/a strings on occupation vectors).
        Both consume identical parts and produce elementwise-identical
        matrices (pinned by tests).
    spherical_core, c_core : see :func:`active_space_integrals`.

    Returns
    -------
    dict: ``"dets"``, ``"h1"``, ``"g"`` plus the matrices ``"H"``,
    ``"Lz"``, ``"Lz2"``, ``"L2"``, ``"S2"``, ``"Sz"``, ``"N"``.
    """
    parts = build_operator_parts(mol, mf, c_act, ncore, na, nb,
                                 spherical_core, c_core)
    dets = parts.pop("dets")
    integrals = (parts.pop("h1"), parts.pop("g"))

    if method == "slater_condon":
        def build(ob, tb):
            return slater_condon_matrix(dets, one_body=ob, two_body=tb)
    elif method == "secondquant":
        def build(ob, tb):
            return secondquant_matrix(dets, one_body=ob, two_body=tb)
    else:
        raise ValueError(f"unknown method: {method!r}")

    out = {"dets": dets, "h1": integrals[0], "g": integrals[1]}
    for name, (ob, tb) in parts.items():
        out[name] = build(ob, tb)
    return out
