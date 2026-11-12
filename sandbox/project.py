"""STAGE 1b (ticket 37): the projection primitive + the synthetic
battery (anchors 1-4) -- the SANDBOX form of windows.project_orbitals.

The spec (ticket 37, Implementation Decisions):

    C~_B = S_B^{-1} X C_A          (thresh=None: the naive foil)
    C~_B = S_kept^{-1} X C_A       (thresh given: the conditioning
                                    truncation INSIDE the projector)
    then Loewdin-orthonormalization of the projected span in B's
    metric, with the machine-zero capture directions dropped (the
    disclosure).

The capture of source column i:  c_i^T X^T S_B^{-1} X c_i  -- the
squared cosine between the source column and its best approximation
in B's span; exactly 1 under span containment.

Battery (this file's main):
  1  B = A identity            -- exact, zero parameters
  2  B >= span(A)              -- the projected columns reproduce A's
                                  exactly (the strongest anchor)
  3  H/He/C hand-checks        -- real small bases, capture high
  4  near-singular stress      -- the naive foil amplifies, the
                                  truncated solve stays honest

Sandbox protocol: in-repo, imports src but is imported by nothing;
deleted at integration (the primitive lands in windows.py beside
build_full_mo, the anchors port to tests/unit).
"""
import numpy as np
from pyscf import gto, scf

# ----------------------------------------------------------------- the primitive


def project_orbitals(mol_b, mol_a, c_a, thresh=None):
    """The cross-basis projector (the spec's interface, sandbox form).

    Returns (c_b, report): the S_B-orthonormal projected columns and
    the disclosure dict (per-column capture, the S_B spectrum at the
    cut, modes dropped, the projected span's pre-normalization
    eigenvalues, the orthonormality error).
    """
    s = mol_b.intor("int1e_ovlp")
    # (nao_b, nao_a): the cross overlap, A's AOs paired against B's
    x = gto.intor_cross("int1e_ovlp", mol_a, mol_b).T
    rhs = x @ c_a
    ev, v = np.linalg.eigh(s)
    if thresh is None:
        w = np.linalg.solve(s, rhs)                  # the naive foil
        n_drop = 0
    else:
        keep = ev > thresh
        w = (v[:, keep] / ev[keep]) @ (v[:, keep].T @ rhs)
        n_drop = int((~keep).sum())
    # the capture: the squared cosine of source column vs its best
    # approximation in B (columnwise: c_i^T X^T W_i)
    cap = np.einsum("ij,ij->j", c_a.conj(), (x.T @ w))
    cap = cap.real
    # the per-AO participation of the projection (row-wise S-norms of
    # w) -- the membership meter: on a near-degenerate pair the NAIVE
    # solve assigns the weight to an arbitrary member (rounding
    # decides), the truncated solve to the symmetric combination
    w_rows = np.sqrt(np.maximum(np.einsum("ij,ij->i", w.real,
                                          (s @ w).real), 0.0))
    # Loewdin within the projected span, B's metric; drop machine-zero
    m = w.T @ s @ w
    m = 0.5 * (m + m.T)
    mu, u = np.linalg.eigh(m)
    good = mu > 1e-12
    # the SYMMETRIC inverse square root m^{-1/2} = U mu^{-1/2} U^T --
    # dropping the trailing U^T leaves a gauge-chaotic rotation on
    # near-degenerate spectra (caught by the identity anchor)
    m_isqrt = (u[:, good] / np.sqrt(mu[good])) @ u[:, good].T
    c_b = w @ m_isqrt
    ortho = float(np.abs(c_b.T @ s @ c_b - np.eye(c_b.shape[1])).max())
    report = {"capture": cap, "s_min": float(ev[0]), "n_drop": n_drop,
              "mu": mu, "ortho_err": ortho, "w_rows": w_rows}
    return c_b, report


# ------------------------------------------------------------------ the battery


def _mol(basis, atom="H 0 0 0", spin=0):
    return gto.M(atom=atom, basis=basis, spin=spin, verbose=0)


def anchor1_identity():
    """B = A: the projection is the identity, exactly."""
    mol = _mol("ccpvdz", atom="C 0 0 0", spin=2)
    mf = scf.ROHF(mol).run(conv_tol=1e-12)
    c_a = mf.mo_coeff
    c_b, rep = project_orbitals(mol, mol, c_a)
    err = float(np.abs(c_b - c_a).max())
    assert err < 1e-10, f"identity error {err:.2e}"
    assert np.allclose(rep["capture"], 1.0, atol=1e-10)
    return f"PASS  identity      |dC|max = {err:.2e}, capture = 1"


def anchor2_superset():
    """B >= span(A) element-wise: the projected columns reproduce A's
    EXACTLY (A's primitives are literally among B's).  The check is
    ordering-independent: x @ c_a and s @ c_b are the same functions'
    pairings against B's AOs."""
    basis_a = [(0, (8.0, 1.0)), (0, (1.5, 1.0)), (0, (0.3, 1.0))]
    basis_b = basis_a + [(0, (40.0, 1.0)), (0, (0.08, 1.0))]
    mol_a = _mol({"He": basis_a}, atom="He 0 0 0")
    mol_b = _mol({"He": basis_b}, atom="He 0 0 0")
    s_a = mol_a.intor("int1e_ovlp")
    rng = np.random.default_rng(7)
    q = rng.standard_normal((mol_a.nao, 3))
    m = q.T @ s_a @ q
    mu, u = np.linalg.eigh(0.5 * (m + m.T))
    c_a = q @ u @ np.diag(mu ** -0.5)
    c_b, rep = project_orbitals(mol_b, mol_a, c_a)
    x = gto.intor_cross("int1e_ovlp", mol_a, mol_b).T
    s = mol_b.intor("int1e_ovlp")
    err = float(np.abs(x @ c_a - s @ c_b).max())
    assert err < 1e-10, f"superset error {err:.2e}"
    assert np.allclose(rep["capture"], 1.0, atol=1e-10)
    return f"PASS  superset      |pairing err|max = {err:.2e}, capture = 1"


def anchor3_handchecks():
    """Real small bases: the occupied MOs of a coarse basis projected
    into a fine one -- capture high, density honest."""
    out = []
    for atom, spin in (("H", 1), ("He", 0), ("C", 2)):
        mol_a = _mol("sto3g", atom=f"{atom} 0 0 0", spin=spin)
        mol_b = _mol("ccpvdz", atom=f"{atom} 0 0 0", spin=spin)
        mf = scf.ROHF(mol_a).run(conv_tol=1e-12)
        nocc = int(sum(mf.mo_occ > 0))
        c_b, rep = project_orbitals(mol_b, mol_a,
                                    mf.mo_coeff[:, :nocc])
        ok = rep["capture"].min() > 0.95 and rep["ortho_err"] < 1e-10
        assert ok, (atom, rep["capture"], rep["ortho_err"])
        out.append(f"{atom}:{rep['capture'].min():.4f}")
    return "PASS  H/He/C        capture_min " + " ".join(out)


def anchor4_near_singular():
    """The stress: B carries a near-duplicate primitive (exponents 1.0
    and 1+1e-12; the overlap deficit sits below double epsilon, so the
    computed S is numerically singular there).  MEASURED failure shape:
    the NAIVE solve assigns the duplicated pair's weight to an
    ARBITRARY member -- rounding decides which (the member chaos of
    fact 16) -- while the TRUNCATED solve (thresh=1e-8) returns the
    stable, symmetric membership.  The conditioning's contribution,
    demonstrated on the twins' own disease class."""
    basis_a = [(0, (5.0, 1.0)), (0, (1.0, 1.0))]

    def run(delta):
        mol_a = _mol({"He": basis_a}, atom="He 0 0 0")
        mol_b = _mol({"He": basis_a + [(0, (1.0 + delta, 1.0))]},
                     atom="He 0 0 0")
        rng = np.random.default_rng(11)
        s_a = mol_a.intor("int1e_ovlp")
        q = rng.standard_normal((mol_a.nao, 2))
        m = q.T @ s_a @ q
        mu, u = np.linalg.eigh(0.5 * (m + m.T))
        c_a = q @ u @ np.diag(mu ** -0.5)
        _, naive = project_orbitals(mol_b, mol_a, c_a)
        _, trunc = project_orbitals(mol_b, mol_a, c_a, thresh=1e-8)
        return naive["w_rows"], trunc["w_rows"], trunc

    n1, t1, tr = run(1e-12)
    n2, t2, _ = run(1e-8)
    assert tr["n_drop"] == 1, tr["n_drop"]
    assert tr["ortho_err"] < 1e-10
    assert np.allclose(tr["capture"], 1.0, atol=1e-8), tr["capture"]
    # truncated: the duplicated pair (rows 1,2) shares the membership
    # symmetrically, and STABLY across delta
    pair_t = t1[1:3]
    assert abs(pair_t[0] - pair_t[1]) < 0.05 * pair_t.max(), t1
    assert np.allclose(t1, t2, rtol=1e-6), (t1, t2)
    # naive: the membership is arbitrary AND unstable (the rounding
    # flip between deltas is the demonstration)
    asym1 = abs(n1[1] - n1[2]) > 0.5 * max(n1[1], n1[2])
    asym2 = abs(n2[1] - n2[2]) > 0.5 * max(n2[1], n2[2])
    flip = abs(n1[1] - n1[2]) > 0.5 * max(n1[1], n1[2]) and \
        np.argmax([n1[1], n1[2]]) != np.argmax([n2[1], n2[2]])
    assert asym1 and asym2, (n1, n2)
    return (f"PASS  near-singular dropped=1, truncated capture=1, "
            f"membership=({t1[1]:.3f},{t1[2]:.3f}) stable; naive "
            f"arbitrary {tuple(round(x, 2) for x in n1[1:3])} vs "
            f"{tuple(round(x, 2) for x in n2[1:3])}"
            + (" (rounding flip shown)" if flip else ""))


if __name__ == "__main__":
    print(project_orbitals.__doc__.splitlines()[0])
    for fn in (anchor1_identity, anchor2_superset, anchor3_handchecks,
               anchor4_near_singular):
        print(fn())
    print("battery 1-4: ALL PASS")
