"""Unit pins of atomterm.windows: the selection/meter layer (pure
numpy: span_projector, select_span, shell_weights, occupied_core_idx,
the primitive chain) and build_full_mo's two representations (AO rows,
column vectors) on a synthetic substrate, plus the STO-3G degenerate
substrates (null core, whole space as window) and the frozen-carbon
selection pins (tests/carbon_reference.py -- no SCF)."""

import numpy as np
import pytest

from pyscf import gto

from atomterm.windows import (NOBLE_CORES, Window, build_full_mo,
                              canonical_order, closed_electrons,
                              complete_orthonormal,
                              fold_shells, lm_sectors, lowdin,
                              occupied_core_idx, primitive_chain,
                              select_primitives, select_span,
                              select_valence_ao, select_shell_ao,
                              select_valence_mo,
                              shell_weights, span_projector)
from helpers import uncontracted_carbon_mole

pytestmark = pytest.mark.fast


# ------------------------------------------------------ the selection layer

def test_span_projector_weights_one_in_span_zero_orthogonal():
    mol = gto.M(atom="C 0 0 0", basis="ccpvdz", spin=2, verbose=0)
    s = mol.intor_symmetric("int1e_ovlp")
    rows = select_valence_ao(mol, ["2p"])
    p = span_projector(s, rows)
    assert np.abs(p - p.T).max() < 1e-12          # symmetric form
    c_in = lowdin(np.eye(mol.nao_nr())[:, rows], s)[:, 0]
    assert abs(float(np.real(c_in @ p @ c_in)) - 1.0) < 1e-10
    c_d = np.eye(mol.nao_nr())[:, select_valence_ao(mol, ["3d"])[0]]
    assert abs(float(c_d @ p @ c_d)) < 1e-12      # different l: exactly 0


def test_select_span_l_vocabulary_and_label_prefix():
    import re

    mol = gto.M(atom="C 0 0 0", basis="ccpvdz", spin=2, verbose=0)
    labels = [line.split()[-1] for line in mol.ao_labels()]
    # spherical p labels carry the member suffix (2px, 2py, 2pz, ...)
    l1 = [i for i, lbl in enumerate(labels) if re.match(r"\d+p", lbl)]
    assert select_span(mol, "l1") == sorted(l1)  # the whole l subspace
    assert select_span(mol, "2p") == select_valence_ao(mol, ["2p"])
    with pytest.raises(ValueError, match="matched no AO labels"):
        select_span(mol, "9x")


def test_select_primitives_rejects_contracted_shells():
    """Exponent-ranked primitive selection is undefined for
    contracted shells (cc-pVDZ carbon p is one 3-exponent
    contraction) -- the primitive vocabulary is for uncontracted
    bases, where shell labels lie about n."""
    mol = gto.M(atom="C 0 0 0", basis="ccpvdz", spin=2, verbose=0)
    with pytest.raises(ValueError, match="contraction slot"):
        select_primitives(mol, 1, count=1)


def test_shell_weights_rows_sum_to_one_for_s_orthonormal_columns():
    mol = gto.M(atom="C 0 0 0", basis="ccpvdz", spin=2, verbose=0)
    s = mol.intor_symmetric("int1e_ovlp")
    c = lowdin(np.eye(mol.nao_nr()), s)  # a complete S-orthonormal set
    names, w = shell_weights(mol, c, s)
    assert names[0] == "s"
    assert np.abs(w.sum(axis=1) - 1.0).max() < 1e-12


class _FakeMf:
    """The synthetic reference: an SCF-free substrate provider (mol +
    S-orthonormal "canonical" MOs + their occupations)."""

    def __init__(self, mol, mo, ncore):
        self.mol, self.mo_coeff = mol, mo
        self.mo_occ = np.zeros(mo.shape[1])
        self.mo_occ[:ncore] = 2.0

    def get_ovlp(self):
        return self.mol.intor_symmetric("int1e_ovlp")


def _substrate():
    """C/cc-pvdz + a synthetic reference: per-l S-orthonormal
    "canonical" MOs and a hand-built exactly-closed pur core (one p
    radial replicated on every member), theta-rotated against each
    other so canonical columns carry cross-span content at a planted
    scale (sin 0.3) -- the physical situation the projections exist
    for."""
    mol = gto.M(atom="C 0 0 0", basis="ccpvdz", spin=2, verbose=0)
    s = mol.intor_symmetric("int1e_ovlp")
    n = mol.nao_nr()
    sec = lm_sectors(mol)
    # NOTE: two independent l-pure blocks would make the projected
    # window Gram STRUCTURALLY singular -- the theta-rotation against
    # the complement is what keeps the substrate non-degenerate
    rows = {l: sorted(r for (ll, _), rs in sec.items() if ll == l
                      for r in rs)
            for l in sorted({ll for ll, _ in sec})}
    blocks = {l: lowdin(np.eye(n)[:, rows[l]], s) for l in rows}

    # the pur core: 2 s radials (s has one member) + one p radial
    # replicated on every member -- E (x) m-set, exactly closed (each
    # member's chain carries the same radial: its first element)
    p_pur = lowdin(np.column_stack(
        [np.eye(n)[sec[(1, m)][0]] for m in range(3)]), s)
    c_pur = np.hstack([blocks[0][:, :2], p_pur])
    ncore = c_pur.shape[1]

    # canonical MOs: [pur | complement] with the p radials theta-rotated
    # against the first complement directions (an S-orthonormal pair
    # stays S-orthonormal under a 2x2 rotation)
    full = complete_orthonormal(c_pur, s)
    theta = 0.3
    c, w = np.cos(theta), np.sin(theta)
    core_can = full[:, :ncore].copy()
    win_can = full[:, ncore:].copy()
    for i in range(3):
        core_can[:, 2 + i], win_can[:, i] = (
            c * full[:, 2 + i] + w * full[:, ncore + i],
            -w * full[:, 2 + i] + c * full[:, ncore + i])
    mo = np.hstack([core_can, win_can])
    return _FakeMf(mol, mo, ncore), mol, s, c_pur, rows


def test_build_full_mo_ao_rows_equal_ao_columns():
    mf, mol, s, c_pur, rows = _substrate()
    ncore = c_pur.shape[1]
    act = select_valence_ao(mol, ["2p"])
    a = build_full_mo(mf, list(range(ncore)), act)
    b = build_full_mo(mf, list(range(ncore)), np.eye(mol.nao_nr())[:, act])
    assert np.array_equal(a, b)


def test_build_full_mo_window_columns_on_pur_core():
    """Canonical columns on the pur core: orthonormal assembly, the
    pur block intact, the planted cross-span leak removed."""
    mf, mol, s, c_pur, rows = _substrate()
    ncore = c_pur.shape[1]
    mo = build_full_mo(mf, list(range(ncore)), mf.mo_coeff[:, ncore:],
                       c_core=c_pur)
    assert mo.shape == (mol.nao_nr(), mol.nao_nr())
    err = np.abs(mo.conj().T @ s @ mo - np.eye(mo.shape[1])).max()
    assert err < 1e-12
    # the core block is the pur core (Loewdin of an S-orthonormal block
    # is itself, to machine precision)
    assert np.allclose(mo[:, :ncore], c_pur, atol=1e-12)
    # the window columns are exactly S-orthogonal to the core span
    w = mo[:, ncore:ncore + 3]
    assert np.abs(w.conj().T @ s @ c_pur).max() < 1e-12
    # and the projection did real work: the canonical window columns
    # carried pur-core-span content at the planted scale
    w_can = mf.mo_coeff[:, ncore:ncore + 3]
    leak = np.linalg.norm(w_can.conj().T @ s @ c_pur)
    assert leak > 1e-3, leak      # the fact-17-class component existed
    assert np.linalg.norm(w.conj().T @ s @ c_pur) < 1e-12  # ...and is gone


def test_window_columns_assemble_through_build_full_mo():
    """Window.columns(c_win) == build_full_mo directly; columns already
    orthogonal to the core pass through unchanged."""
    mf, mol, s, c_pur, rows = _substrate()
    ncore = c_pur.shape[1]
    # an S-orthonormal window orthogonal to the pur core: the d block
    c_win = lowdin(np.eye(mol.nao_nr())[:, rows[2]], s)[:, :3]
    assert np.abs(c_win.conj().T @ s @ c_pur).max() < 1e-12
    win = Window("synthetic no", None, c_win, False)
    mo = win.columns(mf, ncore, c_core=c_pur)
    closed = list(range(ncore))
    assert np.array_equal(mo, build_full_mo(mf, closed, c_win, c_core=c_pur))
    err = np.abs(mo.conj().T @ s @ mo - np.eye(mo.shape[1])).max()
    assert err < 1e-12
    assert np.allclose(mo[:, ncore:ncore + 3], c_win, atol=1e-12)


# ---------------------------------------------------------------------
# degenerate substrates (STO-3G): null core, the whole AO space as the
# window -- rotation-closed by construction, so the certificates are
# exact; any deviation is a bug.

def test_h_sto3g_null_core_full_window():
    """Empty core + full window: the zero-column paths and an
    exactly-zero closure certificate."""
    from pyscf import scf
    from atomterm.operators import window_closure

    mol = gto.M(atom="H 0 0 0", basis="sto-3g", spin=1, verbose=0)
    mf = scf.RHF(mol).run(verbose=0)
    n = mol.nao_nr()
    s = mol.intor_symmetric("int1e_ovlp")
    mo = build_full_mo(mf, [], list(range(n)))
    assert mo.shape == (n, n)
    assert np.abs(mo.conj().T @ s @ mo - np.eye(n)).max() < 1e-12
    cl = window_closure(mol, mo)
    assert cl["lz"] < 1e-14 and cl["lperp"] < 1e-14


def test_c_sto3g_null_core_full_window():
    """Null core, whole AO space as window (CAS(6,5), sector (4,2)):
    certificates exactly zero, [H, Lz]/[H, L2] at machine precision,
    3P ground."""
    from pyscf import scf
    from atomterm.operators import build_operators, window_closure

    mol = gto.M(atom="C 0 0 0", basis="sto-3g", spin=2, verbose=0)
    mf = scf.ROHF(mol).run(verbose=0)
    n = mol.nao_nr()
    mo = build_full_mo(mf, [], list(range(n)))
    cl = window_closure(mol, mo)
    assert cl["lz"] < 1e-13 and cl["lperp"] < 1e-13

    ops = build_operators(mol, mf, mo, 0, 4, 2)
    h, lz, l2, s2 = ops["H"], ops["Lz"], ops["L2"], ops["S2"]
    assert np.linalg.norm(h @ lz - lz @ h) < 1e-11
    assert np.linalg.norm(h @ l2 - l2 @ h) < 1e-11
    ev, evec = np.linalg.eigh(h)
    v = evec[:, 0]
    assert abs(float(np.real(v.conj() @ l2 @ v)) - 2.0) < 1e-8  # 3P ground
    assert abs(float(np.real(v.conj() @ s2 @ v)) - 2.0) < 1e-8


# ------------------------------------------- the frozen-carbon selection pins

def test_select_valence_ao_concatenates_in_pattern_order(carbon):
    _, mol, _, act_idx = carbon
    assert select_valence_ao(mol, ["2p"]) == act_idx
    assert select_valence_ao(mol, ["2s", "2p"]) == \
        select_shell_ao(mol, "2s") + select_shell_ao(mol, "2p")


def test_select_valence_ao_rejects_unknown_pattern(carbon):
    _, mol, _, _ = carbon
    with pytest.raises(ValueError, match="9f"):
        select_valence_ao(mol, ["9f"])


def test_shell_weights_pure_l_for_raw_shell_aos(carbon):
    _, mol, mf, act_idx = carbon
    s = mf.get_ovlp()
    c = np.eye(mol.nao_nr())[:, act_idx]  # raw 2p AOs: exactly l = 1
    names, w = shell_weights(mol, c, s)
    assert names == ["s", "p", "d"]
    assert np.allclose(w[:, names.index("p")], 1.0, atol=1e-12)


# ---------------------- closed-block and primitive-chain pins (synthetic)

def test_occupied_core_idx_selects_the_closed_block_not_canonical_order():
    """The core fold takes the DOUBLY-OCCUPIED block, not canonical
    order; the contracted-basis coincidence is checked, not assumed."""
    import warnings

    class FakeMF:
        def __init__(self, mo_occ):
            self.mo_occ = np.asarray(mo_occ, dtype=float)

    # U-style: canonical [0:5] are virtual junk (eps ~ 0.02) while the
    # closed block sits at 22/61/88/120/143 -- folding range(ncore)
    # corrupted the fold at max|dE| = 14.5 Eh
    occ = [0.0] * 10
    for k in (2, 5, 6, 8, 9):
        occ[k] = 2.0
    mf = FakeMF(occ)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        core = occupied_core_idx(mf, 5)
    assert core == [2, 5, 6, 8, 9]
    assert any("not the first" in str(w.message).lower() for w in rec)

    # contracted-style coincidence: closed block IS the canonical head
    mf = FakeMF([2.0, 2.0, 2.0, 1.0, 0.0])
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        assert occupied_core_idx(mf, 3) == [0, 1, 2]
    assert not any("not the first" in str(w.message).lower() for w in rec)

    # asking for more closed orbitals than exist is refused loudly
    with pytest.raises(ValueError, match="closed block has"):
        occupied_core_idx(mf, 4)


def test_canonical_order_reorders_into_closed_open_virtual():
    """The start='mo' substrate takes the canonical matrix in the RUN
    convention -- closed block | open shell | virtuals, energy-ordered
    within the classes: the raw PySCF order is the m-block listing and
    would fold virtual junk on uncontracted bases. Identity (silent)
    when the canonical order already coincides."""
    import warnings

    class FakeMF:
        def __init__(self, mo_occ, mo_energy):
            self.mo_occ = np.asarray(mo_occ, dtype=float)
            self.mo_energy = np.asarray(mo_energy, dtype=float)

    # U-style scramble: closed at 4, 7; open at 0, 5; virtuals elsewhere
    # with energies that order differently than the indices
    occ = [1.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 2.0, 0.0, 0.0]
    eps = [0.05, 0.30, -1.10, 0.10, 0.02, -0.40, 0.01, -0.90, 0.03, 0.20]
    mf = FakeMF(occ, eps)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        perm = canonical_order(mf, 2)
    assert perm[:2] == [2, 7]                     # the closed block first
    assert perm[2:4] == [5, 0]                    # open, energy-ordered
    assert perm[4:] == [6, 4, 8, 3, 9, 1]         # virtuals, energy-ordered
    assert any("reordering" in str(w.message).lower() for w in rec)

    # contracted-style coincidence: closed head, sorted tail -- silent
    occ = [2.0, 2.0, 1.0, 0.0, 0.0]
    eps = [-1.0, -0.5, -0.2, 0.1, 0.4]
    mf = FakeMF(occ, eps)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        assert canonical_order(mf, 2) == [0, 1, 2, 3, 4]
    assert not any("reordering" in str(w.message).lower() for w in rec)

    # more core orbitals than closed orbitals is refused loudly
    with pytest.raises(ValueError, match="closed block has"):
        canonical_order(mf, 3)


def test_select_valence_mo_ranks_by_named_band_character():
    """The start='moc' window: the SCF's own valence rides (the closed
    tail, then the open shells), the virtual complement is ranked by
    Mulliken character on the named bands -- a higher-energy on-band
    column beats a lower-energy off-band one (different-l rows are
    exactly orthogonal, so an s-band virtual scores ZERO on a p-band
    pattern) -- and the window size is the AO label's own count."""
    mol = gto.M(atom="C 0 0 0", basis="ccpvdz", spin=2, verbose=0)
    nao = mol.nao_nr()
    r_1s = select_shell_ao(mol, "1s")[0]
    r_2s = select_shell_ao(mol, "2s")[0]
    r_2p = select_shell_ao(mol, "2p")
    r_3p = select_shell_ao(mol, "3p")[0]
    r_3s = select_shell_ao(mol, "3s")[0]
    # columns: [1s, 2s] closed core, two 2p opens, then virtuals where
    # the ENERGY order (the 3s row first, eps .05) disagrees with the
    # character order -- the high virtual is a 3p+2p combination whose
    # content sits partly ON the named 2p rows (eps .30); a pure other-
    # band row would score exactly zero (block-restricted metric)
    v_lo = np.zeros(nao); v_lo[r_3s] = 1.0
    v_hi = np.zeros(nao); v_hi[r_3p] = 1.0; v_hi[r_2p[0]] = 0.6
    mo = np.column_stack([np.eye(nao)[:, [r_1s, r_2s, *r_2p[:2]]],
                          v_lo, v_hi])
    occ = [2.0, 2.0, 1.0, 1.0, 0.0, 0.0]
    eps = [-11.0, -2.5, -0.4, -0.3, 0.05, 0.30]

    class FakeMF:
        pass

    mf = FakeMF()
    mf.mol, mf.mo_coeff = mol, mo
    mf.mo_occ, mf.mo_energy = np.asarray(occ), np.asarray(eps)
    mf.get_ovlp = lambda: mol.intor_symmetric("int1e_ovlp")

    ncas = len(select_valence_ao(mol, ["2p"]))
    assert ncas == 3
    perm = select_valence_mo(mf, ["2p"], 2)
    assert sorted(perm) == list(range(mo.shape[1]))   # a permutation
    assert perm[:2] == [0, 1]                        # the core prefix
    assert perm[2:5] == [2, 3, 5]                    # opens, then the
    # ... on-band virtual (col 5, eps .30) over the off-band one
    # (col 4, eps .05) -- the energy order would pick [2, 3, 4]
    assert perm[5] == 4

    # the closed tail rides the window (a promoted pair's columns):
    # ncore=1 leaves the 2s column to the valence block
    perm1 = select_valence_mo(mf, ["2p"], 1)
    assert perm1[:4] == [0, 1, 2, 3]

    # a label smaller than the SCF's own valence is refused loudly
    with pytest.raises(ValueError, match="already exceeds"):
        select_valence_mo(mf, ["2s"], 2)


def test_select_valence_mo_on_frozen_carbon_matches_the_true_argmax(carbon):
    """On the frozen reference the picked virtual is the virtual with
    the largest 2p-row S-projection weight, scored independently of
    the selector (an explicit block quadratic), and the permutation
    covers every column exactly once."""
    _cfg, mol, mf, act_idx = carbon
    s = mf.get_ovlp()
    ss = s[np.ix_(act_idx, act_idx)]
    cb = mf.mo_coeff[act_idx, :]
    score = np.einsum("ij,ij->j", cb, ss @ cb).real
    virt = [i for i in range(mf.mo_coeff.shape[1]) if mf.mo_occ[i] == 0]
    best = max(virt, key=lambda i: score[i])
    perm = select_valence_mo(mf, ["2p"], 2)
    assert sorted(perm) == list(range(mf.mo_coeff.shape[1]))
    assert perm[2 + 2] == best          # the single complement slot


def test_primitive_chain_and_select_primitives():
    """The primitive chain: exponent-ordered m-complete slots; band and
    count select whole primitives; contractions are rejected."""
    from atomterm.operators import window_closure

    mol = uncontracted_carbon_mole()

    # the chain: exponent ordering DIRECT, m-complete slots, row identity
    chain = primitive_chain(mol, 1)
    assert [e for e, _ in chain] == [0.2, 2.0]
    assert chain[0][1] == [6, 7, 8]          # the exp-0.2 p shell
    assert chain[1][1] == [3, 4, 5]          # the exp-2.0 p shell

    # band semantics: exactly the exponents inside the interval
    assert select_primitives(mol, 1, exponent_band=(0.5, 5.0)) == [3, 4, 5]
    assert select_primitives(mol, 1, exponent_band=(0.1, 0.3)) == [6, 7, 8]
    # count = the most diffuse end
    assert select_primitives(mol, 1, count=1) == [6, 7, 8]
    assert select_primitives(mol, 1, count=2) == [3, 4, 5, 6, 7, 8]

    # the error surface: both/neither selector, over-count, empty band
    with pytest.raises(ValueError, match="exactly one"):
        select_primitives(mol, 1)
    with pytest.raises(ValueError, match="exactly one"):
        select_primitives(mol, 1, count=1, exponent_band=(0.0, 9.0))
    with pytest.raises(ValueError, match="chain has"):
        select_primitives(mol, 1, count=3)
    with pytest.raises(ValueError, match="no l=1 primitive"):
        select_primitives(mol, 1, exponent_band=(100.0, 200.0))

    # the physics pin: a selected primitive m-set is exactly closed
    cl = window_closure(mol, np.eye(mol.nao_nr())[
        :, select_primitives(mol, 1, count=1)])
    assert cl["lz"] < 1e-12 and cl["lperp"] < 1e-12

    # contracted shells have no single exponent: clear error, and the
    # message points at the label/MO routes instead
    ccd = gto.Mole()
    ccd.atom = "C 0 0 0"
    ccd.basis = "cc-pvdz"
    ccd.verbose = 0
    ccd.build()
    with pytest.raises(ValueError, match="contractions"):
        primitive_chain(ccd, 1)


# -------------------------------------------------------- the fold vocabulary


def test_fold_shells_rules():
    """The three rule classes resolve: the universal noble-gas table,
    the semi class rule (the noble core ending at np, minus np, plus
    (n+1)s -- the open-d period's), and the explicit list (concatenated
    or comma-separated -- it must consume entirely)."""
    assert fold_shells("he") == ("1s",)
    assert fold_shells("AR") == NOBLE_CORES["ar"]
    assert fold_shells("kr") == ("1s", "2s", "2p", "3s", "3p",
                                 "4s", "3d", "4p")
    assert fold_shells("semi", open_d="3d") == ("1s", "2s", "2p", "3s", "4s")
    assert fold_shells("semi", open_d="4d") == ("1s", "2s", "2p", "3s",
                                                "3p", "4s", "3d", "5s")
    assert fold_shells("1s2s2p4s") == ("1s", "2s", "2p", "4s")
    assert fold_shells("1s, 2s") == ("1s", "2s")


def test_fold_shells_refusals():
    """Unknown names, junk lists, and semi without a usable open-d
    shell all refuse at the resolver (the seam that owns the
    vocabulary)."""
    for name, open_d in (("xe2", None), ("np", None), ("1s2q", None),
                         ("s", None), ("", None),
                         ("semi", None), ("semi", "3p"), ("semi", "7d")):
        with pytest.raises(ValueError):
            fold_shells(name, open_d=open_d)


def test_closed_electrons():
    """Two electrons per m-slot: the noble cores' shell closures and
    the ladder's Ti folds (semi leaves 8 of 22 active; the M-shell
    folds leave 10)."""
    assert closed_electrons(NOBLE_CORES["he"]) == 2
    assert closed_electrons(NOBLE_CORES["ne"]) == 10
    assert closed_electrons(NOBLE_CORES["ar"]) == 18
    assert closed_electrons(NOBLE_CORES["rn"]) == 86
    assert closed_electrons(fold_shells("semi", open_d="3d")) == 14
    assert 22 - closed_electrons(("1s", "2s", "2p", "4s")) == 10
    assert 22 - closed_electrons(("1s", "2s", "2p", "3s")) == 10
