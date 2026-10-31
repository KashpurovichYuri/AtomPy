"""The SO(3) -> D2h dictionary: which Abelian irreps an atomic term's
components occupy, and which irrep every determinant of a sector
carries.

Two routes to L²-symmetric states meet at this dictionary and certify
each other in the tests: ours (diagonalize the symmetry-blind sector,
MEASURE symmetry with the exact operators) and PySCF's (CONSTRUCT
symmetry -- one ``wfnsym`` solver per D2h component). The routes are
compared, never merged; this module is the shared vocabulary both
sides speak.

Contracts
---------
- **D2h is the working dialect.** PySCF's atom ``symmetry="SO3"``
  mode is Abelianized there (m-sum blocks) -- its CI machinery cannot
  pin L; D2h pins the |M_L| classes properly.
- **Irrep ids are PySCF's D2h XOR encoding** (``Ag=0, B1g=1, B2g=2,
  B3g=3, Au=4, B1u=5, B2u=6, B3u=7``): the irrep of a product of
  objects is the XOR of their ids, so a determinant's irrep is the
  XOR over its occupied spatial orbitals and a doubly-occupied pair
  cancels to Ag.
- ``branch_d2h`` multiplicities sum to ``2L + 1``; an irrep may host
  more than one component (for L >= 2) -- the multiplicity is the
  state-average weight and the per-solver root count.
- The determinant ordering has ONE home (``detbasis``):
  ``det_irreps`` is aligned with ``detbasis.determinant_basis``.

Layering: character algebra beneath ``branch_d2h`` (the tables are
derived from the reflection-flip structure, never transcribed), the
XOR product beneath ``det_irreps``, and the PySCF symmetry boundary
in ``orb_irreps`` alone.
"""

from __future__ import annotations

import numpy as np
from pyscf.symm import label_orb_symm

from atomterm.detbasis import determinant_basis

__all__ = [
    "IRREP_IDS",
    "IRREP_NAMES",
    "branch_d2h",
    "det_irreps",
    "orb_irreps",
]


# ---------------------------------------------------------------- tables

#: The eight D2h irrep names, indexed by PySCF's id (the XOR encoding).
IRREP_NAMES = ("Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u")

#: Name -> PySCF's D2h irrep id (Z2^3 bits: irrep product = id XOR).
IRREP_IDS = {name: i for i, name in enumerate(IRREP_NAMES)}

# Each D2h operation flips a fixed set of Cartesian coordinates; each
# irrep is the sign triple (s_x, s_y, s_z) of a monomial basis function
# (z -> B1u, yz -> B3g, xyz -> Au, ...). The character of an irrep at
# an operation is the product of the signs of the coordinates the
# operation flips; the operation order matches PySCF's character table.
_FLIPS = {
    "E": (), "C2x": (1, 2), "C2y": (0, 2), "C2z": (0, 1),
    "i": (0, 1, 2), "sx": (0,), "sy": (1,), "sz": (2,),
}
_SIGNS = {  # id -> the monomial basis function's sign triple
    0: (1, 1, 1), 1: (-1, -1, 1), 2: (-1, 1, -1), 3: (1, -1, -1),
    4: (-1, -1, -1), 5: (1, 1, -1), 6: (1, -1, 1), 7: (-1, 1, 1),
}
_OPS = ("E", "C2x", "C2y", "C2z", "i", "sx", "sy", "sz")

#: (8 irreps, 8 operations) characters, rows in id order -- derived.
_CHARACTERS = np.array([
    [(-1) ** sum(1 for c in _FLIPS[op] if signs[c] < 0) for op in _OPS]
    for signs in _SIGNS.values()
])


# ------------------------------------------------------------- branching

def branch_d2h(L: int, parity: str) -> dict[str, int]:
    """The D2h irreps occupied by one term's ``2L + 1`` components.

    The SO(3) -> D2h subduction by character projection: a term of
    total angular momentum ``L`` and inversion parity ``parity``
    reduces to a multiset of D2h irreps. Equal-energy components of
    the free atom distribute over these irreps -- one ``wfnsym``
    solver per irrep, the multiplicity as its root count, realizes
    the term's complete component set.

    Parameters
    ----------
    L : the term's total orbital angular momentum (integer >= 0).
    parity : ``"g"`` or ``"u"`` -- the term's inversion parity
        (for a one-electron shell, ``(-1)**l``).

    Returns
    -------
    Multiplicity per irrep NAME, in id order, omitting zeros; the
    values sum to ``2L + 1``.

    Raises
    ------
    ValueError : if ``L`` is negative or ``parity`` is not g/u.

    Examples
    --------
    >>> branch_d2h(1, "u")
    {'B1u': 1, 'B2u': 1, 'B3u': 1}
    >>> branch_d2h(2, "g")
    {'Ag': 2, 'B1g': 1, 'B2g': 1, 'B3g': 1}
    """
    if L < 0:
        raise ValueError(f"L must be >= 0 (got {L})")
    if parity not in ("g", "u"):
        raise ValueError(f"parity must be 'g' or 'u' (got {parity!r})")
    p = 1 if parity == "g" else -1
    # characters of the (L, parity) representation over the 8 classes:
    # rotations by pi contribute (-1)**L; reflections (i * C2) also
    # carry the inversion parity.
    chi = np.array(
        [2 * L + 1, (-1) ** L, (-1) ** L, (-1) ** L,
         p * (2 * L + 1), p * (-1) ** L, p * (-1) ** L, p * (-1) ** L]
    )
    counts = _CHARACTERS @ chi / 8
    if np.any(np.abs(counts - np.round(counts)) > 1e-9):
        raise AssertionError("subduction produced fractional counts")
    return {IRREP_NAMES[i]: int(round(c))
            for i, c in enumerate(counts) if round(c)}


# ------------------------------------------------------- determinants

def det_irreps(orbsyms, na: int, nb: int) -> np.ndarray:
    """The D2h irrep id of every determinant of the sector.

    A determinant carries the irrep product of its singly-occupied
    spatial orbitals; a doubly-occupied orbital contributes the
    identity (XOR of a pair cancels). The result identifies the
    symmetry blocks of any sector matrix built in detbasis order --
    the block a ``wfnsym`` solver addresses.

    Parameters
    ----------
    orbsyms : the D2h irrep ID per active spatial orbital (PySCF
        ids; see :func:`orb_irreps`).
    na, nb : the sector's alpha and beta electron counts.

    Returns
    -------
    Integer id per determinant, aligned with
    :func:`atomterm.detbasis.determinant_basis` order.

    Raises
    ------
    ValueError : if any id falls outside 0..7.

    Examples
    --------
    The p-shell triplet sector -- three determinants, the three
    |M_L| classes of the one L=1 term (gerade: u x u = g):

    >>> det_irreps([5, 6, 7], 2, 0).tolist()
    [3, 2, 1]
    """
    ids = np.asarray(orbsyms, dtype=int)
    if ids.ndim != 1 or np.any(ids < 0) or np.any(ids > 7):
        raise ValueError(f"orbsyms must be D2h ids in 0..7 (got {orbsyms!r})")
    norb = len(ids)
    out = np.array(
        [np.bitwise_xor.reduce([ids[i % norb] for i in det], axis=None)
         for det in determinant_basis(norb, na, nb)],
        dtype=int,
    )
    return out


# ------------------------------------------------------------ boundary

def orb_irreps(mol, mo_coeff) -> np.ndarray:
    """The D2h irrep id per column of an MO block (the PySCF boundary).

    Wraps :func:`pyscf.symm.label_orb_symm` and converts to PySCF's
    numeric ids -- the form :func:`det_irreps` and ``wfnsym``
    solvers consume. This is the module's only PySCF symmetry call.

    Parameters
    ----------
    mol : a molecule built with ``symmetry="D2h"``.
    mo_coeff : the MO columns to label (AO rows x columns).

    Returns
    -------
    Integer id per column.

    Raises
    ------
    ValueError : if the molecule's group is not D2h.

    Examples
    --------
    >>> mol = gto.M(atom="C 0 0 0", basis="sto-3g",
    ...             symmetry="D2h", verbose=0)
    >>> orb_irreps(mol, np.eye(mol.nao_nr())).tolist()  # doctest: +SKIP
    [0, 0, 5, 6, 7]
    """
    if mol.groupname != "D2h":
        raise ValueError(
            f"expected a D2h molecule (got group {mol.groupname!r}) -- "
            "the dictionary's ids are PySCF's D2h encoding"
        )
    return np.asarray(
        label_orb_symm(mol, mol.irrep_id, mol.symm_orb, mo_coeff)
    )
