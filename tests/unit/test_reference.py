"""Unit pins of atomterm.reference: the ReferenceSpec contract and the
character-forced occupation hook (numpy-only given a Mole -- no SCF
runs here; the SCF-bearing reference tests live in the integration
tier)."""

import numpy as np
import pytest
from pyscf import gto

from atomterm.reference import ReferenceSpec, _forced_rohf_occ
from atomterm.windows import lm_sectors

pytestmark = pytest.mark.fast


def test_reference_spec_validation():
    """ReferenceSpec owns the occ-pattern contract (fact 21a): the two
    pattern lists come as a pair, and the basin gate needs them."""
    ReferenceSpec()  # plain reference: all-None is the default run
    ReferenceSpec(occ_closed=("4s",), occ_open=("3d",))
    with pytest.raises(ValueError, match="together"):
        ReferenceSpec(occ_closed=("4s",))
    with pytest.raises(ValueError, match="needs the"):
        ReferenceSpec(occ_e_ref=-848.40648875)


def test_forced_occ_l_specifiers_and_open_premium():
    """The occupation hook's span vocabulary and ranking: 'l0'..'l9'
    select whole angular-momentum subspaces (the honest vocabulary on
    uncontracted bases), and the open block is ranked by 2 w_open +
    w_closed -- with whole-l spans a 1:1 sum lets closed-character
    leftovers tie the open MOs (measured on U: DIIS then closed s/p
    garbage, E > +1e3), while the premium keeps the early-cycle leakage
    case (Sc/Ti) findable."""
    mol = gto.Mole()
    mol.atom = "C 0 0 0"
    mol.basis = gto.parse("""
        C S
             8.0  1.0
        C S
             1.0  1.0
        C S
             0.1  1.0
        C P
             2.0  1.0
        C P
             0.2  1.0
        C D
             0.5  1.0
    """)
    mol.spin = 2   # nelec = 6: closed = 2, open = 2
    mol.verbose = 0
    mol.build()

    rows_l = {l: sorted(r for (ll, m), rws in lm_sectors(mol).items()
                        if ll == l for r in rws) for l in range(4)}
    hook = _forced_rohf_occ(mol, ["l0", "l1"], ["l2"])
    mo = np.eye(mol.nao_nr())
    occ = hook(np.zeros(mol.nao_nr()), mo)
    closed = set(np.flatnonzero(occ == 2.0))
    open_ = set(np.flatnonzero(occ == 1.0))
    assert all(i in rows_l[0] + rows_l[1] for i in closed)  # s/p-pure
    assert open_ <= set(rows_l[2])                          # d wins the
    # 2:1 premium even though closed leftovers tie at 1.0
    with pytest.raises(ValueError, match="matched no AO labels"):
        _forced_rohf_occ(mol, ["4z"], ["l2"])
