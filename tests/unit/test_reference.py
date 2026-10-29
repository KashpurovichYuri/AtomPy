"""Unit pins of atomterm.reference: the ReferenceSpec contract and the
character-forced occupation hook (numpy-only given a Mole -- no SCF
runs here; the SCF-bearing reference tests live in the integration
tier)."""

import numpy as np
import pytest

from atomterm.reference import ReferenceSpec, _forced_rohf_occ
from atomterm.windows import lm_sectors
from helpers import uncontracted_carbon_mole

pytestmark = pytest.mark.fast


def test_reference_spec_validation():
    """The occ-pattern contract: the two pattern lists come as a pair,
    and the basin gate needs them."""
    ReferenceSpec()  # plain reference: all-None is the default run
    ReferenceSpec(occ_closed=("4s",), occ_open=("3d",))
    with pytest.raises(ValueError, match="together"):
        ReferenceSpec(occ_closed=("4s",))
    with pytest.raises(ValueError, match="needs the"):
        ReferenceSpec(occ_e_ref=-848.40648875)


def test_forced_occ_l_specifiers_select_whole_subspaces():
    """'l0'..'l9' force occupations onto whole angular-momentum
    subspaces (the honest vocabulary on uncontracted bases)."""
    mol = uncontracted_carbon_mole()  # nelec = 6: closed = 2, open = 2
    rows_l = {l: sorted(r for (ll, m), rws in lm_sectors(mol).items()
                        if ll == l for r in rws) for l in range(4)}
    hook = _forced_rohf_occ(mol, ["l0", "l1"], ["l2"])
    occ = hook(np.zeros(mol.nao_nr()), np.eye(mol.nao_nr()))
    closed = set(np.flatnonzero(occ == 2.0))
    assert closed <= set(rows_l[0] + rows_l[1])       # s/p-pure block
    with pytest.raises(ValueError, match="matched no AO labels"):
        _forced_rohf_occ(mol, ["4z"], ["l2"])


def test_forced_occ_open_block_carries_a_strict_premium():
    """The open block is ranked by 2 w_open + w_closed: open-shell
    character beats closed leftovers even where they tie."""
    mol = uncontracted_carbon_mole()
    rows_l = {l: sorted(r for (ll, m), rws in lm_sectors(mol).items()
                        if ll == l for r in rws) for l in range(4)}
    hook = _forced_rohf_occ(mol, ["l0", "l1"], ["l2"])
    occ = hook(np.zeros(mol.nao_nr()), np.eye(mol.nao_nr()))
    open_ = set(np.flatnonzero(occ == 1.0))
    # on identity MOs the unchosen s/p rows (w_closed=1, w_open=0) tie
    # the d rows (w_open=1, w_closed=0) under a 1:1 sum; the 2:1
    # premium breaks it toward d. Measured on U: a 1:1 sum let DIIS
    # close s/p garbage (E > +1e3); the closed term still catches
    # early-cycle open character hiding in core-like MOs (Sc/Ti).
    assert open_ <= set(rows_l[2])
