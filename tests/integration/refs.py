"""Campaign reference data: the pinned free-ion term tables and the
gated ReferenceSpec constants. A plain module (NOT conftest) so test
files can import it by name at any time -- a ``conftest`` basename is
claimed per-directory by pytest and direct imports of it are
order-dependent across tiers."""

from atomterm.reference import ReferenceSpec


# Sc CAS(3,6)/{3d,4s}, Ms=1/2 sector: the complete free-ion term table of
# 3d1 4s2 + 3d2 4s1 + 3d3 over the 90 determinants (ticket 11).
SC_EXACT_TERMS = {"2D": 20, "4F": 14, "2F": 14, "4P": 6, "2G": 18,
                  "2P": 6, "2S": 1, "2H": 11}

# Ti CAS(4,6)/{3d,4s}, Ms=1 sector: the complete term table of 3d2 4s2 +
# 3d3 4s1 + 3d4 over the 120 determinants (ticket 12). S >= 1 only (the
# sector has no singlet component); dimension check: d2s2 gives 10 =
# 3F x7 + 3P x3, d3s1 gives 60, d4 gives 50.
TI_EXACT_TERMS = {"3D": 15, "3F": 35, "3G": 18, "3H": 22, "3P": 15,
                  "5D": 5, "5F": 7, "5P": 3}


# Ti occupations are forced by character (fact 21 / ticket 12): the default
# energy-ordered ROHF is an aufbau runaway (+123.07 mEh). [Ar]+4s closed,
# 3d open. occ_e_ref gates the SCF basin: with the same forced occupations
# a metastable +11.2 mEh member-gauge basin exists whose polarized core
# smears the input window at the 1e-2 level.
TI_REF = ReferenceSpec(
    occ_closed=("1s", "2s", "3s", "4s", "2p", "3p"),
    occ_open=("3d",),
    occ_e_ref=-848.40648875,
)
