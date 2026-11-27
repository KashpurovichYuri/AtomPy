"""The legs' assertion bands (declared constants)."""

L_L2 = 72.0  # the 5L parent term: L = 8 -> L(L+1)
S_S2 = 6.0   # spin 4 -> S = 2 -> S(S+1)

# the assertion bands (declared constants; each cites the observed
# fresh values -- logs/u14e_vv_*.log):
BAND_L2 = 1e-5     # the Casimir meters are gauge-invariant; the
                   # seeded readout's observed dev 3.0e-06 (free)
                   # and 1.9e-07 (pinned, all 72.000000)
BAND_E2_MEH = 2.0  # a fresh SA draw is a new orbital gauge and E(2)
                   # is state-specific; observed delta 0.0006-0.0029
                   # mEh (the engine wobble is ~0.3 uEh run-to-run)
BAND_ROUTE = 1e-6  # Eh; the ENSEMBLE energy is the route comparison
                   # surface (variational, gauge-invariant; E(2) is
                   # orbital-sensitive and printed, never asserted);
                   # observed 3.3e-09 Eh
BAND_CLASS = 0.01  # the certified-class meters (the (6,13) pinned
                   # components: 16 at L(L+1)=90 + one at 72)
