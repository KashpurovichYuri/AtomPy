# Ticket 23: U SA wave — orbital relaxation vs window growth

Status: OPEN (write-only; execute after ticket 21 part 2; pairs with 22)

Origin: the settled U plan (2026-10-21). The campaign's open physics
question, and one of ticket 22's comparisons in disguise.

## The question

The window-convergence series (1f: f1d3 -> 2f: f2d2 -> 3f: 5L° odd
ground) established that RADIAL CONTENT decides the configuration
ordering. SA-CASSCF optimizes window orbitals freely — can relaxation
recover the missing radial character INSIDE a smaller window and reach
the converged ordering without window growth? If yes, the expensive
40k-detector sectors become unnecessary for identification; if no,
window growth is irreducible and the scan protocol is mandatory.

## The runs (maximal Ms=2 sector, (4,0))

1. **17-root SA over the 5L° M_L set** (L=8: M_L = -8..+8; weights
   1/17 each) in the 2f and 3f windows. Read: does the ensemble window
   borrow the missing f-radial character (occupancies/weights per
   rung); does the per-root table (our CASCI pass in the SA window)
   reproduce the 3f ordering from the 2f window?
2. **relax / pur-relax rungs** (frozen window, core relaxes; pur
   substrate): TM precedent says machine-exact closure expected —
   these are the honest symmetric rungs for U.
3. **Free single-root CASSCF** only as an escape probe (chaos: branch
   on window weights, never pin; the 7s/6p<->window swap class from
   the TM campaign is the likely attractor on U).

## Acceptance

- The relaxation-vs-growth verdict, measured (ordering comparison
  across windows x {frozen CI, SA, relaxed}).
- Per-rung fingerprints in the ticket + results/u_14e; chaos branches
  documented; §3 entries for friction.
- Tests: only deterministic objects (frozen-window rung closure on
  pur; the SA-window CASCI cross-checks) — optimizer reruns are not
  pinned (fact 10 discipline).
