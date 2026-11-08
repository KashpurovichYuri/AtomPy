# Ticket 34 — the folded reference + the window characterization matrix (Ti, then Sc)

Opened 2026-11-08 (user-approved round; supersedes nothing -- new
work).  One session for the Ti leg; the Sc port is the closing step.

## Why (the problem, measured)

The unc lanes' totals sit +118.4 mEh (Tiunc) / +38.7 (Scunc) above the
free rows, and the deep dive showed why: the fold's confined core on a
wide span IS THE FREE ION (Ti4+ RHF reproduced to 0.016 mEh; the
confined core -845.190591 vs the neutral pur core -845.051227 =
-139.4 mEh), and the valence CI runs on that ion field with a window
cut from the NEUTRAL field's MOs -- the frozen window cannot respond
(the 2x2: fixed-wavefunction field cost +119.7 mEh, CI re-relaxation
recovery 0.02 mEh).  The contracted bases sit in the other regime
(their fold spans are the core's own contractions: +18.6 core / -17.1
valence = +1.5 net) -- one construction, span-width-driven regime
flip.  The inconsistency is prerequisite work for the term-ID round:
term splittings are correlation differences and a mismatched field
distorts them.  Elements/Scunc/logs/d2h_probe.log carries the probes.

## Step 1 -- the Tiunc basin census, then PROMOTE (user-decided)

Cycle ROHF init guesses on the unc basis at conv_tol 1e-12, spin 2,
ALL at the character-forced 3d2 4s2 filling where forcing applies
(the plain aufbau draws are basin-chaotic -- the zoo: -848.39639 plain
shallow / -848.40648875 the current trusted record / -848.40744343 the
deepest known, found by the D2h-seeded draw, open pair degenerate at
-0.090219, reverse-seed-stable on C1):

    for guess in (default, "atom", "huckel", "minao"):
        run ROHF(seed or guess) -> collect converged E
    plus the SEEDED starts: from elements/Ti/rungs/rohf.chk (contracted,
    cross-basis lstsq+Loewdin -- reproduces -848.40648875 exactly, the
    38393f8 probe) and from the D2h draw at -848.40744343.

ADOPT the deepest converged same-filling basin as the trusted
reference: re-store elements/Tiunc/rungs/rohf.chk, declare
ReferenceSpec(occ_closed=("1s","2s","3s","4s","2p","3p"),
occ_open=("3d",), occ_e_ref=<the adopted E>) on the case (through the
shell layer -- the forcing spans resolve via shell_realization),
re-run the Tiunc free cells (minutes at the compressed sizes),
regenerate the table.  If the census certifies -848.40744343 the
anchor and every dE shift ~1 mEh -- expected and fine.

## Step 2 -- the folded-reference utility (reference.py)

Frozen-core ROHF: the core density FIXED at the confined solution
(the lane's fold construction), the valence electrons (4 for Ti)
self-consistent in its field, the character forcing + basin
discipline as every runner.  One home; the substrate signature in,
(mol, mf) out.  Tests: (a) the frozen core block is BIT-UNCHANGED by
the run; (b) the folded total >= the full ROHF total (variational,
both same filling); (c) the recovery measurement: the folded valence
re-contracts on the ion field -- predicted to recover ~+120 mEh of
the lane-vs-free gap (a PREDICTION to measure, not assume; if it
does not reconcile, the decomposition machinery says why).

## Step 3 -- the characterization matrix (NO winner; MO default)

On ONE fold (3d4s@ar) and ONE field (the folded reference), run the
window variants and CERTIFY them side by side -- E, <L2>, ptp,
closure ([P,Lz]/[P,Lperp]), window purity, capture:

(i)   AO-compressed k=1, re-cut from the FOLDED reference's MOs (the
      math unchanged; the signal field is the change), report gains
      the PER-RADIAL CAPTURE SPECTRUM (see below);
(ii)  MO-PICKED m-complete -- THE DEFAULT (user-decided): the whole
      reference MOs whose set covers every m-member (5 d partners +
      the 4s), the span closed by set construction, window_closure
      verifying; the core projection owns the core tails; the
      valence-family tails reported by the purity meters;
(iii) NO-informed: one CI pass -> the term-averaged density ->
      natural.occupation_groups -> the correlation-weighted radials
      (the m-set-safe machinery exists).

The support continuum: (i) and (ii) are the ends of the
allowed-support dial (AO combinations restricted to the family = (i);
unrestricted AO combinations near MOs collapse to (ii); the middle =
wider-span radials via a declared support set -- only if the numbers
ask).  The k decision is data-driven: the per-radial spectrum + the
NO occupations (k=2 for the 4s stays one declaration away,
CAS(4,7)).

The MO default is DECLARED after its meters pass -- if a meter fails,
the default stays AO and this ticket records why.  The AO class
remains in the toolbox and the tables: it is the symmetry-preserving
class (band truncation is m-complete by construction; the MO echo's
<L2>=0.005 failure was the m-incomplete SET, not MOs per se).

## Step 4 -- the lane re-runs + the disclosure

The lane legs re-run on the default variant + the folded substrate
(the stored 3F/2D records stand as the constructions that validated
the multi-column fold; the field caveat documented).  The 3F/2D
EXACTNESS re-verified on the new field.  The totals reconciled; the
tables carry the ZERO-DISCLOSURE note (the lane's fold zero is the
confined core = the free ion on a wide span; the free rows' zero is
the neutral core).

## Step 5 -- the Sc port (T-A-prime)

Scunc needs no census (its trusted record is already the deeper of
its two known basins).  The folded reference + the matrix + the lane
re-run, same acceptance.

## Riding along (same commits)

- The ANISOTROPY COLUMN: render core_anisotropy in the tables (the
  lane rows show the machine floor -- itself the certificate).
- The PER-RADIAL CAPTURE SPECTRUM in window_compressed's report:
  p_ij = |e_j^dag S c_i|^2 against each S-orthonormal radial -- the
  k-decision sensor (bounded fractions; rank by G, report by
  projectors).
- The DOCSTRING REVIEW (review, not append -- user-directed): core.py
  (the m-scalar theorem -- the span Loewdin and the closed-shell Fock
  both block-diagonal with identical m-blocks -> replicated
  eigenvectors; the RHF-for-ion framing, brief); operators.py
  (spherical_core: the field defect is quadrupole-like = totally
  symmetric in D2h -- irreps cannot see it; only the SO(3) l=0
  projection isolates it); windows.py (the Gram is the SELECTION
  object -- member-summed ranking, S-metric necessity; the member
  chains); casscf.py (the state_average_mix_ sentence); the twins'
  TABLE preambles once the numbers are final.
- map.md: the projected-Gram DECLINE record (revisit condition: if
  the per-radial spectrum shows ranking distortions).

## Acceptance

The exact term certificates hold on the consistent construction;
the field-window consistency measured (the recovery number); the
captures/spectra/purity recorded; the contracted comparison DIAGNOSTIC
ONLY (evidence, never the anchor -- user-directed: the picking
procedure is judged by our own meters).

## Operational

One regen at a time (the U14e route batch holds its store until it
exits; the watcher completes that table automatically -- do not
interfere).  The fast tier is the per-commit gate.  NO /tmp paths (an
earlier slip): probe chks are in-repo *.tmp beside their targets,
removed after.  Tables presented pipe-checked.  U14e work (T-C) stays
GATED on this ticket + the Sc port.

## After this ticket

T-B the contraction-recipe review (provenance lookup for PySCF's
Ti/Sc ccpvdz; our radials vs the authors' contraction coefficients --
verified structure: Ti (20s16p8d2f) general contractions, s = 5
columns sharing 19 primitives + 1 free, nao 43); T-C U14e (gated);
T-D docs /code-review (anti-inflation); T-E the term-ID round (NIST
as the acceptance test).  Stop after T-E's review -- no physics waves
before it.
