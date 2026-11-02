# Ticket 33 — The symmetry-breaking anatomy deliverable

## Goal

Deliver the frozen-landscape anatomy (WHY symmetry breaks: the flat
and descending directions, the JT-like drive, the ±pair blind spots)
as a designed TABLES deliverable — replacing the retired pointer
sections (`scan.md` curves nobody reads) with tables that show the
physics.

## Status: DESIGN-FIRST (nothing runs until the shape is agreed)

The candidate shape (mockup, per atom, one Production-or-Readings
section):

1. The m-symmetric θ-grid: | θ | ΔE | ‖[H,Lz]‖ | — the closed-window
   radial curve (rises on p-shell valence, the d-shell inversion is
   the TM record);
2. The member-tilt ladder at the scan's optimum: the descending
   subsets with ‖[H,Lz]‖ side by side — the ±pair blind spots VISIBLE
   (subsets that descend with [H,Lz] = 0 are exactly the directions
   the sector commutator cannot see — the projector pair's reason);
3. Reading keys in README · Results (one home), the same tag-line
   contract as every energy table.

Open design questions for the session:

- full tilt map (31 subsets) vs the descent ladder + blind pairs only
  (~10 rows);
- one section per atom, or one cross-atom READING in the record (the
  per-atom tables stay atom-specific — the standing rule);
- where it sits in the four-type grammar (Readings is the natural
  home; it carries no variational energies).

## Inputs

- `scans.py` (the module stays: `ScanSpec` per case, `radial_scan`,
  `member_tilt_map`) — everything regenerable; nothing imports the
  deleted artifacts;
- the retired artifacts' content as the requirements list (the
  θ-grids with one/two-body decomposition; the tilt maps with the
  per-subset commutators; the tilt-follow curves from the
  sa_decomposition era);
- the delivered tables lose their pointers in the same commit
  (scan.md files deleted at the C/Si pass of ticket 32's step 8).

## Verification

The scan machinery is deterministic and already pinned; this ticket's
pins are the RENDERER's (byte-determinism, the declared-subset
policy) — the same pattern as the lane/lattice emitters.
