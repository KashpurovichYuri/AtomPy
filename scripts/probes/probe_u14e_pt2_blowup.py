"""RETIRED RECORD (2026-11-27) -- the PT2 "blow-up" probe.  Its
conclusion is VOID; the code is in git history (the last full version
is at commit 330b807^, ``scripts/probes/probe_u14e_pt2_blowup.py``).

What it claimed: that the large E(2) on the pinned@pur SA rows is a
property of the optimizer's gauge -- "the Sr class carries a
misplaced valence orbital", a "pinned-gauge pathology" visible at
the Eh scale.

Why that is void: the probe priced hand-assembled objects in which
the CI vectors of the CONVERGED window were placed on the START
window's columns (``make_blocks(c_read=...)``) -- a non-orthonormal,
inconsistent orbital set.  Its "case 2" also wrote the table's
``sa-casscf (pinned@pur) (4,12)`` row that way (withdrawn at
330b807).  A consistent object prices the same wavefunction at the
free route's value.  The route layer no longer can build the
inconsistent pair (``vv.routes``: a state is born from one solve and
certified).

The record of what ran: ``logs/probe_u14e_pt2_blowup.log`` (kept --
it is the evidence of the defect, not of a mechanism).
The successor: ``probe_u14e_pt2_autopsy.py`` -- the four-object
autopsy that measures what each inconsistency does.
"""
