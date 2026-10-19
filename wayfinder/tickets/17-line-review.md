# Ticket 17: drivers.py line-by-line review — docstring diet, Löwdin dedup

Status: OPEN (findings approved in principle; execution pending)
Origin: user request 2026-10-19 ("deep line-by-line review … inflated
docstrings … unneeded if-else … shouldn't we reduce the code?") after the
ReferenceSpec deepening (commit 9a22ab0). All findings below were read
line-by-line from drivers.py @ 944 lines (post-9a22ab0 state).

## Verdicts already decided (do NOT relitigate)

1. pipeline vs ladder both stay. They are two deliverables over one
   engine: `pipeline.run_atom` = the multiplet census (CASCI in EVERY Ms
   sector → term tables, results/*_report); `ladder.run_ladder` = the
   treatment experiment (primary sector, 17 rungs as DATA in
   `_registry`, results/*_ladder). They already share every deep module
   (run_rohf/run_casci/constrained_casscf/stage_facts) and, since
   9a22ab0, the sector geometry (system.sector_occupations/ncore_of).
2. The run_casscf_* adapters stay as named functions in drivers (called
   only by ladder._registry, but they encode chkfile-key + mask +
   semantics in one documented place — the AGENTS module table is the
   contract). Diet their docstrings, do not collapse them into data.
3. "Many CASCI variants" is an appearance: ONE run_casci + ONE deep
   constrained_casscf(mo_full=, mask=, weights=); everything else is
   adapters/masks/registry entries.

## Findings to execute

A. `__all__` (drivers.py:21) is stale: add "ReferenceSpec".

B. Docstring diet — move physics narrative to its README/fact home,
   keep the operational contract (semantics, units, gotchas, pointers).
   Targets (current → goal lines):
   - run_rohf ref-param block ~40 → ~15 (the Ti basin census lives in
     fact 21a + ReferenceSpec docstring + README §5.5)
   - run_frac_occ_scf ~33 → ~20 (drop the scalar-core essay ¶1 → README
     §3.2/§2.4; keep control mechanism + open_idx semantics)
   - purified_s_core ~30 → ~12 (4-step causal chain → fact 9 pointer;
     keep the two broken preconditions + the return contract)
   - purified_core_orbitals ~45 → ~28 (drop the closure-theorem essay →
     fact 20 + README §5.3; KEEP construction steps 1–3 and both TRAPS)
   - radial_basis ~22 → ~14 (gauge rationale → fact 16e/README §3.3;
     keep the member-preserving contract + the 8e-2/3e-6 cutoff numbers)
   - run_casscf_core_relax + run_casscf_virtual_mix ~35 → ~16 combined
     (physics notes → RUNGS.md/README §3.3; keep mask reference + the
     max|dE| peculiarity note)
   Net: roughly −120 lines of narrative; ZERO code-behavior change.

C. Löwdin dedup (Duplicated Code, 4 sites): the idiom
   `ev,u = eigh(c.T @ s @ c); c = c @ (u/sqrt(ev)) @ u.T` appears in
   purified_s_core (:430), purified_core_orbitals (:535),
   _core_window_block (:590), radial_basis (:658). Extract ONE helper
   `_lowdin(cols, s)` (private, top of module); four call sites.

D. `_run`'s nroots==1 branch (:794): the energies branch is defensive —
   for nroots=1, `_root_energies(mc)` reproduces `float(mc.e_tot)` by
   construction. Unify to `_root_energies` (KEEP the `ci = [mc.ci] if
   nroots == 1` branch — that one is real, mc.ci is a bare vector).
   Guard: assert equality on one CASCI case in the fast tier is NOT
   needed — the slow tier's max|dE| ~1e-13 pins already cover it.

E. Stale wording: run_casscf_virtual_mix docstring "rung4 (core-relax)"
   → drop the historic rung numbering ("the relax rung").

F. Optional (only if trivial): the two chk-cached SCF factories
   (run_rohf fresh_mf / run_frac_occ_scf inline) share the
   load-or-build pattern; a `_cached_scf(mf_factory, chkfile)` helper
   is defensible but LOW value — skip unless it falls out cleanly.

## Acceptance

- `pytest -m fast` green; full suite green (drivers is physics-bearing).
- Line count: drivers.py ≤ ~830 (from 944) with zero interface change
  (public names, signatures, chkfile keys, __all__ + ReferenceSpec).
- No docstring loses: a number, a trap, a gotcha, or a pointer.
- AGENTS module-table row for drivers unchanged except nothing removed.

## Explicitly out of scope

- Renaming public names; moving adapters into the registry; merging
  pipeline/ladder; touching operators/stages/scans docstrings (they were
  not line-read this pass — a second pass can follow if wanted).
