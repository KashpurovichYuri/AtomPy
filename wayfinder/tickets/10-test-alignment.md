# Ticket 10: Test-suite alignment with src (mechanical)

Status: resolved

## Question

The test tree grew organically (9 files) while src grew to 9 modules with
new public helpers. Can the suite be aligned 1:1 with src — filenames and
test names naming the pinned invariant — without touching assertions?

## Resolution (to implement)

Mapping (`test_<module>.py` <-> `src/atomterm/<module>.py`):

- `test_ladder_si.py` MERGED into `test_ladder.py`: one ladder fixture per
  atom (C, Si), shared pattern tests parameterized by case, per-atom
  specifics as case data (chaotic bands: <L^2> observed 0.99..1.92).
- NEW `test_drivers.py` (drivers had no file): `shell_weights` partitions
  unity and groups by angular letter; `select_valence_ao` concatenation
  order + unknown-pattern error; `build_full_mo` block properties
  (S-orthonormal, (nao,nao) shape, core span preserved, window l-purity
  survives the core-overlap correction -- on Si, where the overlap is
  ~0.3); frac-occ ordering assumption (Si: closed block [Ne]+3s, 2/3 on
  the 3p block).
- NEW `test_analysis.py`: `term_symbol`/`quantum_number` over the term
  ladder; `configuration_summary` labels on a known CI vector.
- NEW `test_system.py`: `build_mol` rejects off-origin atoms; label format.
- `test_operators/pipeline/atoms/detbasis/secondquant/angular.py` stay.

Naming convention: `test_<subject>_<invariant>` — physics tests keep
physics names (test_hamiltonian_matches_pyscf_fci style).

## Acceptance

- Pure refactor + the new pre-item tests; every existing assertion
  preserved (same physics, possibly re-homed). Suite green.
- AGENTS.md Tests section updated.

## Resolution (closed)

Done as specified, 51 tests green (stable across repeated runs):

- `test_ladder_si.py` merged into `test_ladder.py`: `LADDER_CASES` (C, Si)
  with per-atom chaotic bands as case data; shared fingerprint tests
  parameterized; sa_pyscf covered where present.
- New `test_drivers.py`: `build_full_mo` block properties (S-orthonormal,
  (nao,nao), core span preserved via projector agreement, window l-purity
  on BOTH C and Si — the ticket-09 regression guard); `shell_weights`
  partitions unity with angular grouping; `select_valence_ao` order +
  error; frac-occ ordering assumption (Si: [Ne]+3s closed, 2/3 on 3p).
- New `test_analysis.py` (`term_symbol`/`quantum_number` ladder moved here
  from test_operators + `configuration_summary` labels incl. s-admixture);
  new `test_system.py` (off-origin rejection, label format).
- Three test-writing findings, fixed in place: the frac-occ cache path
  returned an mf WITHOUT the ensemble `get_occ` hook (hook now set before
  the cache branch — object identity preserved on cache hits); the avg
  band needed `<=` (carbon's averaged commutator is exactly 0.0);
  configuration labels threshold vs amplitude-squared weights (0.1
  admixture = 0.01 weight).

## Post-review correction (ticket 15, step 1)

The spec-axis review found the merge had silently UNIFIED per-atom bands
to the loosest atom, violating "every existing assertion preserved":
carbon's SA comm 1e-9 -> 1e-7, relax comm 1e-9 -> 1e-3, relax <L^2>
abs 1e-6 -> 1e-4, window purity 0.999 -> 0.95, same-basin 1e-4 -> 1e-2,
and the SA-roots count was dropped. All restored as per-case band data
(`exact_comm`, `relax_comm`, `relax_l2_abs`, `relax_purity`, `sa_comm`,
`sa_l2`, `sa_nroots`, `basin_tol`) in the merged parameterized tests.
