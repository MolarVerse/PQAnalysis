# Developer Changelog

Build, CI, test, refactor, and implementation changes are documented here.
User-visible release notes are documented in [CHANGELOG.md](CHANGELOG.md).
This split starts with the next release; older release entries remain unchanged
in `CHANGELOG.md`.

## Next Release

<!-- insertion marker -->
## [v1.8.0](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.8.0) - 2026-11-13

### CI

- Add a Benchmark Smoke workflow (every benchmark runs once, untimed) and a Performance workflow that times base and head of a pull request on one runner and fails on a much slower benchmark and keeps the table in a sticky PR comment (scripts/perf_compare.py).
- Run kernel-dependent benchmarks with compiled and pure Python kernels (kernel_impl fixture) and make the trajectory reader benchmarks independent of the working directory.
- Add changelog fragments in changes/, a fragment check on pull requests to dev, a Create Release PR workflow and DEV-CHANGELOG.md; releases now stamp both changelogs from the fragments instead of running git-changelog after the tag push.
- Run the Test Coverage workflow only when package, test, example fixture or packaging files change, so documentation, changelog and other workflow edits no longer start it.

### Tests

- Make the add_molecules integration tests assert their results: the framework and box must match the reference, the placed molecule is checked for names, rigid shape and the distance cutoff, and the stale reference restart file was regenerated.
- Add benchmarks for the center of mass: the group routine (500, 2000 and 8000 molecules, with the extent check), AtomicSystem.center_of_mass_residues, and traj_to_com_traj for groups of three atoms, by residue (also with a selection) and for all atoms on seeded water boxes of several sizes and frame counts.
- Add RDF benchmarks for the centers of mass of molecules (file fast path and in-memory, compiled and fallback kernel) and a shared placeholder_atom_names helper for the center of mass tools.
- Add an integration test of traj_to_com_traj on the UMCM-9 example and on a 2000-molecule TIP3P water box (from the PQ example h2o_mm-mshake): centers of mass are compared with an independent calculation, written and read back, and the group-size error is checked.

### Internal

- Add center_of_mass_groups (PQAnalysis/atomic_system/_center_of_mass.py), an array function that returns the center of mass of every group of atoms, and use it for AtomicSystem.center_of_mass, center_of_mass_residues and compute_com_atomic_system; unit tests compare it with the previous per-group formula.
