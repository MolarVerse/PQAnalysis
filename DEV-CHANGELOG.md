# Developer Changelog

Build, CI, test, refactor, and implementation changes are documented here.
User-visible release notes are documented in [CHANGELOG.md](CHANGELOG.md).
This split starts with the next release; older release entries remain unchanged
in `CHANGELOG.md`.

## Next Release

<!-- insertion marker -->
## [v1.10.0](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.10.0) - 2026-11-27

### Enhancements

- New `AtomicSystem.iter_fit_atomic_system` yields each placed molecule; `AddMolecule` writes the placed ones on `AtomicSystemError` or `KeyboardInterrupt` (SIGTERM mapped to it in the main thread), then re-raises. The topology extension uses the placed count (#341).
- Pass a `mode` argument from the rdf, msd and vacf CLIs through `rdf()`, `msd()` and `vacf()` to all their data, log, spectrum and windowed writers.
- Add `TrajectoryFormat.ASE` and `AseTrajReader` (lazy `ase.io.Trajectory`), the `ase` extra, and dispatch file-based readers (H5MD, ASE) in `TrajectoryReader` through one table (#364).
- Add `read_energy_table`, which turns an `EnergyFileReader` result into an `AnalysisTable` (`analysis/_energy_table.py` builds the columns and default plot fields), and route `.en` inputs of `convert_analysis_output` through it, so energy files reuse the existing CSV/TSV/XVG/native writers (#343).
- Energy-table plot presets carry stacked graphs (`AnalysisPlotGraph`, built in `_energy_table.py`); the XVG writer puts them in Grace graphs G0..Gn with styled sets and autoscale, stores them as `# XVG_GRAPH` lines and reads them back. `AnalysisSchema.units_in_header` writes and parses `NAME [unit]` CSV/TSV headers (#343).
- Share crystal geometry preparation between the Python API and CLI, and test restart round trips for orthogonal and skewed cells.
- Add `PQAnalysis/io/energy_file_header.py` (`read_energy_header`); `EnergyFileReader` prefers its labels over the `.info` file, warns when they differ (ignoring the `.info` time unit), applies the step-to-ps conversion only to a `SIMULATION-TIME[step]` column, and checks the header width per file (#343).
- `read_energy_header` maps the header column `STEP[-]` to `SIMULATION-TIME` in `step`, so the step-to-ps conversion and the `.info` comparison stay unchanged.
- Add `TrajectoryFormat.H5MD` and `H5MDReader`, to which `TrajectoryReader` dispatches frames, frame counts and cells; `h5py` is imported lazily, is the `h5md` extra and is part of the `test` extra (#333).
- `PQForceFieldSetup` gets a `periodic` flag and `long_range()`, which pick the method by box; `parm2pq(long_range=...)` and the CLI pass an override (#351).
- Add `PQAnalysis.io.force_field` (`parm2pq`, `load_structure`, `convert_structure`, `PQForceFieldSetup`), the `parm2pq` CLI and the `ff` extra; ParmEd is imported lazily and is part of the `test` extra, and the tests build small ParmEd structures and round-trip them through AMBER and CHARMM files (#342).
- Expose the OpenMM preparation recipe as a reusable API with optional dependencies. (#370)
- `traj2extxyz` reads extxyz input (by extension or a `properties=` comment line), accepts `units=pq`, warns without a units key and raises `PQValueError` for other units; `TrajectoryWriter` appends `units=<profile>` to extxyz comment lines (#324).

### Bug Fixes

- `iter_fit_atomic_system` seeds one random generator per run and passes it to `_fit_atomic_system` instead of reseeding per addition (#349).
- `EnergyFileReader` probes the default `.info` path with `os.path.isfile` instead of constructing a `BaseReader`, which logged an error before the exception was caught (#343).
- The `nep` extxyz profile writes the `properties` value without quotes, because GPUMD splits it at `:` without stripping them (#315).
- The `nep` extxyz profile maps `D` and `T` to `H` and raises `TrajectoryWriterError` for `Q`, `X`, `cav`, `sup` and `dum`, because GPUMD only accepts elements of the periodic table in `nep.in` (#321).
- The converter detects O-H-H bond triangles and writes `mshake.dat` from the bond lengths instead of three SHAKE constraints per water (#356).
- `RDFLogWriter.write_before_run` uses the `Å` string directly instead of UTF-8 bytes (#352).
- Input-file values accept a quoted `STRING` terminal, parse errors raise `InputFileError` with line, column and a quoting hint, and selection ranges check their step (#357).
- `TrajectoryWriter.write` gets an `else` branch that raises `TrajectoryWriterError` for `AUTO` and `H5MD` (#334).

### Performance

- `AtomicSystem.fit_atomic_system` checks overlaps with a periodic cell list (`atomic_system/_overlap_grid.py`), extended per placed molecule, instead of all-pairs `distance` per trial rotation. Same image and norm arithmetic, so a fixed `CONSTANT_SEED` gives the same output (#341).

### Build

- Build and test wheels with cibuildwheel in the new `wheels.yml` workflow, used by `release.yml` and by pull requests that change the build. `scripts/test_wheel.py` checks the wheel modules and commands, then runs the tests outside the source tree. MSVC gets no `-ffp-contract=off`.

### CI

- Match Pylint findings whose limit count or McCabe rating changed with the base finding, so a function that grows from 19 to 20 local variables is not reported as a new finding.
- The Pylint report matches duplicate-code findings by the pair of modules, so a moved or shrunk block is no longer reported as a new finding (#292).

### Tests

- Add benchmarks for the molecule centers, the vibrational analysis and atom selections, with seeded in-memory data.
- Add benchmarks for the trajectory writer and the traj2extxyz, traj2qmcfc and traj2box conversions with seeded in-memory trajectories.
- Add scaling benchmarks of RDF, MSD and VACF over the number of atoms and frames at constant density, so changes of the complexity become visible.
- Add benchmarks for building, slicing and querying a `Topology` and for the shake topology generation with seeded water boxes.
- Add end-to-end tests of rst2xyz, xyz2rst, xyz2gen and gen2xyz: restart columns and moldescriptor inference round trip on a 2000 molecule box, velocities, forces, QMCFC, randomize, gen lattice of a triclinic cell, and the writing modes.
- Add end-to-end tests of traj2box, traj2qmcfc and traj2extxyz on the triclinic examples: box files, VMD vertices, the QMCFC dummy atom and round trip, extended xyz read back with ASE, and the writing modes.
- Add end-to-end tests of `build_nep_traj`, `build_spectrum` and `continue_input` through the command line, checked against independent NumPy and text references.

### Internal

- Remove the stray `scratch_install.log` from the repo root and ignore `*.log` files there (#359).

### Documentation

- Add `docs/source/tutorials/solvated-protein.rst` and its files in `examples/parm2pq` (OpenMM build script, PQ and RDF inputs, tleap/sander inputs); the page includes them with `literalinclude`, and the generated outputs are ignored.

## [v1.9.0](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.9.0) - 2026-11-15

### Bug Fixes

- Build the output writers of `msd()`, `vacf()` and `check_momentum()` before the restart file and trajectory are opened.
- Resolve PQ energy timesteps from native headers or an explicit reader argument, reject invalid and conflicting metadata, and cover chained offsets, legacy step units and QMCFC compatibility.
- The `ase` extxyz profile maps the PQ names `D` and `T` to `H`, and `Q`, `cav`, `sup` and `dum` to `X`; other names are written as before (#305).
- Build the RDF data and log writers right after the input file is read in `rdf()`, before the restart file and trajectory are opened.

### CI

- Make Pylint and performance reports concise, with expandable diagnostics instead of large tables. Compare Pylint with the PR base and publish performance comments through an isolated pq-perf-bot workflow.
- Keep raw Pylint diagnostics out of comments, collapse linked finding summaries and passing gate rules, and use readable benchmark names with plain-language timing changes.
- Match duplicate-code findings by file ranges rather than example excerpts, avoiding false new findings when Pylint chooses a different example.

### Tests

- Add end-to-end tests of the rdf, msd, vacf, check_momentum and vibrations commands on examples/water, checked against independent NumPy calculations, with compiled and fallback kernels and the CSV and XVG exports.

### Internal

- Add MoleculeCenters.mean_velocities (vectorized and chunked over frames) and VACF benchmarks for the centers of mass of molecules (file fast path and in-memory, compiled and fallback kernel).
- Move the molecule centers of the RDF to PQAnalysis/analysis/_molecule_centers.py with a from_topology factory shared with the MSD, and add MSD benchmarks for the centers of mass (file fast path and in-memory, compiled and fallback kernel).
- Move the inference of the topology files of the RDF input file reader into a shared mixin (_TopologyFilesMixin) that the MSD input file reader uses as well.

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
