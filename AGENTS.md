# AGENTS.md — Contributor Guide for PQAnalysis

This file documents the codebase structure, development workflows, and conventions that contributors and coding assistants must follow when contributing to PQAnalysis. The full developer documentation lives in `docs/source/developerGuide/`.

---

## Project Overview

PQAnalysis is a **Python (>= 3.12)** package for analysing [PQ](https://github.com/MolarVerse/PQ) molecular dynamics simulations (RDF, MSD, VACF, momentum, vibrations, format conversion). It ships a `pqanalysis` CLI plus a Python API, uses **Cython** kernels with NumPy/Python fallbacks, and enforces types at runtime with **beartype**. Licensed under MIT. Versions are derived from Git tags by `setuptools_scm`; never edit the generated `PQAnalysis/_version.py`.

---

## Git Workflow — CRITICAL RULES

> **Contributors and coding assistants must follow these rules without exception.**

`main` is the release branch and `dev` is the integration branch.

1. **Never push directly to `main` or `dev`.**
2. **Always create a branch from an up-to-date `dev`** for any change: `git fetch origin && git checkout -b <type>/<short-description> origin/dev`.
3. **Open pull requests targeting `dev`**, not `main`. Only release PRs go from `dev` to `main`.
4. **Never amend or force-push** shared branches; add new commits instead.
5. Use `git push -u origin <branch-name>` when pushing a new branch.
6. Do not push or move release tags (see [Releases](#releases)).

```bash
git fetch origin
git checkout -b feat/my-feature origin/dev
# ... make changes, commit ...
git push -u origin feat/my-feature
# Then open a PR targeting dev
```

---

## Commit and PR Title Convention

PRs are squash-merged and **the PR title becomes the commit message**, so it must follow Conventional Commits. This is enforced by `.github/workflows/pr-title.yml` with the pattern:

```
^(feat|fix|revert|ref|refactor|perf|build|chore|ci|deps|doc|docs|style|test|tests)(\([a-z0-9._-]+\))?!?: .+
```

Examples: `feat: add a new analysis command`, `fix(io): handle missing trajectory data`, `refactor: simplify trajectory parsing`. Use the same style for commit subjects.

PRs use `.github/pull_request_template.md` (a `## Summary` section).

---

## Repository Structure

```
PQAnalysis/
├── .github/workflows/     # ci.yml, pylint.yml, docs.yml, pr-title.yml, release.yml
├── .github/scripts/       # parse_pylint.py (pylint score/comment helper)
├── changes/               # Unreleased changelog fragments (user/, developer/)
├── scripts/               # Changelog tooling and its tests (scripts/tests)
├── PQAnalysis/            # The package (see below)
├── tests/                 # pytest suite, mirrors the package layout; tests/data holds fixtures
├── benchmarks/            # pytest-benchmark suite (own pytest.ini); CI smoke-runs it and compares base vs head timings
├── docs/source/           # Sphinx documentation (user guide, developer guide, reference)
├── examples/              # Tutorial fixtures (e.g. examples/water) used by the docs
├── conda-recipes/         # conda packaging
├── pyproject.toml         # Metadata, dependencies, entry points, setuptools_scm
├── pytest.sh              # Runs the suite twice (beartype DEBUG, then RELEASE)
├── pytest.ini             # Markers, doctest options, coverage options
├── .pylintrc              # Pylint configuration
├── .style.yapf            # Formatter configuration
└── .docstr.yaml           # Docstring-coverage configuration
```

### Package boundaries (`PQAnalysis/`)

| Path | Responsibility |
| --- | --- |
| `analysis/` | Scientific estimators, spectra, result models, analysis-table output |
| `cli/` | Argument definitions and dispatch to public API functions |
| `io/` | Simulation-file readers, writers, format conversion |
| `traj/` | Trajectory containers, engine formats, trajectory-wide checks |
| `atomic_system/`, `core/` | Atomic coordinates, cells, atoms, residues |
| `topology/` | Selections, molecular identity, bonded topology |
| `tools/`, `utils/`, `physical_data/` | Helper tools, utilities, physical constants |

An input-file analysis flows: CLI class → public `api.py` function → input reader → trajectory/structure reader → analysis object and kernel → `AnalysisTable` with an `AnalysisSchema` → writers. **The CLI must not contain a second implementation of the scientific method**; it parses arguments and calls the public API.

Adding an analysis: follow `docs/source/developerGuide/adding-analysis.rst` (`PQAnalysis/analysis/<name>/`, `api.py`, `PQAnalysis/cli/<name>.py`, `analysis/_output_schemas.py`, tests in `tests/analysis/<name>/` and `tests/data/<name>/`).

---

## Development Environment

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev,test,docs]"      # add ",benchmark" for benchmarks
```

Optional extras: `test` (pytest, pylint, coverage, docstr-coverage, ase), `dev` (yapf), `docs` (Sphinx toolchain), `benchmark` (pytest-benchmark).

---

## Build, Test and Quality Gates

```bash
bash pytest.sh                       # full suite: beartype DEBUG, then RELEASE
bash pytest.sh tests/analysis/rdf -q # focused run (args are forwarded to pytest)
python -m pylint PQAnalysis --persistent n
docstr-coverage PQAnalysis
```

- **Run both runtime-type-checking configurations** (`pytest.sh` does) before requesting review; `PQANALYSIS_BEARTYPE_LEVEL` is `DEBUG` then `RELEASE`.
- `pytest.ini` collects `tests` and `PQAnalysis` with `--doctest-modules`, so docstring examples are tests.
- Tests use markers (`atomicSystem`, `core`, `topology`, `traj`, `io`, `analysis`, `utils`, `integration`). Tests mirror the package layout under `tests/`.
- Public interfaces use **NumPy-style docstrings** documenting parameters, returns, raised exceptions, units and array shapes.
- Format with `yapf` using `.style.yapf` (79-column limit).

### Numerical kernels

RDF, MSD, VACF and momentum use compiled Cython kernels (`*.pyx`) with NumPy/Python fallbacks. Both implementations must keep the same signature, normalization and edge-case behavior, and a kernel change needs tests for both plus a direct parity test. Keep file parsing and logging out of kernels. Optimized batch paths must preserve the numeric operation order of the streaming fallback and fall back when the memory limit is exceeded. Do not reorder floating-point reductions.

Benchmarks (the `Benchmark Smoke` workflow runs them once with `--benchmark-disable`; the `Performance` workflow compares base and head timings, see `perf_compare.py`). Kernel-dependent benchmarks use the `kernel_impl` fixture and run as `[compiled]` and `[fallback]`; give new ones seeded, fixed-size data and no cwd-relative paths:

```bash
pip install -e ".[test,benchmark]"
pytest -c benchmarks/pytest.ini benchmarks --benchmark-only
```

### Public contracts

Treat as compatibility surfaces: non-underscored names imported by a package `__init__.py`, CLI command names/arguments/input-file keys, result attributes and array shapes, output field identifiers/units/column order, accepted file formats, and exception types for invalid input. Underscore-prefixed names are internal. Changing a public contract requires tests, documentation, and backward compatibility or an explicit deprecation path.

---

## CI/CD Workflows (`.github/workflows/`)

| Workflow | Trigger | Purpose |
| --- | --- | --- |
| `ci.yml` (Test Coverage) | push to `main`/`dev`, all PRs | `pip install -e ".[test]"`, `bash pytest.sh`, upload coverage to Codecov (fails CI on upload error) |
| `pylint.yml` | push to `main`/`dev`, all PRs | Pylint on `PQAnalysis`; **fails if the score is not above 9.75** |
| `docs.yml` | push and PRs to `main`/`dev` | Strict Sphinx build (`-W --keep-going`); deploys to GitHub Pages only from `main` |
| `pr-title.yml` | PR opened/edited/synchronized | Validates the Conventional Commits PR title |
| `release.yml` | any pushed tag | Builds the sdist, publishes to PyPI, creates the signed GitHub release with notes taken from the prepared changelogs |
| `benchmarks.yml` | push to `main`/`dev` and PRs touching the package or benchmarks | Runs every benchmark once without timing so the suite cannot rot |
| `perf.yml` | PRs touching the package or benchmarks, manual | Times base and head on one runner (`scripts/perf_compare.py`), fails if a benchmark is more than 1.5x slower; maintainers can set `perf-regression-accepted` |
| `changelog.yml` | PRs to `dev` | Runs the changelog tooling tests and requires a valid changelog fragment (or the `skip-changelog` label) |
| `create-release-pr.yml` | Manual | Creates `release/vX.Y.Z` from `dev`, stamps the changelogs, opens the release PR to `main` |
| `linked-issue.yml` | PRs, manual | Requires a linked issue (`Closes #N`, `Fixes #N`, `Resolves #N` or `Related to #N`); release PRs are exempt |

Python 3.12 is used throughout CI.

---

## Documentation

Sphinx sources are in `docs/source/`. Build with warnings as errors:

```bash
python -m sphinx -E -W --keep-going -b html docs/source docs/build/html
python -m sphinx -E -W --keep-going -b linkcheck docs/source docs/build/linkcheck
```

- API reference is generated from package modules at Sphinx start; **do not hand-edit generated files under `docs/source/code`**.
- Add public callables to `docs/source/reference/functions.rst`; put implementation guidance in `docs/source/developerGuide/`.
- Executable figures in `docs/source/_plots` must be deterministic.
- User-facing behavior changes need matching updates in `docs/source/analyses/`, `data/`, `reference/` or `userGuide/`.

---

## Changelog

- **`CHANGELOG.md`**: user-facing release notes.
- **`DEV-CHANGELOG.md`**: build, CI, test, refactor and internal changes.
- **`changes/`**: unreleased changelog fragments; see `changes/README.md`.

**Regular PRs never edit `CHANGELOG.md` or `DEV-CHANGELOG.md` directly.** Every regular PR instead adds or updates a fragment:

```
changes/user/<category>.<title>.md
changes/developer/<category>.<title>.md
```

- Use `user` when an installed-PQAnalysis user would notice the change (behavior, results, inputs/outputs, errors, compatibility, runtime performance); use `developer` for build tooling, CI, tests, refactors and internal maintenance. A PR may add fragments for either or both audiences.
- `<category>` is one of the fixed categories in `changes/README.md` (user: `enhancement`, `change`, `bugfix`, `performance`, `compatibility`, `documentation`; developer: `enhancement`, `bugfix`, `performance`, `build`, `ci`, `test`, `internal`, `documentation`). `<title>` is a lowercase-with-hyphens slug that does not appear in the changelog.
- If an unreleased fragment already covers the change, append a bullet to it instead of creating a new file. A fragment holds Markdown bullets only, one per line, each at most 360 characters.
- A PR that needs no entry carries the `skip-changelog` label instead; coding assistants add a fragment and do not set the label themselves.
- Do not delete fragments or the `.gitkeep` placeholders.

`.github/workflows/changelog.yml` enforces this with `scripts/check_changelog_fragment.py`. The release PR (`create-release-pr.yml`) runs `scripts/update_changelog.py`, which moves every fragment into the matching changelog and deletes the consumed fragments; that is the only point at which the changelog files change. Tests for the tooling live in `scripts/tests/` (`python3 -m unittest discover -s scripts/tests -p 'test_*.py'`).

---

## Releases

Releases are cut by maintainers only; coding assistants must not do this.

1. Merge `main` into `dev` if `main` has commits that `dev` lacks (the workflow in step 2 refuses otherwise).
2. Run the **Create Release PR** workflow with the version (`vX.Y.Z`). It creates `release/vX.Y.Z` from `dev`, stamps `CHANGELOG.md` and `DEV-CHANGELOG.md` from the fragments and opens a PR to `main`.
3. Run `bash pytest.sh` and the strict docs build + link check, then merge the release PR.
4. Tag the merge commit on `main` with `vMAJOR.MINOR.PATCH` (annotated tag).

**A tag push is the publication action**: it publishes to PyPI, which is immutable. Never create or push tags. Correct a defective release with a new patch release. Details: `docs/source/developerGuide/release.rst`.

---

## Contribution Checklist

Before submitting any change:

- [ ] Branch created from `origin/dev` (not `main`); PR targets `dev`
- [ ] PR title follows Conventional Commits (`type(scope): description`)
- [ ] `bash pytest.sh` passes (both beartype levels), including focused tests for the modified module
- [ ] Tests added/updated under `tests/` mirroring the changed module; kernel changes test both Cython and fallback plus parity
- [ ] A performance slowdown is intentional (label `perf-regression-accepted`) or fixed; `Performance` check is green
- [ ] `python -m pylint PQAnalysis` score stays above 9.75
- [ ] New/changed public code has NumPy-style docstrings (`docstr-coverage PQAnalysis`)
- [ ] Documentation updated; strict Sphinx build passes (`-W`)
- [ ] CLI contains no duplicated analysis logic; public contracts preserved or deprecated
- [ ] A changelog fragment under `changes/user/` and/or `changes/developer/` (or the `skip-changelog` label set by a maintainer)
- [ ] No edits to `CHANGELOG.md`, `DEV-CHANGELOG.md`, `PQAnalysis/_version.py` or generated `docs/source/code`
- [ ] No tags created or pushed
