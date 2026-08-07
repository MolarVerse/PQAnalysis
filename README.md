<img src="https://raw.githubusercontent.com/MolarVerse/PQAnalysis/main/docs/source/logo/PQAnalysis.png" width="250">

# PQAnalysis

[![CI](https://github.com/MolarVerse/PQAnalysis/actions/workflows/ci.yml/badge.svg)](https://github.com/MolarVerse/PQAnalysis/actions/workflows/ci.yml)
[![Docs](https://github.com/MolarVerse/PQAnalysis/actions/workflows/docs.yml/badge.svg)](https://MolarVerse.github.io/PQAnalysis/)
[![codecov](https://codecov.io/gh/MolarVerse/PQAnalysis/graph/badge.svg?token=IDFK8L6IIQ)](https://codecov.io/gh/MolarVerse/PQAnalysis)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

PQAnalysis reads structures, trajectories, velocities and Hessians produced by
[PQ](https://github.com/MolarVerse/PQ). Its command-line and Python interfaces
share parsers, numerical kernels and schema-defined outputs for RDF, MSD, VACF,
vibrational, spectral and momentum analyses.

Development focuses on validated analysis methods and support for additional
molecular-dynamics engines. The maintainers develop PQAnalysis in their free
time; focused analysis contributions and bug fixes are welcome.

## Installation

Install with pip:

    pip install pqanalysis

## Development

Clone the repository and install the development, test and documentation
dependencies in an isolated environment:

    git clone https://github.com/MolarVerse/PQAnalysis.git
    cd PQAnalysis
    python -m venv .venv
    source .venv/bin/activate
    python -m pip install -e ".[dev,test,docs]"

Run the test suite with both debug and release runtime type checking:

    bash pytest.sh

The [developer documentation](https://molarverse.github.io/PQAnalysis/developerGuide/developerGuide.html)
covers package architecture, adding an analysis, scientific validation and the
tag-driven release process. The
[function index](https://molarverse.github.io/PQAnalysis/reference/functions.html)
lists the supported Python entry points directly.

Use squash merges for pull requests. The pull request title becomes the commit
message on the target branch, so PR titles must follow
[Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/):

    feat: add a new analysis command
    fix(io): handle missing trajectory data

The local commit hook is optional contributor feedback. Enable it with:

    git config core.hooksPath .githooks
