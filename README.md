<img src="https://raw.githubusercontent.com/MolarVerse/PQAnalysis/main/docs/source/logo/PQAnalysis.png" width="250">

# PQAnalysis

[![CI](https://github.com/MolarVerse/PQAnalysis/actions/workflows/ci.yml/badge.svg)](https://github.com/MolarVerse/PQAnalysis/actions/workflows/ci.yml)
[![Docs](https://github.com/MolarVerse/PQAnalysis/actions/workflows/docs.yml/badge.svg)](https://MolarVerse.github.io/PQAnalysis/)
[![codecov](https://codecov.io/gh/MolarVerse/PQAnalysis/graph/badge.svg?token=IDFK8L6IIQ)](https://codecov.io/gh/MolarVerse/PQAnalysis)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The main purpose of this package is to provide useful tools for the analysis of the Molecular Dynamics software package [PQ](https://github.com/MolarVerse/PQ). Furthermore, the intent of this package is to enable straightforward implementations of newly developed analysis tools on top of the provided API.

The future development of this package focuses on two main goals. On the one hand the enhancement of the provided analysis tools and extending its API to be compatible with many other different Molecular Dynamics engines. As this project is only a *hobby* project of the maintainers, any contributions considering enhancement or bug fixes are highly welcomed.

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
