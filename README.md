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

Clone the PQAnalysis GitHub repository and navigate into the directory:

    git clone https://github.com/MolarVerse/PQAnalysis.git
    cd PQAnalysis

Install with pip:

    pip install .

In order to have a nice changelog strategy, please stick to [conventional commit definiton](https://www.conventionalcommits.org/en/v1.0.0-beta.4/). In order to make git accept only valid commit message please make a symlink of the git-hook template `.githooks/commit-msg`:

    ln -s .githooks/commit-msg .git/hooks/

More information on preferred and undesirable commit message will come soon...
