.. _developerGuide:

Development
===========

PQAnalysis integrates changes on ``dev`` and releases from ``main``. The guides
define package boundaries, implementation contracts, validation evidence and
release operations.

Extension path
--------------

.. list-table::
   :class: pq-record-table pq-extension-table
   :header-rows: 1
   :widths: 24 38 38

   * - Stage
     - Primary location
     - Contract
   * - Scientific method
     - ``PQAnalysis/analysis/<name>/``
     - Estimator, normalization, units and result shape
   * - Python interface
     - ``PQAnalysis/analysis/<name>/api.py``
     - Validated orchestration shared with the CLI
   * - Command line
     - ``PQAnalysis/cli/<name>.py``
     - Arguments and dispatch, without duplicate computation
   * - Scientific output
     - ``PQAnalysis/analysis/_output_schemas.py``
     - Stable fields, symbols, units and plot projection
   * - Evidence
     - ``tests/analysis/<name>/`` and ``tests/data/<name>/``
     - Analytical, independent, parity and end-to-end tests

.. toctree::
   :maxdepth: 1

   architecture
   adding-analysis
   validation
   release

Local environment
-----------------

Install the package with development, test and documentation dependencies in an
isolated environment:

.. code-block:: console

   $ git clone https://github.com/MolarVerse/PQAnalysis.git
   $ cd PQAnalysis
   $ python -m venv .venv
   $ source .venv/bin/activate
   $ python -m pip install -e ".[dev,test,docs]"

Quality gates
-------------

``pytest.sh`` runs the suite with debug runtime type checking and repeats it
with release settings:

.. code-block:: console

   $ bash pytest.sh
   $ bash pytest.sh tests/analysis/rdf -q

Run pylint against the package and retain a score above the CI threshold of
9.75:

.. code-block:: console

   $ python -m pylint PQAnalysis --persistent n

Public Python interfaces use NumPy-style docstrings. Document parameters,
returns, raised exceptions, units and array shapes. Inspect coverage with:

.. code-block:: console

   $ docstr-coverage PQAnalysis

Documentation
-------------

Build the complete documentation and check links with warnings treated as
errors:

.. code-block:: console

   $ python -m sphinx -E -W --keep-going \
       -b html docs/source docs/build/html
   $ python -m sphinx -E -W --keep-going \
       -b linkcheck docs/source docs/build/linkcheck

The API reference is generated from package modules when Sphinx starts. Do not
hand-edit generated files under ``docs/source/code``. Add public callables to
:doc:`../reference/functions`, and put implementation-level guidance in this
development section.

Executable figures under ``docs/source/_plots`` must be deterministic. Captions
must distinguish analytic schematics, versioned validation fixtures and
physical benchmark results.

Pull requests
-------------

Feature and fix pull requests normally target ``dev``. Release pull requests
merge ``dev`` into ``main``. Use a Conventional Commits PR title, such as
``feat: add a new analysis command`` or
``fix(io): handle missing trajectory data``; the title becomes the squash-merge
commit message.

Enable the optional local commit-message hook with:

.. code-block:: console

   $ git config core.hooksPath .githooks

Before requesting review, run the focused tests for the modified ownership
boundary and every relevant strict documentation build. Pull requests and
``dev`` pushes build documentation without deploying it; deployment occurs
from ``main``.

Performance validation
----------------------

File-backed VACF, MSD, RDF and momentum analyses use bounded compiled fast
paths. A batch path must preserve the numeric operation order of its streaming
fallback and must return to that fallback when the configured memory limit is
exceeded. Parallel work is restricted to independent lag ranges, frames or
private integer histograms; floating-point reductions within one legacy result
must not be reordered.

Install the benchmark dependency and run the focused benchmark suite with:

.. code-block:: console

   $ pip install -e ".[test,benchmark]"
   $ pytest -c benchmarks/pytest.ini benchmarks --benchmark-only

Store a baseline with ``--benchmark-json=baseline.json`` and compare a changed
branch with ``--benchmark-compare=baseline.json``. Runtime assertions do not
belong in CI because host load is variable. Every optimization must instead
pass the compiled and fallback tests plus the relevant fixed-bit legacy oracle
before its benchmark result is considered.
