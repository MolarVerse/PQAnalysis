.. _developerGuide:

Development
===========

PQAnalysis integrates feature and fix branches on ``dev`` and releases from
``main``. The pages below define package ownership, extension steps,
validation requirements and release operations.

Extension path
--------------

.. list-table:: Analysis implementation path
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
     - Arguments and dispatch without duplicate computation
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

Install development, test and documentation dependencies in an isolated
environment:

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

Run Pylint against the package and retain a score above the CI threshold of
9.75:

.. code-block:: console

   $ python -m pylint PQAnalysis --persistent n

Public Python interfaces use NumPy-style docstrings. Document parameters,
returns, raised exceptions, units and array shapes. Inspect coverage with:

.. code-block:: console

   $ docstr-coverage PQAnalysis

Documentation
-------------

Build HTML and check external links with warnings treated as errors:

.. code-block:: console

   $ python -m sphinx -E -W --keep-going \
       -b html docs/source docs/build/html
   $ python -m sphinx -E -W --keep-going \
       -b linkcheck docs/source docs/build/linkcheck

The API reference is generated from package modules when Sphinx starts. Do not
hand-edit generated files under ``docs/source/code``. Add public callables to
:doc:`../reference/functions`, and keep implementation guidance in this
development section.

Executable figures under ``docs/source/_plots`` must be deterministic. Captions
must identify analytical models, versioned validation fixtures and physical
benchmark results correctly.

Pull requests
-------------

Feature and fix pull requests normally target ``dev``. Release pull requests
merge ``dev`` into ``main``. Use a Conventional Commits title, such as
``feat: add a new analysis command`` or
``fix(io): handle missing trajectory data``; CI validates the title. Keep each
commit scoped and reviewable because multi-commit pull requests may retain
their individual commits.

Enable the optional local commit-message hook with:

.. code-block:: console

   $ git config core.hooksPath .githooks

Before review, run the focused tests for every modified ownership boundary and
the relevant strict documentation builds. Pull requests and ``dev`` pushes
build documentation without deploying it; deployment occurs from ``main``.
