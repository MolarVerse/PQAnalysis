.. _developerGuide:

Development
===========

PQAnalysis uses a ``dev`` integration branch and releases from ``main``.
Feature and fix pull requests normally target ``dev``; release pull requests
merge ``dev`` into ``main``.

Local setup
-----------

Clone the repository and install editable development, test and documentation
dependencies:

.. code-block:: console

   $ git clone https://github.com/MolarVerse/PQAnalysis.git
   $ cd PQAnalysis
   $ python -m venv .venv
   $ source .venv/bin/activate
   $ python -m pip install -e ".[dev,test,docs]"

Keep changes focused and add tests at the same ownership boundary as the
behavior being changed.

Tests
-----

The full test script runs the suite with runtime type checking enabled and
again with release settings:

.. code-block:: console

   $ bash pytest.sh

For a focused iteration, pass ordinary pytest arguments:

.. code-block:: console

   $ bash pytest.sh tests/analysis/rdf -q

Documentation
-------------

Build the complete documentation with warnings treated as errors:

.. code-block:: console

   $ python -m sphinx -W --keep-going \
       -b html docs/source docs/build/html

Check internal and external links separately:

.. code-block:: console

   $ python -m sphinx -W --keep-going \
       -b linkcheck docs/source docs/build/linkcheck

The API reference is generated from package modules when Sphinx starts. Do not
hand-edit generated files under ``docs/source/code`` unless the generator or
its templates are being changed. User-facing scientific conventions belong in
the maintained analysis, data and reference pages.

Documentation structure
-----------------------

* ``getting-started.rst`` provides the shortest working path.
* ``analyses/`` explains physical definitions, inputs and interpretation.
* ``_plots/`` contains executable Matplotlib figures built from documented
  analytic models or versioned validation fixtures.
* ``data/`` covers file grammar, trajectories, selections and conversion.
* ``reference/`` indexes CLI and Python interfaces.
* ``userGuide/analysisOutputFiles.rst`` is the canonical output-schema source.
* ``code/`` is generated API material.

Every analysis guide should state the physical quantity, assumptions, units,
minimal input, output fields and interpretation limits. Keep duplicated option
tables in generated API documentation rather than copying them into several
manual pages. Figure captions must identify their data source and distinguish
analytic schematics, validation fixtures and physical benchmark results.

Pull requests
-------------

Pull requests should be reviewer-readable and use a Conventional Commits title,
for example ``feat: add a new analysis command`` or
``fix(io): handle missing trajectory data``. The repository validates the PR
title and uses it as the squash-merge commit message.

The optional local commit-message hook provides earlier feedback:

.. code-block:: console

   $ git config core.hooksPath .githooks

Before requesting review, run the focused tests for the change and every
relevant strict documentation build. CI publishes documentation only from
``main``; pull requests and ``dev`` pushes build it without deploying.

Docstrings
----------

Public Python interfaces use NumPy-style docstrings. Document parameters,
returns, raised exceptions, units and array shapes precisely. Documentation
coverage can be inspected with:

.. code-block:: console

   $ docstr-coverage PQAnalysis
