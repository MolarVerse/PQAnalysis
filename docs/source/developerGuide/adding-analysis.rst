Adding an Analysis
==================

Use an existing complete analysis, such as RDF or MSD, as the structural
reference. Keep the scientific estimator independent from its CLI and file
format adapters.

1. Define the scientific contract
---------------------------------

State the observable, normalization, units, periodic-boundary treatment,
selection semantics and returned array shapes before implementing the
calculation. Decide which behavior reproduces a legacy tool and which behavior
is a corrected or newly defined method.

The analysis guide and tests must use the same definitions. A numerical result
without its normalization and units is not a complete interface.

2. Implement the analysis package
---------------------------------

A file-driven analysis typically owns these modules:

.. code-block:: text

   PQAnalysis/analysis/<name>/
      __init__.py
      api.py
      <name>.py
      <name>_input_file_reader.py
      <name>_output_file_writer.py
      exceptions.py

The analysis class or numerical function owns the calculation. The input reader
validates configuration and the writer serializes results. Avoid importing CLI
code from the analysis package.

If a compiled kernel is required, provide a Python or NumPy fallback with the
same callable signature. Import the compiled implementation first and fall back
only when it is unavailable, following the RDF, MSD and VACF packages.

3. Define inputs and outputs
----------------------------

Add input keys to the analysis input reader with explicit types, defaults and
validation. Required files should be validated before the trajectory is
processed. Keep aliases only when they preserve an established input contract.

Define output columns in ``PQAnalysis/analysis/_output_schemas.py`` using
:class:`~PQAnalysis.analysis.output.AnalysisColumn` and
:class:`~PQAnalysis.analysis.output.AnalysisSchema`. Field identifiers are ASCII
programmatic names; symbols and units may use Unicode scientific notation.

The data writer should subclass
:class:`~PQAnalysis.analysis.output.AnalysisDataWriter`, create an
:class:`~PQAnalysis.analysis.output.AnalysisTable` from the numerical columns,
and delegate CSV, TSV and XVG exports to the common writer. Preserve a legacy
native row format only when compatibility requires it.

4. Add the public API
---------------------

The function in ``api.py`` is the shared orchestration layer. It should:

1. read and validate the analysis input;
2. construct trajectory, structure or Hessian readers;
3. construct every output writer so path conflicts fail early;
4. instantiate and run the scientific analysis;
5. write the result and return useful in-memory data where appropriate.

Export the function and supported result types from the analysis package
``__init__.py`` and, for a general analysis workflow, from
``PQAnalysis.analysis``. Add the callable to :doc:`../reference/functions`.

5. Add the CLI
--------------

Implement a ``CLIBase`` subclass in ``PQAnalysis/cli/<name>.py``. Its
``add_arguments`` method defines only command-line parsing; ``run`` calls the
public API function. Reuse common arguments from
``PQAnalysis/cli/_argument_parser.py``, including repeatable ``--export`` for
analysis tables.

Register the class in the dispatch dictionary in ``PQAnalysis/cli/main.py``.
Add a ``[project.scripts]`` entry in ``pyproject.toml`` only when a standalone
executable is part of the supported interface.

6. Add evidence and documentation
---------------------------------

The minimum complete change includes:

* analytical or independently computed numerical tests;
* legacy parity tests when compatibility is claimed;
* compiled-kernel and fallback parity where both exist;
* input-reader validation and default tests;
* API and CLI end-to-end tests;
* native output and CSV, TSV and XVG tests;
* existing-file failure tests for every output path;
* an analysis page defining equations, assumptions, inputs and interpretation;
* entries in the function index, command reference and output-schema page.

Follow :doc:`validation` for reference-data provenance and tolerance rules.
