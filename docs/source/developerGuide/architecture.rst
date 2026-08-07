Architecture
============

PQAnalysis separates scientific computation from file orchestration. An
input-file analysis normally follows this path:

.. code-block:: text

   CLI class
      -> public API function
         -> analysis input reader
         -> trajectory, structure or Hessian reader
         -> analysis object and numerical kernel
         -> AnalysisTable with an AnalysisSchema
         -> native writer and optional CSV, TSV or XVG writers

The command line and Python API therefore share the calculation, validation and
output code. A CLI class should parse arguments and call a public API function;
it should not contain a second implementation of the scientific method.

Package boundaries
------------------

.. list-table:: Source ownership
   :class: pq-record-table pq-package-table
   :header-rows: 1
   :widths: 30 70

   * - Path
     - Responsibility
   * - ``PQAnalysis/analysis/``
     - Scientific estimators, spectra, result models and analysis-table output
   * - ``PQAnalysis/cli/``
     - Argument definitions and dispatch to public API functions
   * - ``PQAnalysis/io/``
     - Simulation-file readers, writers and format conversion
   * - ``PQAnalysis/traj/``
     - Trajectory containers, engine formats and trajectory-wide checks
   * - ``PQAnalysis/atomic_system/`` and ``PQAnalysis/core/``
     - Atomic coordinates, cells, atoms and residues
   * - ``PQAnalysis/topology/``
     - Selections, molecular identity and bonded topology

The :doc:`../reference/functions` page lists callable entry points. The
:doc:`../reference/api` page exposes the classes and generated modules behind
them.

Public contracts
----------------

Treat the following as compatibility surfaces:

* non-underscored functions and classes deliberately imported by a package
  ``__init__.py``;
* command names, arguments and input-file keys;
* analysis result attributes and array shapes;
* output field identifiers, symbols, units and column order;
* accepted trajectory and structure formats;
* exception types raised for invalid user input.

Modules, functions and attributes beginning with an underscore are internal.
Changing a public contract requires tests, documentation and either backward
compatibility or an explicit deprecation path.

Numerical kernels
-----------------

RDF, MSD and VACF use compiled Cython kernels when available and NumPy/Python
fallbacks otherwise. The compiled and fallback implementations must keep the
same signature, normalization and edge-case behavior. A kernel change therefore
requires tests of both implementations and a direct parity test between them.

File parsing and logging belong outside numerical kernels. Kernels should accept
validated arrays and scalar parameters and return numerical results without
creating files.

Analysis-table contract
-----------------------

Scientific columns are defined by
``PQAnalysis/analysis/_output_schemas.py``. Each
:class:`~PQAnalysis.analysis.output.AnalysisSchema` records stable ASCII field
identifiers, display symbols, units and an optional xmgrace projection.

Writers convert numerical results into an
:class:`~PQAnalysis.analysis.output.AnalysisTable`. Native formatting may retain
legacy row layout, while CSV, TSV and XVG use the same schema. Construct all
requested writers before the calculation starts so an existing output path
fails before expensive work or partial output occurs.
