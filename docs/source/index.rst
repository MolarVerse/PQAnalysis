PQAnalysis
==========

PQAnalysis provides command-line and Python tools for quantitative analysis of
PQ molecular-dynamics simulations. It reads structures, trajectories,
velocities and Hessians, then produces documented scientific tables for
structural, transport and vibrational observables.

:doc:`Get started <getting-started>` | :doc:`Choose an analysis <analyses/index>` |
:doc:`Python functions <reference/functions>` | :doc:`Develop PQAnalysis <developerGuide/developerGuide>`

Quick start
-----------

PQAnalysis requires Python 3.12 or newer.

.. code-block:: console

   $ python -m pip install pqanalysis
   $ pqanalysis rdf rdf.in

The output filename in an analysis input file selects native text, CSV, TSV or
XVG. Repeat ``--export`` to write several formats in the same run.

.. code-block:: console

   $ pqanalysis rdf rdf.in --export rdf.csv --export rdf.xvg

Analysis methods
----------------

.. list-table:: Implemented observables
   :class: pq-record-table pq-method-table
   :header-rows: 1
   :widths: 24 38 38

   * - Method
     - Required data
     - Reported quantity
   * - :doc:`Radial distribution <analyses/rdf>`
     - Positions and periodic cell
     - :math:`g_{AB}(r)` and cumulative coordination
   * - :doc:`Mean square displacement <analyses/msd>`
     - Positions and periodic cell
     - Cartesian MSD and diffusion fits
   * - :doc:`VACF and spectra <analyses/vacf>`
     - Velocities, sampling interval and optional charges
     - Normalized correlation and wavenumber spectrum
   * - :doc:`Vibrational analysis <analyses/vibrations>`
     - Structure, masses and Cartesian Hessian
     - Normal modes, wavenumbers and optional IR intensities
   * - :doc:`Momentum diagnostic <analyses/momentum>`
     - Velocities and atomic masses
     - Frame-resolved total linear momentum

Python interface
----------------

The public analysis functions use the same validated input readers and
scientific kernels as the command line:

.. code-block:: python

   from PQAnalysis.analysis import rdf, read_analysis_table

   rdf("rdf.in", export_files=["rdf.csv"])
   table = read_analysis_table("rdf.csv")

The :doc:`function index <reference/functions>` groups analysis workflows,
numerical methods, scientific-table operations and simulation-file I/O by
task.

Development
-----------

New methods follow a documented path from estimator and validation evidence to
the public API, CLI and schema-backed output. See
:doc:`Adding an Analysis <developerGuide/adding-analysis>` for the required
implementation steps and :doc:`Architecture <developerGuide/architecture>` for
package ownership boundaries.

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Use PQAnalysis

   getting-started
   analyses/index
   Python Functions <reference/functions>
   Command Line <reference/cli>
   Files and Formats <data/index>
   Package Reference <reference/api>

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Develop PQAnalysis

   developerGuide/developerGuide
