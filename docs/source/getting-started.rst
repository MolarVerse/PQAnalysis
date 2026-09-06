Getting Started
===============

Install PQAnalysis
------------------

PQAnalysis supports Python 3.12 and newer. Install the current release from
PyPI:

.. code-block:: console

   $ python -m pip install pqanalysis

Confirm that the command dispatcher and analysis commands are available:

.. code-block:: console

   $ pqanalysis --help
   $ pqanalysis rdf --help

Run a first analysis
--------------------

From a clone of the repository, the tutorial fixture in :doc:`examples` is
ready to run:

.. code-block:: console

   $ cd examples/water
   $ pqanalysis rdf rdf.in

``rdf.in`` is:

.. code-block:: text

   traj_files = trajectory.xyz
   reference_selection = O
   target_selection = H
   delta_r = 0.5
   r_max = 4.0
   out_file = rdf.dat

``rdf.dat`` contains the bin-center distance, radial distribution function,
cumulative coordination number, density-normalized shell population and
ideal-gas pair-count residual. Its commented metadata header records the field
names, scientific symbols and units. See :ref:`analysis-output-rdf` for the
exact definitions.

Choose output formats
---------------------

The output filename selects the table format. ``.csv`` and ``.tsv`` open
directly in spreadsheet software, ``.xvg`` opens in xmgrace, and any other
extension uses native PQAnalysis text.

Additional outputs do not require another analysis run:

.. code-block:: console

   $ pqanalysis rdf rdf.in \
       --export rdf.csv \
       --export rdf.tsv \
       --export rdf.xvg

Existing analysis tables can be converted later:

.. code-block:: console

   $ pqanalysis convert rdf.dat -o rdf.csv -o rdf.xvg

PQAnalysis refuses to overwrite an existing output file. Conversion and
support tools such as ``convert`` accept ``--mode o`` to request
replacement explicitly; the input-file driven analyses have no overwrite
flag, so move or delete the old output first.

Use the Python API
------------------

The same input file works from Python. Paths below are from the repository
root. The wrapper writes the table;
:func:`~PQAnalysis.analysis.output.read_analysis_table` reloads it with
column metadata:

.. code-block:: python

   from PQAnalysis.analysis import rdf, read_analysis_table

   rdf("examples/water/rdf.in", export_files=["rdf.csv"])
   table = read_analysis_table("rdf.csv")
   print(table.column("r_i")[:5])
   print(table.column("g_r_i")[:5])

For in-memory trajectories, analysis objects such as
:class:`~PQAnalysis.analysis.rdf.rdf.RDF` and recipes for every method, see
:doc:`python-api`.

Next steps
----------

* :doc:`examples` lists the bundled water fixture.
* :doc:`python-api` covers file wrappers, analysis objects and tables.
* :doc:`workflow` runs RDF and MSD in one Python session.
* :doc:`analyses/index` compares the physical observables and required data.
* :doc:`reference/functions` lists public Python workflows and numerical
  functions.
* :doc:`reference/cli` lists commands and options.
* :doc:`reference/api` identifies the Python analysis and I/O entry points.
* :doc:`developerGuide/developerGuide` documents architecture, extension and
  validation.
