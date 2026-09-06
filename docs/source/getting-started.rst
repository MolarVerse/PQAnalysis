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

The ``rdf.dat`` from the previous step converts without rerunning the
analysis:

.. code-block:: console

   $ pqanalysis convert rdf.dat -o rdf.csv -o rdf.xvg

To write several formats in a single analysis run, repeat ``--export``:

.. code-block:: console

   $ rm rdf.dat rdf.csv rdf.xvg
   $ pqanalysis rdf rdf.in \
       --export rdf.csv \
       --export rdf.tsv \
       --export rdf.xvg

The ``rm`` is needed because PQAnalysis refuses to overwrite an existing
output file. ``convert`` and the other support tools accept ``--mode o`` to
request replacement explicitly; the input-file driven analyses have no
overwrite flag, so move or delete the old output first.

Use the Python API
------------------

Still in ``examples/water``,
:func:`~PQAnalysis.analysis.output.read_analysis_table` reloads the table
written above with its column metadata:

.. code-block:: python

   from PQAnalysis.analysis import read_analysis_table

   table = read_analysis_table("rdf.csv")
   print(table.column("r_i")[:5])
   print(table.column("g_r_i")[:5])

The analysis itself runs from Python as ``rdf("rdf.in", export_files=["rdf.csv"])``
once the earlier outputs are removed. The input file can also be given by
path from another directory; filenames inside it resolve relative to the
input file. :doc:`python-api` covers the wrappers and the in-memory analysis
objects such as :class:`~PQAnalysis.analysis.rdf.rdf.RDF`.

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
