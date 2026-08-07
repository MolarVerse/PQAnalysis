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

Create ``rdf.in`` beside a PQ trajectory named ``trajectory.xyz``:

.. code-block:: text

   traj_files = trajectory.xyz
   reference_selection = O
   target_selection = H
   delta_r = 0.05
   out_file = rdf.dat

Run the calculation:

.. code-block:: console

   $ pqanalysis rdf rdf.in

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

PQAnalysis refuses to overwrite an existing output unless replacement is
requested explicitly with ``--mode o``.

Next steps
----------

* :doc:`analyses/index` compares the physical observables and required data.
* :doc:`data/index` defines input grammar, trajectory conventions and table
  formats.
* :doc:`reference/cli` lists commands and options.
* :doc:`reference/api` identifies the Python analysis and I/O entry points.
