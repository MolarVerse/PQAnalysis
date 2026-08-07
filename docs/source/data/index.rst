Files and Formats
=================

PQAnalysis separates simulation data, analysis configuration and output-table
serialization. File extensions select output formats, while input content and
explicit engine options determine how trajectories are read.

This section covers four interfaces:

* :ref:`inputFile` defines the key-value grammar.
* :ref:`analysisOutputFiles` defines table fields, symbols and units.
* :doc:`../reference/cli` documents table, structure and trajectory conversion.
* :doc:`../reference/api` identifies readers, writers and trajectory objects.

Analysis configuration
----------------------

RDF, MSD, VACF and vibrational calculations use key-value input files. Lists
may be written in brackets or as multiline values according to the
:ref:`inputFile` grammar. Relative filenames are resolved by the process
running the command, so reproducible workflows should execute from a known run
directory.

Trajectories and engines
------------------------

Analysis commands default to PQ conventions. Use ``--engine`` when reading a
supported alternative convention. Position analyses require coordinates and a
consistent atom count; MSD additionally needs periodic cells for unwrapping.
VACF and momentum analyses require velocity data. Molecular exclusions in RDF
require topology information from a restart and moldescriptor.

The Python format definitions are documented by
:class:`PQAnalysis.traj.formats.MDEngineFormat` and
:class:`PQAnalysis.traj.formats.TrajectoryFormat`.

Selections
----------

Analysis selections are parsed by :class:`PQAnalysis.topology.selection.Selection`.
Use elemental or atom-name selections for simple systems and full atom
information when residue-aware selection is required. Always verify that the
selection contains the intended atoms; normalization and statistical quality
depend directly on its population.

Output and conversion
---------------------

Native, CSV, TSV and PQAnalysis-generated XVG tables are mutually convertible.
The converter detects input content rather than trusting the extension and can
write several outputs atomically:

.. code-block:: console

   $ pqanalysis convert rdf.xvg \
       -o rdf.dat \
       -o rdf.csv \
       -o rdf.tsv

No output is written if any requested destination already exists. Use
``--mode o`` only when intentional replacement is acceptable.

.. toctree::
   :hidden:
   :maxdepth: 1

   Analysis input files <../userGuide/inputFile>
   Analysis output files <../userGuide/analysisOutputFiles>
