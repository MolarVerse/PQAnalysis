Files and Formats
=================

PQAnalysis separates simulation data, analysis configuration and output-table
serialization. File extensions select output formats, while input content and
explicit engine options determine how trajectories are read.

Four file contracts govern analysis workflows:

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

Analysis keys such as ``reference_selection`` use the string language
documented in :doc:`selections`.

Output and conversion
---------------------

Native, CSV, TSV and PQAnalysis-generated XVG tables are mutually convertible
with ``pqanalysis convert``; see :ref:`analysisOutputFiles`.

.. toctree::
   :hidden:
   :maxdepth: 1

   Analysis input files <../userGuide/inputFile>
   Atom selections <selections>
   Analysis output files <../userGuide/analysisOutputFiles>
