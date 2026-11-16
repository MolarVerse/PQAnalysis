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
:ref:`inputFile` grammar. Relative filenames inside an input file resolve
against the directory that contains the file.

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

PQ energy time axes
-------------------

PQ energy files store integration-step counts in their first column. The
``SIMULATION-TIME`` entry in a PQ ``.info`` file describes the final summary,
not the unit of those step counts. ``EnergyFileReader`` therefore reports
``simulation_time_unit == "step"`` when no integration timestep is known.
The existing ``simulation_time`` attribute keeps its name for compatibility.
Do not infer a timestep from the final summary or the output frequency.

With PQ info labels and a native ``# timestep = 0.5 fs`` energy header,
step counts are converted to ps. PQ can write this header using
``include_output_metadata = true``. For a legacy file, supply an explicit
constant timestep from its simulation input::

    energy = EnergyFileReader("run.en", timestep_fs=0.5).read()

Step 100 then corresponds to 0.05 ps, irrespective of output stride.
Restarted files retain their recorded step offsets. This conversion assumes
that the recorded counts use the supplied constant timestep; it does not
reconstruct histories with changing timestep sizes or unknown time origins.

Conflicting, non-positive or non-finite timestep metadata are rejected.
A multi-file read cannot mix known and unknown timestep metadata unless an
explicit consistent timestep is supplied. QMCFC behavior is unchanged.
Reading without info labels leaves the raw numeric columns unchanged.
