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

PQ H5MD trajectories (``traj_format = h5md``, files ending in ``.h5md`` or
``.h5``) are read like the text trajectories, including positions,
velocities, forces, charges, cells and the potential energy. They need the
optional ``h5py`` package: ``pip install PQAnalysis[h5md]``. Values are
converted to the units of the PQ text files, for example velocities from
Angstrom fs\ :sup:`-1` to Angstrom s\ :sup:`-1`; an unknown unit is an error.

ASE trajectories (``traj_format = ase``, files ending in ``.traj``), as written
by PQ or ASE, are read the same way. They need the optional ``ase`` package:
``pip install PQAnalysis[ase]``. ASE values in eV and its own time unit are
converted to the PQ units, for example forces from eV Angstrom\ :sup:`-1` to
kcal mol\ :sup:`-1` Angstrom\ :sup:`-1`. Atoms are named by their element
symbol.

PQ Zarr trajectories (``traj_format = zarr``, directories ending in
``.zarr``) use the groups, arrays and units of the PQ H5MD files and are read
the same way. They need the optional ``zarr`` package:
``pip install PQAnalysis[zarr]``.

The Python format definitions are documented by
:class:`PQAnalysis.traj.formats.MDEngineFormat` and
:class:`PQAnalysis.traj.formats.TrajectoryFormat`.

Selections
----------

Analysis keys such as ``reference_selection`` use the string language
documented in :doc:`selections`.

Force-field import
------------------

``parm2pq`` writes the PQ setup files (moldescriptor, topology, parameter,
intra non-bonded and start file) for a system parametrized with AMBER,
CHARMM or GROMACS; see :doc:`force-field-import`. The tutorial
:doc:`../tutorials/solvated-protein` goes from a PDB file to a running PQ
simulation.

Output and conversion
---------------------

Native, CSV, TSV and PQAnalysis-generated XVG tables are mutually convertible
with ``pqanalysis convert``; see :ref:`analysisOutputFiles`.

.. toctree::
   :hidden:
   :maxdepth: 1

   Analysis input files <../userGuide/inputFile>
   Atom selections <selections>
   Force-field import <force-field-import>
   Crystal preparation <crystal-preparation>
   Analysis output files <../userGuide/analysisOutputFiles>

PQ energy file header
---------------------

A PQ energy file can name its columns in a header at the top of the file::

    # PQ energy file v1.1.0
    # timestep = 0.5 fs
    #STEP[-] TEMPERATURE[K] PRESSURE[bar] ...

The last comment line above the data has one ``NAME[unit]`` token per data
column; the timestep line is optional. Older PQ versions name the first
column ``SIMULATION-TIME[step]``; ``STEP[-]`` is read the same way, as
``SIMULATION-TIME`` in ``step``. If the first line is
``# PQ energy file``, ``EnergyFileReader`` and ``pqanalysis convert`` take
the column names and units from the header and do not need the ``.info``
file. They match the names of the ``.info`` file. A ``SIMULATION-TIME``
column in ``step`` is converted to ps as described below; other units are
used as written. If an ``.info`` file exists and names different columns or
units, the header is used and a warning is logged. Files without this header
keep using the ``.info`` file.

.. _pq-energy-time-axes:

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
