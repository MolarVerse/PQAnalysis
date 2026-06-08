.. _userGuide:

##########
User Guide
##########

Command Line Interface
======================

The PQAnalysis package does not only provide an API but also a number of different command line tools. These tools can be categorized into two groups primary groups: pure command line tools and tools that are based on an input file. 

Input file based tools
----------------------

For more details on the grammar and syntax of the input file see :ref:`inputFile`.

- :ref:`rdf<cli.rdf>`

RDF input files
^^^^^^^^^^^^^^^

Basic RDF calculations only need a trajectory, selections, bin settings,
and an output file. A restart file is not required for this case:

.. code-block:: text

    reference_selection = H
    target_selection = O
    delta_r = 0.05
    out_file = rdf.out
    traj_files = trajectory.xyz

Restart files and moldescriptor files are only needed when the calculation
requires molecular topology information. For example,
:code:`no_intra_molecular = True` excludes pairs from the same molecule.
If the files are not given explicitly, PQAnalysis tries to infer
:code:`trajectory.rst` from :code:`trajectory.xyz` and
:code:`moldescriptor.dat` from the trajectory directory:

.. code-block:: text

    reference_selection = H
    target_selection = O
    delta_r = 0.05
    out_file = rdf_inter.out
    traj_files = trajectory.xyz
    no_intra_molecular = True

Explicit :code:`restart_file` and :code:`moldescriptor_file` values are used
as-is. If both files are given and :code:`no_intra_molecular` is omitted,
:code:`no_intra_molecular` defaults to :code:`True`. If
:code:`no_intra_molecular` is set to :code:`False`, intra molecular pairs are
included. Inferred and defaulted values are written to the normal PQAnalysis
log output.

Pure command line tools
-----------------------

- :ref:`continue_input<cli.continue_input>`
- :ref:`rst2xyz<cli.rst2xyz>`
- :ref:`traj2qmcfc<cli.traj2qmcfc>`
- :ref:`traj2box<cli.traj2box>`
