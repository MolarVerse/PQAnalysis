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
:code:`no_intra_molecular = True` excludes pairs from the same molecule and
therefore requires both files:

.. code-block:: text

    reference_selection = H
    target_selection = O
    delta_r = 0.05
    out_file = rdf_inter.out
    traj_files = trajectory.xyz
    restart_file = trajectory.rst
    moldescriptor_file = moldescriptor.dat
    no_intra_molecular = True

If :code:`no_intra_molecular` is omitted or set to :code:`False`, intra
molecular pairs are included. File names are not inferred automatically;
provide :code:`restart_file` and :code:`moldescriptor_file` explicitly when
they are required.

Pure command line tools
-----------------------

- :ref:`continue_input<cli.continue_input>`
- :ref:`rst2xyz<cli.rst2xyz>`
- :ref:`traj2qmcfc<cli.traj2qmcfc>`
- :ref:`traj2box<cli.traj2box>`
