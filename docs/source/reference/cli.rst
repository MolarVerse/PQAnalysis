Command-Line Reference
======================

``pqanalysis`` dispatches all supported commands from one executable. Every
subcommand also provides local help:

.. code-block:: console

   $ pqanalysis --help
   $ pqanalysis rdf --help

Analysis commands
-----------------

.. list-table:: Analysis commands
   :class: pq-command-table
   :header-rows: 1
   :widths: 24 50 26

   * - Command
     - Purpose
     - Primary input
   * - :ref:`rdf <cli.rdf>`
     - Radial distribution and cumulative coordination
     - Input file
   * - :ref:`msd <cli.msd>`
     - Mean square displacement and diffusion fits
     - Input file
   * - :ref:`vacf <cli.vacf>`
     - Velocity or charge-flux correlation and spectra
     - Input file
   * - :ref:`vibrations <cli.vibrations>`
     - Hessian normal modes and optional IR intensities
     - Input file
   * - :ref:`check_momentum <cli.check_momentum>`
     - Total linear momentum per velocity frame
     - Trajectory files
   * - :ref:`build_spectrum <cli.build_spectrum>`
     - Gaussian or Lorentzian broadening of discrete lines
     - Line table

Analysis commands accept ``--export FILE`` where applicable. Repeat the option
to produce several output formats without repeating the calculation.

Table conversion
----------------

``pqanalysis convert`` reads native, CSV, TSV or PQAnalysis-generated XVG
analysis tables and writes one or more target formats. It preserves complete
schemas and hidden XVG data sets. See
:doc:`the generated option reference <../code/PQAnalysis.cli.convert>`.

Structure and trajectory conversion
-----------------------------------

.. list-table:: Structure and trajectory commands
   :class: pq-command-table
   :header-rows: 1
   :widths: 28 72

   * - Command
     - Purpose
   * - :ref:`rst2xyz <cli.rst2xyz>`
     - Convert a PQ restart structure to XYZ
   * - :ref:`xyz2rst <cli.xyz2rst>`
     - Convert XYZ coordinates to a PQ restart structure
   * - :ref:`xyz2gen <cli.xyz2gen>`
     - Convert XYZ to DFTB+ GEN
   * - :ref:`gen2xyz <cli.gen2xyz>`
     - Convert DFTB+ GEN to XYZ
   * - :ref:`traj2box <cli.traj2box>`
     - Extract periodic box data from trajectories
   * - :ref:`traj2extxyz <cli.traj2extxyz>`
     - Write extended XYZ trajectories with selected metadata
   * - :ref:`traj2qmcfc <cli.traj2qmcfc>`
     - Convert trajectories to QMCFC conventions

Simulation-support commands
---------------------------

.. list-table:: Simulation-support commands
   :class: pq-command-table
   :header-rows: 1
   :widths: 28 72

   * - Command
     - Purpose
   * - :ref:`continue_input <cli.continue_input>`
     - Continue indexed PQ or QMCFC input/output sequences
   * - :ref:`add_molecules <cli.add_molecules>`
     - Add molecular structures to an existing system
   * - :ref:`build_nep_traj <cli.build_nep_traj>`
     - Assemble a trajectory from nudged-elastic-band data

Commands refuse unsafe output replacement by default. Consult each generated
reference page for its supported writing modes and format-specific options.
