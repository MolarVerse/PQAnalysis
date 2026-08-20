Command-Line Reference
======================

``pqanalysis`` dispatches all supported commands from one executable. Every
subcommand also provides local help:

.. code-block:: console

   $ pqanalysis --help
   $ pqanalysis rdf --help

The command tables on this page are validated against the command registry at
build time: a missing, renamed or undocumented command fails the
documentation build.

Analysis commands
-----------------

.. pq-cli-table::
   :title: Analysis commands

   rdf -- Radial distribution and cumulative coordination -- Input file
   msd -- Mean square displacement and diffusion fits -- Input file
   vacf -- Velocity or charge-flux correlation and spectra -- Input file
   vibrations -- Hessian normal modes and optional IR intensities -- Input file
   check_momentum -- Total linear momentum per velocity frame -- Trajectory files
   build_spectrum -- Gaussian or Lorentzian broadening of discrete lines -- Line table

Analysis commands accept ``--export FILE`` where applicable. Repeat the option
to produce several output formats without repeating the calculation.

Table conversion
----------------

``pqanalysis convert`` reads native, CSV, TSV or PQAnalysis-generated XVG
analysis tables and writes one or more target formats. It preserves complete
schemas and hidden XVG data sets. See
:doc:`the generated option reference <../code/PQAnalysis.cli.convert>`.

.. pq-cli-covered::

   convert

Structure and trajectory conversion
-----------------------------------

.. pq-cli-table::
   :title: Structure and trajectory commands

   rst2xyz -- Convert a PQ restart structure to XYZ
   xyz2rst -- Convert XYZ coordinates to a PQ restart structure
   xyz2gen -- Convert XYZ to DFTB+ GEN
   gen2xyz -- Convert DFTB+ GEN to XYZ
   traj2box -- Extract periodic box data from trajectories
   traj2extxyz -- Write extended XYZ trajectories with selected metadata
   traj2qmcfc -- Convert trajectories to QMCFC conventions

Simulation-support commands
---------------------------

.. pq-cli-table::
   :title: Simulation-support commands

   continue_input -- Continue indexed PQ or QMCFC input/output sequences
   add_molecules -- Add molecular structures to an existing system
   build_nep_traj -- Build Neuroevolution Potential (NEP) training and test trajectories

Commands refuse unsafe output replacement by default. Consult each generated
reference page for its supported writing modes and format-specific options.
