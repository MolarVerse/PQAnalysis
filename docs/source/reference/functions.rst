.. _function-index:

Function Index
==============

Callable Python interfaces are grouped below by task. Analysis wrappers accept
the same input files as their command-line counterparts. Lower-level numerical
functions operate on arrays or PQAnalysis objects and do not parse command-line
arguments.

Analysis workflows
------------------

.. autosummary::

   ~PQAnalysis.analysis.rdf.api.rdf
   ~PQAnalysis.analysis.msd.api.msd
   ~PQAnalysis.analysis.vacf.api.vacf
   ~PQAnalysis.analysis.vibrational.api.vibrations
   ~PQAnalysis.analysis.momentum.api.check_momentum
   ~PQAnalysis.analysis.spectrum_broadening.api.build_spectrum

Numerical methods
-----------------

.. autosummary::

   ~PQAnalysis.analysis.vacf.spectrum.apodization_window
   ~PQAnalysis.analysis.vacf.spectrum.vacf_spectrum
   ~PQAnalysis.analysis.spectrum_broadening.spectrum_broadening.alpha_from_fwhm
   ~PQAnalysis.analysis.spectrum_broadening.spectrum_broadening.fwhm_from_alpha
   ~PQAnalysis.analysis.spectrum_broadening.spectrum_broadening.wavenumber_grid
   ~PQAnalysis.analysis.spectrum_broadening.spectrum_broadening.broaden
   ~PQAnalysis.analysis.vibrational.vibrational_analysis.calculate
   ~PQAnalysis.analysis.vibrational.vibrational_analysis.read_hessian_file
   ~PQAnalysis.analysis.vibrational.vibrational_analysis.select_mode_indices

Analysis tables
---------------

.. autosummary::

   ~PQAnalysis.analysis.output.infer_output_format
   ~PQAnalysis.analysis.output.read_analysis_table
   ~PQAnalysis.analysis.output.write_analysis_table
   ~PQAnalysis.analysis.output.convert_analysis_output

Structure and trajectory I/O
----------------------------

.. autosummary::

   ~PQAnalysis.io.traj_file.api.read_trajectory
   ~PQAnalysis.io.traj_file.api.read_trajectory_generator
   ~PQAnalysis.io.traj_file.api.write_trajectory
   ~PQAnalysis.io.traj_file.api.calculate_frames_of_trajectory_file
   ~PQAnalysis.io.restart_file.api.read_restart_file
   ~PQAnalysis.io.restart_file.api.write_restart_file
   ~PQAnalysis.io.gen_file.api.read_gen_file
   ~PQAnalysis.io.gen_file.api.write_gen_file
   ~PQAnalysis.io.topology_file.api.read_topology_file
   ~PQAnalysis.io.topology_file.api.write_topology_file
   ~PQAnalysis.io.box_reader.read_box
   ~PQAnalysis.io.optimizer_file_reader.read_optimizer_file
   ~PQAnalysis.io.write_api.write
   ~PQAnalysis.io.write_api.write_box
   ~PQAnalysis.traj.api.check_trajectory_pbc
   ~PQAnalysis.traj.api.check_trajectory_vacuum

Format conversion
-----------------

.. autosummary::

   ~PQAnalysis.io.conversion_api.rst2xyz
   ~PQAnalysis.io.conversion_api.xyz2rst
   ~PQAnalysis.io.conversion_api.xyz2gen
   ~PQAnalysis.io.conversion_api.gen2xyz
   ~PQAnalysis.io.conversion_api.traj2box
   ~PQAnalysis.io.conversion_api.traj2extxyz
   ~PQAnalysis.io.conversion_api.traj2qmcfc

See :doc:`api` for classes, enums, exceptions and the generated module
hierarchy.
