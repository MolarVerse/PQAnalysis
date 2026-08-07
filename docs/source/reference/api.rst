Python API
==========

The public analysis wrappers accept the same input files as the command line
and are the simplest integration points:

.. list-table:: Analysis entry points
   :header-rows: 1
   :widths: 34 66

   * - Function
     - Purpose
   * - :func:`PQAnalysis.analysis.rdf.api.rdf`
     - Radial distribution analysis
   * - :func:`PQAnalysis.analysis.msd.api.msd`
     - Mean square displacement analysis
   * - :func:`PQAnalysis.analysis.vacf.api.vacf`
     - Velocity or charge-flux correlation analysis
   * - :func:`PQAnalysis.analysis.vibrational.api.vibrations`
     - Vibrational analysis from a structure and Hessian
   * - :func:`PQAnalysis.analysis.momentum.api.check_momentum`
     - Frame-resolved total linear momentum

Package areas
-------------

.. grid:: 1 2 3 3
   :gutter: 2

   .. grid-item-card:: Analysis API
      :link: ../code/PQAnalysis.analysis
      :link-type: doc

      Calculations, input readers, result models and output writers.

   .. grid-item-card:: Input and output
      :link: ../code/PQAnalysis.io
      :link-type: doc

      Trajectory, restart, topology and simulation-file readers and writers.

   .. grid-item-card:: Trajectories
      :link: ../code/PQAnalysis.traj
      :link-type: doc

      Engine formats, trajectory containers and high-level operations.

   .. grid-item-card:: Atomic systems
      :link: ../code/PQAnalysis.atomic_system
      :link-type: doc

      Atomic coordinates, cells and topology-bearing systems.

   .. grid-item-card:: Topology and selection
      :link: ../code/PQAnalysis.topology
      :link-type: doc

      Selections, residues, bonded topology and SHAKE definitions.

   .. grid-item-card:: Complete package index
      :link: ../code/PQAnalysis
      :link-type: doc

      Every generated module, class, function and exception.

Use the curated analysis guides for physical conventions and the generated
reference for signatures and implementation-level details.
