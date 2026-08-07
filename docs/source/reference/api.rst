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

.. list-table:: Generated package reference
   :header-rows: 1
   :widths: 34 66

   * - Package
     - Scope
   * - :doc:`Analysis <../code/PQAnalysis.analysis>`
     - Calculations, input readers, result models and output writers
   * - :doc:`Input and output <../code/PQAnalysis.io>`
     - Trajectory, restart, topology and simulation-file readers and writers
   * - :doc:`Trajectories <../code/PQAnalysis.traj>`
     - Engine formats, trajectory containers and high-level operations
   * - :doc:`Atomic systems <../code/PQAnalysis.atomic_system>`
     - Coordinates, cells and topology-bearing systems
   * - :doc:`Topology and selection <../code/PQAnalysis.topology>`
     - Selections, residues, bonded topology and SHAKE definitions
   * - :doc:`Package index <../code/PQAnalysis>`
     - Generated module hierarchy

Use the analysis guides for physical conventions and the generated reference
for signatures and implementation details.
