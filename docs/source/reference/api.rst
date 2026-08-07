Package Reference
=================

Start with the :doc:`functions` page for callable analysis, numerical and I/O
interfaces. This page maps the principal data types and generated module
reference.

Core types
----------

.. list-table:: Principal data types
   :class: pq-record-table pq-types-table
   :header-rows: 1
   :widths: 34 66

   * - Type
     - Role
   * - :class:`PQAnalysis.atomic_system.AtomicSystem`
     - One structure with coordinates, cell and topology
   * - :class:`PQAnalysis.traj.Trajectory`
     - Ordered atomic-system frames
   * - :class:`PQAnalysis.topology.Topology`
     - Atoms, residues, molecular identity and bonded topology
   * - :class:`PQAnalysis.topology.Selection`
     - Atom selection parser and index resolution
   * - :class:`PQAnalysis.analysis.output.AnalysisTable`
     - Numerical analysis data coupled to scientific column metadata
   * - :class:`PQAnalysis.analysis.output.AnalysisSchema`
     - Stable fields, symbols, units and plot defaults

Package areas
-------------

.. list-table:: Generated package reference
   :class: pq-record-table pq-package-reference-table
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

Use the analysis guides for physical conventions, the function index for
callable workflows and the generated reference for implementation details.
