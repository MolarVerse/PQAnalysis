Learn the Python API
====================

PQAnalysis exposes the same scientific kernels from Python that the command
line uses. There are two layers:

* **File wrappers** such as :func:`~PQAnalysis.analysis.rdf.api.rdf` read the
  same ``.in`` files as ``pqanalysis rdf`` and write the same tables.
* **Analysis objects** such as :class:`~PQAnalysis.analysis.rdf.rdf.RDF` take a
  trajectory or reader and return NumPy arrays. Use them when the data are
  already in memory or when you want to post-process without writing a file.

Both layers are single-use for a given analysis object: call ``run()`` once,
then construct a new object for another calculation.

File wrappers and analysis tables
---------------------------------

The shortest Python path mirrors the command line. From the repository root,
using :doc:`examples`:

.. code-block:: python

   from PQAnalysis.analysis import rdf, read_analysis_table

   rdf("examples/water/rdf.in", export_files=["rdf.csv"])
   table = read_analysis_table("rdf.csv")

   r = table.column("r_i")
   g = table.column("g_r_i")
   print(table.schema.fields)
   print(r[:5], g[:5])

:func:`~PQAnalysis.analysis.output.read_analysis_table` accepts native text,
CSV, TSV or XVG and reconstructs the scientific column schema. Stable field
names (``r_i``, ``g_r_i``, ``lag``, ``normalized_correlation``, …) are listed
in :ref:`analysisOutputFiles`.

The same pattern works for the other file-driven analyses:

.. code-block:: python

   from PQAnalysis.analysis import msd, vacf, vibrations, check_momentum

   msd("examples/water/msd.in", export_files=["msd.csv"])
   vacf("examples/water/vacf.in", export_files=["vacf.csv"])
   vibrations("examples/water/vibrations.in", export_files=["wavenumbers.csv"])
   norms = check_momentum(
       "examples/water/trajectory.vel",
       output="momentum.dat",
       selection="all",
   )

``check_momentum`` is the exception: it takes trajectory paths directly and
returns the scaled momentum norms as a NumPy array.

Analysis objects
----------------

When frames are already loaded, construct the analysis class and call
``run()``. The example below builds a two-frame orthorhombic water-like cell
and computes a short RDF without an input file:

.. code-block:: python

   import numpy as np

   from PQAnalysis.analysis import RDF
   from PQAnalysis.atomic_system import AtomicSystem
   from PQAnalysis.core import Atom, Cell
   from PQAnalysis.traj import Trajectory

   cell = Cell(10.0, 10.0, 10.0)
   atoms = [Atom("O"), Atom("H"), Atom("H")]
   frames = [
       AtomicSystem(
           atoms=atoms,
           pos=np.array([
               [0.0, 0.0, 0.0],
               [1.0, 0.0, 0.0],
               [0.0, 1.0, 0.0],
           ]),
           cell=cell,
       ),
       AtomicSystem(
           atoms=atoms,
           pos=np.array([
               [0.1, 0.0, 0.0],
               [1.1, 0.0, 0.0],
               [0.0, 1.1, 0.0],
           ]),
           cell=cell,
       ),
   ]

   r, g, coordination, shell, residual = RDF(
       Trajectory(frames),
       reference_species="O",
       target_species="H",
       delta_r=0.1,
       r_max=4.0,
   ).run()

Readers from :func:`~PQAnalysis.io.traj_file.api.read_trajectory` and
:class:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader` are valid
``traj`` arguments as well. File-backed orthorhombic inputs can still take the
legacy-compatible fast paths documented on each method page.

Reading trajectories
--------------------

.. code-block:: python

   from PQAnalysis.io import read_trajectory, TrajectoryReader

   traj = read_trajectory("examples/water/trajectory.xyz")  # loads all frames
   reader = TrajectoryReader("examples/water/trajectory.xyz")  # streams frames

Use ``read_trajectory`` for small systems and interactive work. Prefer
``TrajectoryReader`` (or the analysis file wrappers) for long trajectories so
frames are not held in memory at once.

Atom strings such as ``"O"``, ``"0..2"`` and ``"*|H"`` are documented in
:doc:`data/selections`. A continuous load → RDF → MSD → figure session is
:doc:`workflow`.

Where each method is documented
-------------------------------

Method pages give a focused Python recipe next to the CLI input:

* :doc:`analyses/rdf` — ``rdf()`` and :class:`~PQAnalysis.analysis.rdf.rdf.RDF`
* :doc:`analyses/msd` — ``msd()`` and :class:`~PQAnalysis.analysis.msd.msd.MSD`
* :doc:`analyses/vacf` — ``vacf()``, :class:`~PQAnalysis.analysis.vacf.vacf.VACF`
  and :func:`~PQAnalysis.analysis.vacf.spectrum.vacf_spectrum`
* :doc:`analyses/vibrations` — ``vibrations()`` and
  :func:`~PQAnalysis.analysis.vibrational.vibrational_analysis.calculate`
* :doc:`analyses/momentum` — ``check_momentum()`` and
  :class:`~PQAnalysis.analysis.momentum.momentum.Momentum`

The full callable index is :doc:`reference/functions`. Generated class pages
live under :doc:`reference/api`.
