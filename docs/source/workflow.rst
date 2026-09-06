From trajectory to figure
=========================

This page runs one continuous Python session on ``examples/water``: load the
trajectory, inspect it, compute an RDF and an MSD, write a table, and plot.
The figure is the tutorial fixture, not a material property.

Load and inspect
----------------

.. code-block:: python

   from PQAnalysis.io import read_trajectory

   traj = read_trajectory("examples/water/trajectory.xyz")
   frame = traj[0]
   print(len(traj))                 # 25
   print(frame.n_atoms)             # 3
   print([atom.name for atom in frame.atoms])  # ['O', 'H', 'H']
   print(frame.cell.box_lengths)    # [10. 10. 10.]

Check: 25 frames, three atoms, cubic 10 Å cell. If the atom count is wrong,
stop; every later normalization will be wrong.

RDF of oxygen around hydrogen
-----------------------------

.. code-block:: python

   from PQAnalysis.analysis import RDF
   from PQAnalysis.io import TrajectoryReader

   r, g, coordination, shell, residual = RDF(
       TrajectoryReader("examples/water/trajectory.xyz"),
       reference_species="O",
       target_species="H",
       delta_r=0.5,
       r_max=4.0,
   ).run()
   print(r[0], g[0])                # 0.25  0.0
   print(coordination[-1])          # 2.0  (both hydrogens)

Check: ``r_max = 4.0`` is below half the shortest box edge (5 Å), so the
minimum-image clamp is not in play. The first bin is empty (no O–H contact
inside 0.5 Å); coordination saturates at 2.

MSD of the oxygen
-----------------

.. code-block:: python

   from PQAnalysis.analysis import MSD

   analysis = MSD(
       TrajectoryReader("examples/water/trajectory.xyz"),
       target_species="O",
       window=8,
       gap=2,
       time_step=0.001,
       fit_window=4,
   )
   lags, msd_x, msd_y, msd_z, msd_tot = analysis.run()
   print(analysis.fit_results)

Check: one oxygen is selected. The Einstein fit uses the last four lag points
of an eight-frame window. That is enough to exercise the API; it is not a
diffusion measurement. See :doc:`analyses/msd`.

Write a table and plot
----------------------

.. code-block:: python

   from PQAnalysis.analysis import rdf, read_analysis_table

   rdf("examples/water/rdf.in", export_files=["rdf.csv"])
   table = read_analysis_table("rdf.csv")
   r = table.column("r_i")
   g = table.column("g_r_i")

The figure below is produced from the same fixture during the documentation
build:

.. plot:: _plots/workflow.py
   :alt: RDF and MSD of the bundled isolated-water tutorial fixture
   :caption: Isolated water molecule, 25 frames, 10 Å box. The RDF peak is the
      intramolecular O–H distance. The MSD is the oxygen of that one molecule.
      Neither panel is a bulk-liquid result.

``read_trajectory`` loads every frame. For long files use
:class:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader` instead,
as the RDF and MSD constructors do above. Next: :doc:`data/selections` and the
method pages under :doc:`analyses/index`.
