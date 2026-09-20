Mean Square Displacement
========================

The mean square displacement measures translational motion over a lag time
[Allen2017]_. For Cartesian component :math:`\alpha`, PQAnalysis averages over
multiple time origins [Rahman1964]_ according to

.. math::

   \mathrm{MSD}_{\alpha}(\tau) =
   \left\langle [r_{i\alpha}(t+\tau)-r_{i\alpha}(t)]^2 \right\rangle_{i,t}.

Coordinates are unwrapped with the periodic cell before displacements are
accumulated. In an isotropic diffusive regime the Einstein relation
[Einstein1905]_ links the long-time slope to the self-diffusion coefficient,

.. math::

   D = \frac{1}{6}\frac{d}{dt}\mathrm{MSD}_{\mathrm{total}}(t),

and each Cartesian component is fitted with the one-dimensional factor
:math:`\tfrac{1}{2}`. Below, the bundled oxygen-atom fixture is plotted with
an assumed 0.5 ps frame interval; the dashed line fits the last 20 total-MSD
points.

.. plot:: _plots/msd.py
   :alt: Cartesian and total mean square displacement with a linear fit

Input
-----

.. code-block:: text

   traj_files = trajectory.xyz
   target_selection = O
   out_file = msd.dat
   window = 1000
   gap = 10
   time_step = 0.001
   fit_window = 200

Saved as ``msd.in``, it runs with:

.. code-block:: console

   $ pqanalysis msd msd.in

``window`` is the largest lag in frames and must be divisible by ``gap``, the
spacing between time origins. ``time_step`` in ps enables the diffusion fit;
``fit_window`` is the number of trailing points it uses (default
``max(2, window // 5)``). The bundled :doc:`../examples` fixture uses ``window = 8``
because it has only 25 frames. The full key table is on
:class:`~PQAnalysis.analysis.msd.msd_input_file_reader.MSDInputFileReader`.

Output
------

``out_file`` holds the lag index and the x, y and z components in Å²; their
sum is the total MSD (:ref:`analysis-output-msd`). Diffusion coefficients,
their standard errors and :math:`R^2` go to the log file in m²·s⁻¹.

Python
------

.. code-block:: python

   from PQAnalysis.analysis import MSD
   from PQAnalysis.io import TrajectoryReader

   analysis = MSD(
       TrajectoryReader("examples/water/trajectory.xyz"),
       target_species="O",
       window=8,
       gap=2,
       time_step=0.001,
       fit_window=4,
   )
   lags, msd_x, msd_y, msd_z, msd_tot = analysis.run()
   print(analysis.fit_results)  # D components in m²·s⁻¹ when time_step is set

:func:`~PQAnalysis.analysis.msd.api.msd` runs an input file instead.

Before you trust it
-------------------

* Fit only where the MSD is linear in time. A high :math:`R^2` on a ballistic
  or caged segment still yields a finite, meaningless :math:`D`.
* The fit is anchored at the largest lag. Shrink ``fit_window`` to skip the
  short-time regime; shrink ``window`` to move away from the noisy tail.
* Every lag is averaged over ``(n_frames - window) // gap`` origins. Keep the
  trajectory much longer than ``window``.
* The printed uncertainty is a lower bound, and no finite-size (Yeh-Hummer)
  correction is applied.
* Atoms must not move more than half a box edge between written frames, or
  unwrapping fails silently.

Each point is derived in :doc:`msd-details`.

.. toctree::
   :hidden:

   msd-details
