Mean Square Displacement
========================

The mean square displacement measures translational motion over a lag time.
For Cartesian component :math:`\alpha`, PQAnalysis evaluates multiple time
origins according to

.. math::

   \mathrm{MSD}_{\alpha}(\tau) =
   \left\langle [r_{i\alpha}(t+\tau)-r_{i\alpha}(t)]^2 \right\rangle_{i,t}.

Coordinates are unwrapped with the periodic cell before displacements are
accumulated.

Minimal input
-------------

.. code-block:: text

   traj_files = trajectory.xyz
   target_selection = O
   out_file = msd.dat
   window = 1000
   gap = 10
   time_step = 0.001
   fit_window = 200

.. code-block:: console

   $ pqanalysis msd msd.in

``window`` is the largest lag in frames and must be divisible by ``gap``.
``gap`` controls the spacing between time origins. ``time_step`` is expressed
in ps and enables diffusion fitting; ``fit_window`` selects the trailing
points used by that fit.

Interpretation
--------------

The output contains the lag index and the x, y and z components in
Angstrom squared. Their sum is the total three-dimensional MSD. In an
isotropic diffusive regime,

.. math::

   D = \frac{1}{6}\frac{d}{dt}\mathrm{MSD}_{\mathrm{total}}(t).

PQAnalysis also fits each Cartesian component with the corresponding
one-dimensional factor. The resulting coefficients, uncertainties and
:math:`R^2` values are written to the log file in m\ :sup:`2`/s. A fit is
physically meaningful only over a linear diffusive interval; short-time
ballistic motion and poorly sampled long lags should not be included blindly.

Output and API
--------------

See :ref:`analysis-output-msd` for the exact table layout. The input-file entry
point is :func:`PQAnalysis.analysis.msd.api.msd`; direct workflows can use
:class:`PQAnalysis.analysis.msd.msd.MSD` and inspect its total MSD and fit
results.
