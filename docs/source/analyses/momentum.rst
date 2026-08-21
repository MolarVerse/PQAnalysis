Total Linear Momentum
=====================

For every velocity frame, PQAnalysis evaluates the total linear momentum of the
selected atoms,

.. math::

   \mathbf{P}(t) = \sum_i m_i\mathbf{v}_i(t) ,

and writes its scaled norm

.. math::

   p(t) = \sigma\,\lVert\mathbf{P}(t)\rVert
        = \sigma\left\lVert\sum_i m_i\mathbf{v}_i(t)\right\rVert ,

one value per frame. Here :math:`m_i` is the atomic mass in amu taken from the
topology, :math:`\mathbf{v}_i` is the atomic velocity read from the trajectory,
and :math:`\sigma` is the ``scale`` factor. This is a diagnostic for
center-of-mass drift and momentum conservation [Allen2017]_, not a substitute
for inspecting the thermostat, constraints or integration scheme.

Units
-----

PQ velocity trajectories store velocities in Å·s⁻¹, so :math:`\mathbf{P}` is in
amu·Å·s⁻¹. The default :math:`\sigma = 10^{-15}` converts that to
amu·Å·fs⁻¹, the unit of the second output column. Any other value of
``--scale`` simply multiplies the norm, and it is then the user's
responsibility to make :math:`\sigma` match the velocity convention of the
input trajectory.

Run the diagnostic
------------------

.. code-block:: console

   $ pqanalysis check_momentum velocity.vel \
       --selection all \
       --output momentum.dat

The output contains a one-based frame index and the scaled momentum norm. Use
``--scale`` when the input convention differs from Å·s⁻¹, and ``--selection``
to restrict the sum to a subset of atoms.

Precision and the noise floor
-----------------------------

:math:`\mathbf{P}` is a heavily cancelling sum. In a well-behaved simulation
the individual terms :math:`m_i\mathbf{v}_i` are large and nearly cancel, so
the surviving norm is smaller than any single term by many orders of magnitude.
The smallest norm that still carries information is therefore set by the
relative precision :math:`\varepsilon` of the velocity values, not by the
accumulator, which is always float64:

.. math::

   \lVert\mathbf{P}\rVert \gtrsim
   \varepsilon\sum_i m_i\lVert\mathbf{v}_i\rVert .

Two code paths set :math:`\varepsilon` differently:

* File-backed PQ and QMCFC velocity trajectories — files recognized as ``.vel``
  or ``.velocs`` and read through ``check_momentum`` — are parsed directly as
  float64. Here :math:`\varepsilon` is the precision of the text itself, that
  is, the number of significant digits the MD engine wrote.
* Other xyz-family trajectory formats read from file keep the single-precision
  arrays produced by the general frame reader, giving
  :math:`\varepsilon\approx 1.2\times10^{-7}` before the values are widened to
  float64 for the sum.
* Trajectory objects built in memory keep the precision of the velocity arrays
  they were given, so a float64 array is summed without any loss.

As a rule of thumb, a scaled norm below roughly :math:`10^{-7}` of
:math:`\sum_i m_i\lVert\mathbf{v}_i\rVert` is parsing and round-off noise
rather than physical drift. Compare against that scale before reading anything
into an absolute value.

The compatibility path multiplies and sums atoms in the same order as the
legacy ``equipartition.jl`` calculation [thhTools]_, so residuals near the
float64 noise floor are reproduced bit for bit. Native output uses 17
significant digits, so reloading it as float64 preserves each calculated value
exactly.

Validity and interpretation
---------------------------

A trustworthy result is a flat trace: :math:`p(t)` fluctuating around the noise
floor with no trend over the whole trajectory. The shape of the series carries
the information, not any single value.

* **A systematic increase** indicates center-of-mass drift — the classic
  signature of an integration time step that is too large, of accumulated
  round-off, or of a thermostat that adds momentum without removing it.
* **A step** at one frame usually marks a restart, a velocity reassignment or a
  change of ensemble rather than a physical process.
* **A flat trace at a large value** means the simulation started with non-zero
  total momentum. It is conserved, but the center of mass is translating, and
  MSD or diffusion coefficients from that trajectory are biased unless the
  drift is removed.

The diagnostic does not apply, or must be read differently, in these cases:

* **Partial selections.** Momentum conservation is a statement about the whole
  system. The momentum of a subset of atoms obeys no conservation law and
  fluctuates by construction, so ``--selection`` is useful for locating which
  species carries a drift, not for testing conservation.
* **Systems with external forces.** Walls, position restraints, frozen atoms,
  external fields and momentum-removing thermostats break translational
  invariance on purpose. A non-conserved momentum is then the expected result.
* **Unknown masses.** Every selected atom must have a known mass; the analysis
  refuses to run otherwise, because a missing mass would silently change the
  sum.
* **Equipartition and temperature.** This is a single vector sum over the
  system. It says nothing about how kinetic energy is distributed over degrees
  of freedom, and a conserved total momentum is no evidence of a correct
  temperature or of proper thermostatting.

Output and API
--------------

See :ref:`analysis-output-momentum` for the output schema. Python workflows
can call :func:`PQAnalysis.analysis.momentum.api.check_momentum` or use
:class:`PQAnalysis.analysis.momentum.momentum.Momentum` directly.

References
----------

* [Allen2017]_ describes conserved quantities in a molecular-dynamics run and
  the removal of center-of-mass motion.
* [thhTools]_ is the legacy program whose summation order the compatibility
  path reproduces.

Full entries are listed in :doc:`../references`.
