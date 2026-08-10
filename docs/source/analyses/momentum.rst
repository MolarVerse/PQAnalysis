Total Linear Momentum
=====================

For every velocity frame, PQAnalysis evaluates the selected atoms' total
linear momentum,

.. math::

   \mathbf{P}(t) = \sum_i m_i\mathbf{v}_i(t),

and writes its scaled norm. This is a diagnostic for center-of-mass drift and
momentum conservation, not a substitute for inspecting the thermostat,
constraints or integration scheme.

Run the diagnostic
------------------

.. code-block:: console

   $ pqanalysis check_momentum velocity.vel \
       --selection all \
       --output momentum.dat

The default scale of ``1e-15`` converts PQ velocity-trajectory values from
amu·Å·s⁻¹ to amu·Å·fs⁻¹. Use ``--scale`` when the input convention differs.

Interpretation
--------------

The output contains a one-based frame index and the scaled momentum norm. A
systematic increase can indicate center-of-mass drift. Oscillatory or noisy
behavior must be interpreted relative to the total mass, velocity scale and
numerical precision.

File-backed PQ and QMCFC velocity trajectories are parsed directly as float64.
The compatibility path multiplies and sums atoms in the same order as the
legacy ``equipartition.jl`` calculation. Native output uses 17 significant
digits, so reloading it as float64 preserves each calculated value exactly.
Trajectory objects use the numerical precision already stored in the object.

See :ref:`analysis-output-momentum` for the output schema. Python workflows
can call :func:`PQAnalysis.analysis.momentum.api.check_momentum` or use
:class:`PQAnalysis.analysis.momentum.momentum.Momentum` directly.
