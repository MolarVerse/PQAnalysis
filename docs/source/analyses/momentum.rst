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
amu Angstrom s\ :sup:`-1` to amu Angstrom fs\ :sup:`-1`. Use ``--scale`` when
the input convention differs.

Interpretation
--------------

The output contains a one-based frame index and the scaled momentum norm. A
systematic increase can indicate center-of-mass drift. Oscillatory or noisy
behavior must be interpreted relative to the total mass, velocity scale and
numerical precision.

PQ velocity trajectories are parsed in single precision. Norms below roughly
``1e-7 * sum_i(m_i * |v_i|) * scale`` are therefore parsing noise rather than
resolved physical drift.

See :ref:`analysis-output-momentum` for the output schema. Python workflows
can call :func:`PQAnalysis.analysis.momentum.api.check_momentum` or use
:class:`PQAnalysis.analysis.momentum.momentum.Momentum` directly.
