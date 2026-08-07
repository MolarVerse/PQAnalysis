Vibrational Analysis
====================

Vibrational analysis diagonalizes the mass-weighted Cartesian Hessian. Its
eigenvectors define normal modes and its eigenvalues determine signed
wavenumbers. Negative wavenumbers represent imaginary modes associated with
negative curvature of the potential-energy surface.

Minimal input
-------------

.. code-block:: text

   structure_file = structure.rst
   hessian_file = hessian.dat
   out_file = wavenumbers.dat
   normal_modes_file = normal_modes.dat
   modes_file = modes.xyz
   modes = positive
   unit = kcal
   hessian_sign = auto

.. code-block:: console

   $ pqanalysis vibrations vibrations.in

``structure_file`` may be a PQ restart or a single-frame XYZ file. ``unit``
describes the Hessian energy unit and accepts ``kcal``, ``hartree`` or ``ev``.
``hessian_sign = auto`` evaluates both supported sign conventions and chooses
the one with more non-negative vibrational modes.

Scientific checks
-----------------

* A stable, fully optimized minimum should not contain genuine imaginary
  internal modes. Small values can arise from incomplete optimization or
  numerical noise.
* Translational and rotational near-zero modes depend on boundary conditions,
  molecular geometry and numerical precision.
* IR intensities require a ``moldescriptor_file`` containing partial charges.
* The Hessian coordinate order, structure atom order and selected unit must
  agree exactly.

Mode output
-----------

``normal_modes_file`` stores the dimensionless Cartesian mode matrix.
``modes_prefix`` writes sinusoidal multi-frame XYZ animations, while
``modes_file`` writes one extended-XYZ image per selected mode with vectors and
metadata. Explicit mode numbers are one-based.

See :ref:`analysis-output-vibrations` for every table and file schema. The main
entry point is :func:`PQAnalysis.analysis.vibrational.api.vibrations`; direct
calculations use
:func:`PQAnalysis.analysis.vibrational.vibrational_analysis.calculate_from_system`.
