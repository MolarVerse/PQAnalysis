Vibrational Analysis
====================

Vibrational analysis diagonalizes the mass-weighted Cartesian Hessian of a
single structure [Wilson1955]_. Its eigenvectors are the harmonic normal modes;
its eigenvalues give signed wavenumbers, force constants and reduced masses.
When partial charges are supplied, point-charge infrared intensities are
reported as well. The whole calculation is the harmonic approximation applied
to one isolated structure: no dynamics, temperature or anharmonicity enters it.

.. plot:: _plots/vibrations.py
   :alt: Infrared stick spectrum for the water validation fixture
   :caption: IR stick spectrum calculated by PQAnalysis from the bundled H₂O
      structure, Hessian and partial-charge fixtures. Only internal modes above
      100 cm⁻¹ are shown; translational and rotational modes are omitted.

Input
-----

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
names the Hessian energy unit: ``kcal`` and ``ev`` expect Å⁻², ``hartree``
expects bohr⁻². ``hessian_sign`` is ``positive``, ``negative`` or ``auto``
(majority vote over internal curvatures). IR intensities are written only when
a ``moldescriptor_file`` supplies partial charges. The bundled
:doc:`../examples` fixture is ``examples/water``; the full key table is on
:class:`~PQAnalysis.analysis.vibrational.vibrational_input_file_reader.VibrationalAnalysisInputFileReader`.

``modes`` selects which modes reach ``modes_prefix`` and ``modes_file``,
compared against ``modes_threshold`` :math:`\theta` (default
:math:`10^{-8}` cm⁻¹):

.. list-table:: ``modes`` selection rules
   :class: pq-record-table
   :header-rows: 1
   :widths: 24 76

   * - Value
     - Selected modes
   * - ``all`` (default)
     - Every mode, in order of increasing wavenumber
   * - ``nonzero``
     - :math:`\lvert\tilde{\nu}_j\rvert > \theta`, keeping imaginary modes
   * - ``positive``
     - :math:`\tilde{\nu}_j > \theta`, dropping imaginary modes
   * - Numbers
     - Explicit one-based mode numbers: an integer, a list or a range

Raise ``modes_threshold`` to a few cm⁻¹ when the intent is "internal modes
only"; at the default, ``positive`` still keeps near-zero translations.

Output
------

``out_file`` lists wavenumber (negative means imaginary), optional IR
intensity, reduced mass and force constant per mode
(:ref:`analysis-output-vibrations`). ``normal_modes_file`` stores the
dimensionless Cartesian mode matrix, ``modes_file`` one extended-XYZ image per
selected mode, and ``modes_prefix`` one sinusoidal animation per mode scaled by
``modes_amplitude`` (default 0.25 Å) or ``modes_temperature``.

Python
------

.. code-block:: python

   from PQAnalysis.analysis import vibrations, read_analysis_table
   from PQAnalysis.analysis.vibrational import calculate, read_hessian_file
   from PQAnalysis.io import read_restart_file

   vibrations("examples/water/vibrations.in", export_files=["wavenumbers.csv"])
   table = read_analysis_table("wavenumbers.csv")
   print(table.column("wavenumber")[:10])

   system = read_restart_file("examples/water/structure.rst")
   hessian = read_hessian_file("examples/water/hessian.dat")
   result = calculate(
       system.atomic_masses,
       system.pos,
       hessian,
       unit="kcal",
   )
   print(result.wavenumbers[:10])
   print(result.force_constants[:10])

Before you trust it
-------------------

* Expect :math:`3N` modes: six external ones (five for linear molecules) near
  zero, the rest positive at a minimum. Residual external modes of tens of
  cm⁻¹ mean an unconverged geometry or a noisy Hessian.
* A wrong ``unit`` or atom ordering yields a plausible spectrum on the wrong
  scale without an error. Set ``hessian_sign`` explicitly when the producing
  code's convention is known.
* Harmonic wavenumbers exceed observed fundamentals; no scaling factor is
  applied. Point-charge IR intensities are qualitative.
* The external-mode construction assumes an isolated molecule.

Each point is derived in :doc:`vibrations-details`.

.. toctree::
   :hidden:

   vibrations-details
