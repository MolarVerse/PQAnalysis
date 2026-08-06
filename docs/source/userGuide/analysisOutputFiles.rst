.. _analysisOutputFiles:

#####################
Analysis Output Files
#####################

PQAnalysis analysis commands write whitespace-separated text files. The
filename extension does not select or change the analysis output format;
``.out``, ``.dat``, ``.txt`` and extensionless names are all accepted. Column
numbers in this reference are one-based. Every tabular output starts with a
compact ``#`` metadata block. ``FIELDS`` lists one whitespace-free scientific
symbol per numeric column and ``UNITS`` lists the corresponding units in the
same order. Compound units use ``*`` so each unit remains one token. Symbols
use ASCII spellings such as ``nu_tilde`` and ``DeltaV`` to keep the files
portable; the tables below give their full quantities and definitions. Numeric
rows retain the legacy ordering and formatting. Readers such as
``numpy.loadtxt`` ignore the comment block automatically and continue to work
without special options.

Log files are labeled, human-readable text and are not described as columnar
data here.

.. _analysis-output-rdf:

RDF
===

The ``rdf`` command writes five columns to ``out_file``. For bin :math:`i`, let
:math:`H_i` be the number of eligible reference-target pairs accumulated over
all frames, :math:`N_R` the number of reference atoms, :math:`N_F` the number
of frames and :math:`\rho_T` the effective target number density. The bin edges
and spherical-shell volume are

.. math::

   r_i^- = r_{\min} + i\,\Delta r, \qquad
   r_i^+ = r_i^- + \Delta r, \qquad
   \Delta V_i = \frac{4\pi}{3}\left((r_i^+)^3 - (r_i^-)^3\right).

The ideal-gas pair count for the shell is
:math:`E_i = \rho_T N_R N_F \Delta V_i`.

.. list-table:: RDF ``out_file`` columns
   :header-rows: 1
   :widths: 8 24 48 20

   * - Column
     - Quantity
     - Definition
     - Unit
   * - 1
     - Bin-center distance
     - :math:`r_i = (r_i^- + r_i^+) / 2`
     - Angstrom
   * - 2
     - Radial distribution function
     - :math:`g_i = H_i / E_i`
     - Dimensionless
   * - 3
     - Cumulative coordination number
     - :math:`\sum_{j=0}^{i} H_j / (N_R N_F)`, the average number of
       eligible target atoms per reference atom up to :math:`r_i^+`
     - Dimensionless
   * - 4
     - Density-normalized shell population
     - :math:`H_i / (\rho_T N_R N_F) = g_i\Delta V_i`
     - Angstrom\ :sup:`3`
   * - 5
     - Ideal-gas pair-count residual
     - :math:`H_i - E_i`; positive values are an excess and negative values
       are a deficit relative to the ideal-gas shell count
     - Pair count

Self pairs are excluded. When ``no_intra_molecular`` is enabled, pairs from the
same molecule are excluded from :math:`H_i` and from the effective target
density.

.. _analysis-output-msd:

MSD
===

The ``msd`` command writes the legacy Diffcalc layout to ``out_file``.

.. list-table:: MSD ``out_file`` columns
   :header-rows: 1
   :widths: 8 32 40 20

   * - Column
     - Quantity
     - Definition
     - Unit
   * - 1
     - Lag index
     - Frame lag from zero through ``window``
     - Frames
   * - 2
     - :math:`\mathrm{MSD}_x`
     - Mean squared displacement along x
     - Angstrom\ :sup:`2`
   * - 3
     - :math:`\mathrm{MSD}_y`
     - Mean squared displacement along y
     - Angstrom\ :sup:`2`
   * - 4
     - :math:`\mathrm{MSD}_z`
     - Mean squared displacement along z
     - Angstrom\ :sup:`2`

The total MSD is the sum of columns 2 through 4. It is returned by the Python
API but is not repeated in the file. If ``time_step`` is provided, multiply
column 1 by it to obtain lag time in ps. The fitted diffusion coefficients are
written to ``log_file``, not ``out_file``.

.. _analysis-output-vacf:

VACF and charge-flux correlation
================================

The ``vacf`` command can write three two-column files. Correlations are
normalized by their zero-lag value, including charge-weighted correlations.

.. list-table:: VACF ``out_file`` columns
   :header-rows: 1
   :widths: 8 42 30 20

   * - Column
     - Quantity
     - Definition
     - Unit
   * - 1
     - Lag time
     - Frame lag multiplied by ``time_step``
     - ps
   * - 2
     - Normalized correlation
     - VACF, or charge-flux autocorrelation when charges are supplied
     - Dimensionless

.. list-table:: VACF ``spectrum_file`` columns
   :header-rows: 1
   :widths: 8 42 30 20

   * - Column
     - Quantity
     - Definition
     - Unit
   * - 1
     - Wavenumber
     - Legacy cosine-transform frequency axis
     - cm\ :sup:`-1`
   * - 2
     - Spectrum amplitude
     - Absolute cosine-transform amplitude of the optionally windowed
       normalized correlation
     - Arbitrary units

.. list-table:: VACF ``windowed_out_file`` columns
   :header-rows: 1
   :widths: 8 42 30 20

   * - Column
     - Quantity
     - Definition
     - Unit
   * - 1
     - Lag time
     - Same axis as ``out_file``
     - ps
   * - 2
     - Windowed correlation
     - Column 2 of ``out_file`` multiplied by the selected apodization window
     - Dimensionless

.. _analysis-output-spectrum:

Broadened spectrum
==================

The ``build_spectrum`` command writes two columns to ``--output`` or standard
output.

.. list-table:: ``build_spectrum`` output columns
   :header-rows: 1
   :widths: 8 42 30 20

   * - Column
     - Quantity
     - Definition
     - Unit
   * - 1
     - Wavenumber
     - Regular output grid from ``--min`` to the exclusive ``--max``
     - cm\ :sup:`-1`
   * - 2
     - Broadened intensity
     - Sum of the Gaussian or Lorentzian peak-height profiles
     - Same as the input intensities

The broadening uses a peak-height convention and does not normalize peak area.

.. _analysis-output-momentum:

Total linear momentum
=====================

The ``check_momentum`` command writes two columns to ``--output`` or standard
output.

.. list-table:: ``check_momentum`` output columns
   :header-rows: 1
   :widths: 8 42 30 20

   * - Column
     - Quantity
     - Definition
     - Unit
   * - 1
     - Frame index
     - One-based index across all input trajectory files
     - Frames
   * - 2
     - Scaled momentum norm
     - ``scale`` multiplied by :math:`\left|\sum_i m_i\mathbf{v}_i\right|`
     - Set by ``--scale``; default is amu Angstrom fs\ :sup:`-1`

.. _analysis-output-vibrations:

Vibrational analysis
====================

The ``vibrations`` command writes one required table and up to three optional
normal-mode representations.

``out_file``
------------

Without a ``moldescriptor_file``, the table has three columns. With partial
charges, the IR-intensity column is inserted as column 2 and the table has four
columns. The file's ``FIELDS`` and ``UNITS`` lines reflect the selected layout.

.. list-table:: Vibrational ``out_file`` columns
   :header-rows: 1
   :widths: 8 34 38 20

   * - Column
     - Quantity
     - Availability
     - Unit
   * - 1
     - Signed wavenumber
     - Always; negative values represent imaginary modes
     - cm\ :sup:`-1`
   * - 2
     - IR intensity
     - Only with partial charges
     - km mol\ :sup:`-1`
   * - 2 or 3
     - Force constant
     - Always
     - mdyn Angstrom\ :sup:`-1`
   * - 3 or 4
     - Reduced mass
     - Always
     - amu

``normal_modes_file``
---------------------

This file contains the normal-mode matrix. Its metadata identifies the matrix
element as :math:`e_{\alpha j}`, the row order as ``x_1,y_1,z_1,...``, and each
column as ``mode_1``, ``mode_2`` and so on. Columns are modes in the same order
as the rows of ``out_file``. Values are dimensionless normalized Cartesian mode
components.

``modes_prefix``
----------------

One multi-frame XYZ animation named ``<prefix>-<mode>.xyz`` is written per
selected mode. Each atom row contains species, x, y and z in Angstrom. The XYZ
comment records the one-based mode number, wavenumber in cm\ :sup:`-1`, frame
number and sinusoidal phase.

``modes_file``
--------------

This extended XYZ file contains one image per selected mode. The atom columns
are species, equilibrium x/y/z coordinates in Angstrom and normalized mode
x/y/z components. The comment declares
``Properties=species:S:1:pos:R:3:mode:R:3`` and records the one-based mode
number, wavenumber and optional IR intensity.
