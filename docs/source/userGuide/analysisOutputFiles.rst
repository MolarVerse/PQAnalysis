.. _analysisOutputFiles:

Analysis Output Files
=====================

PQAnalysis analysis commands can write native text, CSV, TSV or XVG tables. The
output filename selects the format:

.. list-table:: Analysis output formats
   :header-rows: 1
   :widths: 18 30 52

   * - Extension
     - Format
     - Intended use
   * - ``.csv``
     - Comma-separated values
     - Excel, LibreOffice Calc, pandas and other table tools
   * - ``.tsv``
     - Tab-separated values
     - Spreadsheets and tables whose values may later contain commas
   * - ``.xvg``
     - Grace XY data with plot directives
     - Direct inspection with xmgrace
   * - Any other extension
     - Native PQAnalysis text
     - Self-describing scientific data and legacy workflows

This means that ``out_file = table.csv`` in an RDF, MSD, VACF or vibrations input
file writes CSV directly. Names ending in ``.dat``, ``.out``, ``.txt`` or no
extension retain the native format.

Native format
=============

Column numbers in this reference are one-based. Every native tabular output
starts with a compact UTF-8 ``#`` metadata block. ``FIELDS`` lists one stable
ASCII identifier per numeric column, ``SYMBOLS`` gives the corresponding
scientific notation in Unicode, and ``UNITS`` lists the units in the same
order. Values remain single whitespace-free tokens; compound units use a
middle dot. Examples include ``ν̃``, ``Å``, ``Å³`` and ``cm⁻¹``. The tables
below give the full quantities and definitions. Numeric rows retain the legacy
ordering and formatting. Readers such as ``numpy.loadtxt`` ignore the comment
block automatically and continue to work without special options.

For example, an RDF data file begins with

.. code-block:: text

   # PQAnalysis: Radial distribution function
   # FIELDS r_i g_r_i N_r_i g_r_i_dV_i H_i_minus_E_i
   # SYMBOLS rᵢ g(rᵢ) N(rᵢ) g(rᵢ)ΔVᵢ Hᵢ−Eᵢ
   # UNITS Å 1 1 Å³ pairs
   0.5 0.0 0.0 0.0 -0.05026548245743666

Log files are labeled, human-readable text and are not described as columnar
data here.

CSV and TSV
===========

CSV and TSV files contain one header row followed by numeric rows. The header
uses the stable identifiers from the native ``FIELDS`` line, for example:

.. code-block:: text

   r_i,g_r_i,N_r_i,g_r_i_dV_i,H_i_minus_E_i
   0.5,0.0,0.0,0.0,-0.05026548245743666

These files open directly in Excel and LibreOffice Calc and can be read as
ordinary CSV or TSV. Scientific symbols, units and definitions are documented
by the field identifier in the tables below.

XVG and xmgrace
===============

XVG files contain Grace ``@`` directives and standard XY data-set blocks. They
can be opened directly with ``xmgrace rdf.xvg``; PQAnalysis does not launch the
GUI itself. Each original analysis-table column is stored as one Grace data set
with the selected x axis. Columns outside the selected quick plot are retained
as hidden Grace sets. PQAnalysis metadata records the original schema and plot
projection, making its XVG output fully convertible back to native, CSV, TSV or
XVG without losing analysis columns. The default quick plots are:

.. list-table:: Default XVG plots
   :header-rows: 1
   :widths: 34 28 38

   * - Analysis table
     - x axis
     - y data sets
   * - RDF
     - ``r_i``
     - ``g_r_i``
   * - MSD
     - ``lag``
     - ``msd_x``, ``msd_y``, ``msd_z``
   * - VACF or windowed correlation
     - ``time``
     - Correlation value
   * - VACF or broadened spectrum
     - ``wavenumber``
     - Amplitude or intensity
   * - Total momentum
     - ``frame``
     - ``scaled_momentum_norm``
   * - Vibrations without IR intensities
     - One-based mode index
     - ``wavenumber``
   * - Vibrations with IR intensities
     - ``wavenumber``
     - ``ir_intensity``

Additional outputs
==================

All tabular analysis CLIs accept repeatable ``--export FILE`` options. The
primary output is still controlled by ``out_file`` or ``--output``; every
export filename independently selects its format.

.. code-block:: console

   $ pqanalysis rdf rdf.in \
       --export rdf.csv \
       --export rdf.tsv \
       --export rdf.xvg

The same option is available for ``msd``, ``vacf``, ``build_spectrum``,
``check_momentum`` and ``vibrations``. For VACF, ``--export`` represents the
main correlation table. The configured ``spectrum_file`` and
``windowed_out_file`` select their own formats from their filenames.

Converting existing output
==========================

The converter reads native, CSV, TSV or PQAnalysis-generated XVG analysis
tables and can create several outputs in one command:

.. code-block:: console

   $ pqanalysis convert rdf.dat \
       -o rdf.csv \
       -o rdf.tsv \
       -o rdf.xvg

The conversion also works in the other direction:

.. code-block:: console

   $ pqanalysis convert rdf.xvg \
       -o rdf.out \
       -o rdf.csv

Input format is detected from the file content rather than its extension, so a
CSV table named ``table.dat`` or an XVG table named ``table.out`` can still be
converted. Exact ``FIELDS`` headers restore the known scientific schema and XVG
plot preset. ``--x FIELD`` and repeatable ``--y FIELD`` options override the
default XVG projection without removing unplotted data from the XVG file. Input
and output paths, and all output paths, must be distinct.

By default, conversion stops with an error naming the existing file if any
output path already exists. No requested output is written in that case. Use
``--mode o`` only to request replacement explicitly.

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
   :class: analysis-output-columns
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
   :class: analysis-output-columns
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
   :class: analysis-output-columns
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
   :class: analysis-output-columns
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
   :class: analysis-output-columns
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
   :class: analysis-output-columns
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
   :class: analysis-output-columns
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
The accompanying ``SYMBOLS`` line provides the Unicode scientific notation.

.. list-table:: Vibrational ``out_file`` columns
   :class: analysis-output-columns
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
