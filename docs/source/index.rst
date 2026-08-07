PQAnalysis
==========

PQAnalysis provides command-line and Python tools for quantitative analysis of
PQ molecular-dynamics simulations. It reads structures, trajectories,
velocities and Hessians, then produces documented scientific tables for
structural, transport and vibrational observables.

:doc:`Get started <getting-started>` | :doc:`Choose an analysis <analyses/index>` |
:doc:`Work with data <data/index>` | :doc:`Command reference <reference/cli>`

Quick start
-----------

PQAnalysis requires Python 3.12 or newer.

.. code-block:: console

   $ python -m pip install pqanalysis
   $ pqanalysis rdf rdf.in

The output filename in an analysis input file selects native text, CSV, TSV or
XVG. Repeat ``--export`` to write several formats in the same run.

.. code-block:: console

   $ pqanalysis rdf rdf.in --export rdf.csv --export rdf.xvg

Analysis methods
----------------

.. list-table:: Implemented observables
   :header-rows: 1
   :widths: 24 38 38

   * - Method
     - Required data
     - Reported quantity
   * - :doc:`Radial distribution <analyses/rdf>`
     - Positions and periodic cell
     - :math:`g_{AB}(r)` and cumulative coordination
   * - :doc:`Mean square displacement <analyses/msd>`
     - Positions and periodic cell
     - Cartesian MSD and diffusion fits
   * - :doc:`VACF and spectra <analyses/vacf>`
     - Velocities, sampling interval and optional charges
     - Normalized correlation and wavenumber spectrum
   * - :doc:`Vibrational analysis <analyses/vibrations>`
     - Structure, masses and Cartesian Hessian
     - Normal modes, wavenumbers and optional IR intensities
   * - :doc:`Momentum diagnostic <analyses/momentum>`
     - Velocities and atomic masses
     - Frame-resolved total linear momentum

Scientific output
-----------------

Native analysis tables retain their established numeric layout and add a
compact UTF-8 metadata header. Stable ASCII field names support scripts while
Unicode symbols and units describe the physical quantities.

.. code-block:: text

   # PQAnalysis: Radial distribution function
   # FIELDS r_i g_r_i N_r_i g_r_i_dV_i H_i_minus_E_i
   # SYMBOLS rᵢ g(rᵢ) N(rᵢ) g(rᵢ)ΔVᵢ Hᵢ−Eᵢ
   # UNITS Å 1 1 Å³ pairs
   0.5 0.0 0.0 0.0 -0.05026548245743666

See :ref:`analysisOutputFiles` for column definitions, normalization
conventions and conversion behavior.

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Documentation

   getting-started
   analyses/index
   data/index
   reference/index
   Development <developerGuide/developerGuide>
