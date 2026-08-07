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

Documentation
-------------

.. grid:: 1 2 3 3
   :gutter: 2

   .. grid-item-card:: Getting started
      :link: getting-started
      :link-type: doc

      Install PQAnalysis and run a first radial-distribution calculation.

   .. grid-item-card:: Analyses
      :link: analyses/index
      :link-type: doc

      RDF, MSD, VACF, spectra, normal modes and momentum diagnostics.

   .. grid-item-card:: Data and conversion
      :link: data/index
      :link-type: doc

      Input syntax, trajectories, selections and scientific table formats.

   .. grid-item-card:: Command line
      :link: reference/cli
      :link-type: doc

      Analysis, conversion and trajectory command reference.

   .. grid-item-card:: Python API
      :link: reference/api
      :link-type: doc

      Curated entry points and the complete generated package reference.

   .. grid-item-card:: Development
      :link: developerGuide/developerGuide
      :link-type: doc

      Branching, tests, documentation checks and contribution conventions.

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

See :ref:`analysisOutputFiles` for every column, normalization convention and
conversion path.

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Documentation

   getting-started
   analyses/index
   data/index
   reference/index
   Development <developerGuide/developerGuide>
