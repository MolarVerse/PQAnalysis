Radial Distribution Function
============================

The radial distribution function measures the probability of finding a target
atom at distance :math:`r` from a reference atom relative to an ideal gas at
the same effective target density [Hansen2013]_. For histogram bin :math:`i`,
PQAnalysis uses the standard simulation estimator [Allen2017]_,

.. math::

   g_i = \frac{H_i}{\rho_T N_R N_F \Delta V_i},

where :math:`H_i` is the eligible pair count, :math:`\rho_T` the target number
density, :math:`N_R` the number of reference atoms, :math:`N_F` the number of
frames and :math:`\Delta V_i` the spherical-shell volume. Peaks mark preferred
pair separations, minima separate coordination shells, and
:math:`g(r) \approx 1` is uncorrelated bulk-like pair density.

.. plot:: _plots/rdf.py
   :alt: Radial distribution function and cumulative coordination number
   :caption: Analytic schematic rather than simulation output. The shaded
      interval ends at the first minimum. The lower panel evaluates
      :math:`N(r) = 4\pi\rho\int_0^r g(s)\,s^2\,\mathrm{d}s` with
      :math:`\rho = 0.0334` Å⁻³.

Input
-----

.. code-block:: text

   traj_files = trajectory.xyz
   reference_selection = O
   target_selection = H
   delta_r = 0.05
   r_max = 8.0
   out_file = rdf.dat

Saved as ``rdf.in``, it runs with:

.. code-block:: console

   $ pqanalysis rdf rdf.in

The keys above are typical for a bulk trajectory. The bundled
:doc:`../examples` fixture uses ``delta_r = 0.5`` and ``r_max = 4.0`` because
it is one molecule in a 10 Å box. ``restart_file`` and ``moldescriptor_file``
are needed only with ``no_intra_molecular = True``, which drops pairs inside
the same molecule; PQAnalysis infers the usual PQ companion filenames when
they sit beside the trajectory. Selection strings are described in
:doc:`../data/selections`; the full key table is on
:class:`~PQAnalysis.analysis.rdf.rdf_input_file_reader.RDFInputFileReader`.

Output
------

``out_file`` has five columns: bin center, :math:`g(r)`, running coordination
number, shell population and pair-count residual. Exact definitions are in
:ref:`analysis-output-rdf`. The coordination number is read off at the row
where :math:`g` has its first minimum.

Python
------

.. code-block:: python

   from PQAnalysis.analysis import RDF
   from PQAnalysis.io import TrajectoryReader

   r, g, n, shell, residual = RDF(
       TrajectoryReader("examples/water/trajectory.xyz"),
       reference_species="O",
       target_species="H",
       delta_r=0.5,
       r_max=4.0,
   ).run()

:func:`~PQAnalysis.analysis.rdf.api.rdf` runs an input file instead.

Before you trust it
-------------------

* ``r_max`` above half the shortest box edge is clamped; on triclinic cells
  set it to half the smallest perpendicular width yourself.
* Bins near the first peak need about :math:`10^4` pair counts for 1 % noise.
  Add frames rather than widening ``delta_r``.
* The plateau at large :math:`r` should sit at 1. If it does not, look at the
  volume distribution (NPT) or the selections.
* Intramolecular pairs are included unless excluded explicitly.

Each point is derived in :doc:`rdf-details`.

.. toctree::
   :hidden:

   rdf-details
