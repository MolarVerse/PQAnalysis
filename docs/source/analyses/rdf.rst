Radial Distribution Function
============================

The radial distribution function measures the probability of finding a target
atom at distance :math:`r` from a reference atom relative to an ideal gas at
the same effective target density. For histogram bin :math:`i`, PQAnalysis
uses

.. math::

   g_i = \frac{H_i}{\rho_T N_R N_F \Delta V_i},

where :math:`H_i` is the eligible pair count, :math:`\rho_T` the target number
density, :math:`N_R` the number of reference atoms, :math:`N_F` the number of
frames and :math:`\Delta V_i` the spherical-shell volume.

Minimal input
-------------

.. code-block:: text

   traj_files = trajectory.xyz
   reference_selection = O
   target_selection = H
   delta_r = 0.05
   r_max = 8.0
   out_file = rdf.dat

.. code-block:: console

   $ pqanalysis rdf rdf.in

``restart_file`` and ``moldescriptor_file`` are unnecessary for a basic
species RDF. They are required when ``no_intra_molecular = True`` is used to
exclude pairs belonging to the same molecule. PQAnalysis can infer the usual
PQ companion filenames when they are beside the trajectory.

Interpretation
--------------

* Peaks mark preferred pair separations; minima separate coordination shells.
* :math:`g(r) \approx 1` indicates bulk-like, uncorrelated pair density at that
  distance.
* The cumulative coordination column gives the mean number of eligible target
  atoms per reference atom inside the current upper bin edge.
* Self pairs are excluded. Intramolecular pairs are included unless molecular
  topology is supplied and explicitly excluded.

Normalization, finite-size effects, selection definitions and trajectory
sampling should be considered before comparing RDFs from different systems.

Output and API
--------------

See :ref:`analysis-output-rdf` for the five output columns and their exact
normalization. The main Python entry point is
:func:`PQAnalysis.analysis.rdf.api.rdf`; lower-level calculations use
:class:`PQAnalysis.analysis.rdf.rdf.RDF`.

The complete input-key table is documented with
:class:`PQAnalysis.analysis.rdf.rdf_input_file_reader.RDFInputFileReader`.
