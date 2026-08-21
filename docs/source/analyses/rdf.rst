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
frames and :math:`\Delta V_i` the spherical-shell volume.

Structural interpretation
-------------------------

.. plot:: _plots/rdf.py
   :alt: Radial distribution function and cumulative coordination number
   :caption: Analytic schematic, not simulation output. The shaded interval
      ends at the first minimum. The lower panel evaluates
      N(r) = 4πρ∫₀ʳ g(s)s² ds with ρ = 0.0334 Å⁻³.

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

Legacy-compatible arithmetic
----------------------------

For a file-backed periodic orthorhombic trajectory, specifying ``delta_r``
alone with the default ``r_min = 0`` selects the legacy-compatible RDF path.
Coordinates are parsed as float64; histogram binning and all five output
columns preserve the corrected operation order of the legacy ``RDF`` C code
[thhTools]_. Explicit ``r_max`` or
``n_bins``, triclinic or vacuum cells, and intramolecular exclusion use the
general PQAnalysis path. The minimal example above sets ``r_max`` explicitly
and therefore uses the general path.

Interpretation
--------------

* Peaks mark preferred pair separations; minima separate coordination shells.
* :math:`g(r) \approx 1` indicates bulk-like, uncorrelated pair density at that
  distance.
* The cumulative coordination column gives the mean number of eligible target
  atoms per reference atom inside the current upper bin edge.
* Self pairs are excluded. Intramolecular pairs are included unless molecular
  topology is supplied and explicitly excluded.

Validity and interpretation
---------------------------

r_max and the minimum-image limit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Distances are evaluated under the minimum-image convention, which is only
meaningful while the sphere of radius :math:`r` fits inside the periodic cell.
Beyond that radius the shell is no longer fully covered by nearest images, the
pair count falls short of the ideal-gas expectation used for normalization, and
:math:`g(r)` is biased low for purely geometric reasons. The hard limit is half
the shortest perpendicular width of the cell, which for an orthorhombic box is
half the shortest box vector,

.. math::

   r_{\max} \le \tfrac{1}{2}\min(a, b, c).

PQAnalysis enforces this on the general path: a requested ``r_max`` — whether
given directly, or implied by ``n_bins`` and ``delta_r`` — that exceeds
:math:`\tfrac{1}{2}\min(a,b,c)` is clamped down to it, with a warning in the
log. The bound is taken over *all* frames, so for a variable cell the smallest
box in the whole trajectory sets the limit. If a run covers a shorter range
than requested, this clamp is why.

.. warning::

   Two cases escape the bound and must be checked by hand.

   **Non-cubic boxes on the legacy path.** With ``delta_r`` alone on a periodic
   orthorhombic trajectory, the bin count is derived from *half the longest*
   box vector, while the legacy kernel discards every pair beyond half the
   *shortest* one. A 10 × 14 × 30 Å box with ``delta_r = 0.1`` therefore
   produces bins out to 14.9 Å, of which everything past 5.0 Å is exactly
   zero — an artifact, not a depletion zone. Ignore all bins beyond
   :math:`\tfrac{1}{2}\min(a,b,c)`, or set ``r_max`` explicitly to take the
   general path.

   **Triclinic cells.** The bound uses box-vector lengths, not the
   perpendicular widths of the cell, and for a skewed cell the inscribed sphere
   is smaller than half the shortest vector. For :math:`a=b=c=10` Å with
   :math:`\gamma = 60^\circ`, the clamp allows :math:`r_{\max} = 5.0` Å while
   the true limit is 4.33 Å; an ideal gas in that cell already shows
   :math:`g(r) \approx 0.92` between 4.33 and 4.7 Å and :math:`0.79` between
   4.7 and 5.0 Å. On triclinic cells, set ``r_max`` to half the smallest
   perpendicular width yourself.

Sampling: how smooth is smooth enough
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:math:`H_i` is a raw integer pair count, so a bin carries a relative counting
noise of roughly :math:`1/\sqrt{H_i}`. Reaching 1 % noise in a bin needs about
:math:`10^4` counts in that bin. The expected count grows with the shell
volume,

.. math::

   \langle H_i\rangle = \rho_T N_R N_F \Delta V_i \propto r_i^2\,\Delta r,

so the small-:math:`r` bins are always the poorest and a first peak looks
ragged long before the plateau does. The three levers are the number of frames
:math:`N_F`, the size of the reference and target selections, and ``delta_r``:
halving ``delta_r`` halves the counts per bin and raises the relative noise by
:math:`\sqrt{2}`. Choose
``delta_r`` fine enough to locate the first peak and the first minimum — 0.02
to 0.05 Å is typical — and then buy the smoothness back with frames, not by
widening bins. Curves that still wobble around 1 in the plateau region are not
converged, whatever the first peak looks like.

Coordination numbers and the first minimum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first minimum of :math:`g(r)` is the conventional boundary of the first
coordination shell: it is where the shell population is lowest, so the
coordination number is least sensitive to exactly where the cut is placed.
Column 3 of the output is the running coordination number

.. math::

   N(r_i) = \frac{1}{N_R N_F}\sum_{j \le i} H_j,

the mean number of eligible target atoms per reference atom. Read it off at the
row whose :math:`g` is minimal. Two details matter when quoting the number:

* Column 1 is the bin *center* while :math:`N(r_i)` accumulates through the
  *upper* edge of that bin, so the quoted radius and the integration limit
  differ by :math:`\Delta r / 2`.
* A shallow or ill-defined minimum means the coordination number is not
  well-defined either. Quote the cutoff radius alongside the number, and check
  how much it changes when the cutoff moves by one bin.

Density normalization uses the average box volume
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The target density is :math:`\rho_T = N_T / \langle V\rangle`, where
:math:`\langle V\rangle` is the arithmetic mean of the per-frame cell volumes
over the whole trajectory. One density is used for all frames and all bins.

For NVT, NVE and any other fixed-cell trajectory this is exact. For NPT it is
not: the correct normalization would divide each frame by its own volume, and
using a single mean instead scales the whole curve by
:math:`\langle V\rangle\langle 1/V\rangle`, which is greater than one whenever
the volume fluctuates.

.. important::

   The bias is uniform in :math:`r`, so it moves the plateau away from 1
   without changing peak positions. A deliberately extreme test — an ideal gas
   in a trajectory alternating between 1000 Å³ and 2197 Å³ — returns
   :math:`g(r) \approx 1.16` where the exact answer is 1.00, matching
   :math:`\langle V\rangle\langle 1/V\rangle = 1.163`.

   Ordinary NPT volume fluctuations are far smaller and the effect is usually
   negligible, but it is worth checking: the plateau of :math:`g(r)` at large
   :math:`r` should sit at 1. If it does not, and the sampling is converged,
   the volume distribution is the first thing to look at. Note also that
   ``r_max`` is bounded by the *smallest* box in the trajectory, so a
   fluctuating cell shortens the usable range.

Selections, exclusions and comparability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Self pairs are always excluded. Intramolecular pairs are *included* by default,
which for a molecular liquid puts intramolecular bond and angle distances into
the first bins of the histogram; they are real distances, but they are not the
intermolecular structure most RDFs are meant to show. Set
``no_intra_molecular = True``, with the topology files it requires, to remove
them.

.. note::

   With ``no_intra_molecular = True`` the effective target density is derived
   from the exclusion list of the *first* reference atom only. The
   normalization is therefore correct as long as every reference atom excludes
   the same number of target atoms, which holds for a single molecular species
   and fails for a reference selection spanning molecules of different sizes.

Because :math:`g(r)` is normalized to an effective density defined by the
selections and the average volume, two RDFs are only comparable when the
selections, the exclusion setting, the thermodynamic state and the sampled
:math:`r` range all match. State them with the curve.

Output and API
--------------

See :ref:`analysis-output-rdf` for the five output columns and their exact
normalization. The main Python entry point is
:func:`PQAnalysis.analysis.rdf.api.rdf`; lower-level calculations use
:class:`PQAnalysis.analysis.rdf.rdf.RDF`.

The complete input-key table is documented with
:class:`PQAnalysis.analysis.rdf.rdf_input_file_reader.RDFInputFileReader`.

References
----------

* [Hansen2013]_ defines :math:`g(r)` in classical liquid-state theory and
  relates it to the structure factor and to thermodynamic averages.
* [Allen2017]_ and [Frenkel2002]_ give the histogram estimator, its
  spherical-shell normalization and the finite-size caveats that apply to a
  periodic simulation cell.
* [thhTools]_ is the legacy program whose operation order the
  legacy-compatible path reproduces.

Full entries are listed in :doc:`../references`.
