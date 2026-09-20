MSD: Theory and Validity
========================

Background for :doc:`msd`: when the Einstein fit is legitimate, how to choose
the fit window, how much averaging the curve carries, what the printed
uncertainty means, finite-size effects and the unwrapping and normalization
traps of the legacy-compatible path.

Compatibility path
------------------

File-backed orthorhombic trajectories use a bounded compatibility path that
preserves the operation order of the legacy ``Diffcalc`` program [thhTools]_.
Unsupported cells or inputs return to the general streaming implementation.

Validity
--------

Confirm the diffusive regime before trusting a fit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Einstein relation applies only where the MSD is linear in time. The
diagnostic is the local log-log slope

.. math::

   \beta(t) = \frac{\mathrm{d}\log\mathrm{MSD}_{\mathrm{total}}(t)}
                   {\mathrm{d}\log t},

which equals 2 in the ballistic short-time regime, can sit well below 1 on a
sub-diffusive plateau in cage-forming liquids, glasses and confined systems,
and approaches 1 only once the motion has become diffusive. Fit only over an
interval where :math:`\beta \approx 1`, and quote that interval together with
:math:`D`.

PQAnalysis does not compute :math:`\beta` and performs no linearity test. The
diagnostic has to be evaluated from ``out_file``: multiply column 1 by
``time_step`` to obtain the lag time, and sum columns 2 through 4 to obtain
:math:`\mathrm{MSD}_{\mathrm{total}}`.

.. warning::

   A high :math:`R^2` is not evidence of diffusion. A purely ballistic
   :math:`\mathrm{MSD}\propto t^2` fitted over the trailing points of a window
   still yields :math:`R^2 \approx 0.999` and a finite, entirely meaningless
   :math:`D`. The :math:`R^2` in the log file measures how straight the
   selected points are, not whether they belong to a diffusive regime.

Choosing the fit window
^^^^^^^^^^^^^^^^^^^^^^^

``fit_window`` is the number of *trailing* MSD points used by the fit, so the
fit always ends at the largest lag. With :math:`W` for ``window``, :math:`F`
for ``fit_window`` and :math:`\Delta t` for ``time_step``, the fitted lag range
is

.. math::

   \left[(W - F + 1)\,\Delta t,\; W\,\Delta t\right].

The default is
``max(2, window // 5)``, the last 20 % of the window. Because the fit is
anchored at the end, two consequences follow:

* The ballistic short-time region is excluded by making ``fit_window``
  *smaller*, which moves the start of the fit to later lag times.
* The noisy long-lag tail cannot be excluded through ``fit_window`` at all; it
  is always inside the fit. To move the fit away from the tail, reduce
  ``window`` so that the largest lag is itself still well sampled.

How much averaging the curve actually carries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Time origins spawn every ``gap`` frames up to
``stop_frame = (n_frames - window) // gap * gap``, so their number is

.. math::

   N_{\mathrm{origins}} = \left\lfloor
   \frac{N_{\mathrm{frames}} - W}{G}\right\rfloor,

for ``window`` :math:`= W` and ``gap`` :math:`= G`. It is printed as
``Number of origins`` in the log file. Every origin covers the full
window, so in this implementation *every* lag bin from 0 to ``window`` is
averaged over exactly this many origins. The origin count does not decay with
lag, as it does in estimators that use every frame as a time origin.

That does not make the tail reliable. The origins overlap heavily, and at lag
:math:`\tau` a trajectory of total length :math:`T` contains at most
:math:`T/\tau` statistically independent displacement windows. The effective
sample size still falls as :math:`1/\tau`, which is why the long-lag tail
remains the noisiest part of the curve even though its nominal origin count is
unchanged.

The same expression sets the trajectory length you need. Choosing ``window``
close to the trajectory length starves the entire curve, not only its tail:
1100 frames with ``window = 1000`` and ``gap = 10`` leave 10 origins for every
lag. A defensible :math:`D` requires a trajectory much longer than the longest
fitted lag: that lag must already lie beyond the velocity correlation time so
that the motion is diffusive, and the trajectory must then be long enough to
contain many independent windows of that length.

What the reported uncertainty is, and is not
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``+/-`` value in the log file is the standard error of the fitted slope
returned by ``scipy.stats.linregress``, scaled by the same
:math:`10^{-8}/(2d)` factor as the coefficient itself, with :math:`d = 1` for
the Cartesian components and :math:`d = 3` for the total. It is the statistical
error of a straight-line fit, computed as if the fitted MSD points were
independent samples.

They are not. Neighbouring lag bins share time origins and atoms and are
strongly correlated, so the quoted error systematically understates the true
uncertainty. PQAnalysis performs no block averaging, no averaging over
independent trajectories and no correction for the correlation between lags.
Treat the printed uncertainty as a lower bound, and obtain a realistic error
bar from independent runs, or at least from the spread of :math:`D` under
variation of the fit window and between the three Cartesian components.

Finite-size effects on D
^^^^^^^^^^^^^^^^^^^^^^^^

Self-diffusion coefficients computed under periodic boundary conditions are
systematically too small, because the periodic images suppress hydrodynamic
backflow. The leading correction for a cubic box is

.. math::

   D_0 = D_{\mathrm{PBC}} + \frac{\xi\,k_{\mathrm{B}} T}{6\pi\eta L},
   \qquad \xi = 2.837297,

with shear viscosity :math:`\eta` and box length :math:`L`.

.. important::

   PQAnalysis does **not** apply the Yeh-Hummer correction, or any other
   finite-size correction. The value written to the log file is the raw
   periodic :math:`D_{\mathrm{PBC}}` at the simulated box size, and no
   viscosity enters the code anywhere. Apply the correction externally, or
   extrapolate :math:`D` to :math:`1/L \to 0` across several box sizes.

   Yeh, I.-C.; Hummer, G. System-Size Dependence of Diffusion Coefficients and
   Viscosities from Molecular Dynamics Simulations with Periodic Boundary
   Conditions. *The Journal of Physical Chemistry B* **2004**, *108* (40),
   15873-15879. `doi:10.1021/jp0477147
   <https://doi.org/10.1021/jp0477147>`__

Coordinates and periodic images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Displacements must be free of periodic jumps, and PQAnalysis handles this
itself: each per-frame displacement is folded into the minimum image of the
current cell and accumulated into a running shift, so a trajectory written with
coordinates wrapped into the box yields exactly the same MSD as the
corresponding unwrapped trajectory. No pre-unwrapping step is required.

The scheme is exact only while no selected atom travels further than half the
shortest box vector between two *written* frames. A sparse output stride breaks
that condition silently: the misassigned images truncate every affected
displacement to at most half a box vector, which biases the MSD and therefore
:math:`D` downwards with no warning. For a vacuum cell no unwrapping is applied
and the coordinates are used as they are.

Two normalization traps
^^^^^^^^^^^^^^^^^^^^^^^

.. warning::

   A late start frame (``first_frame`` or ``start`` in the input file,
   ``n_start`` in Python) shrinks the MSD. Earlier frames are read so that the
   unwrapping stays continuous, but they spawn no time origins, while the
   divisor keeps its legacy value ``stop_frame // gap``, which still counts
   them. The whole curve, and therefore :math:`D`, is scaled by the ratio of
   origins that actually spawned to that divisor. With 60 frames,
   ``window = 20`` and ``gap = 5``, a start frame of 20 leaves 5 of 8 counted
   origins and returns 62.5 % of the correct MSD. This legacy Diffcalc
   convention is deliberate. To start later without the bias, truncate the
   trajectory file instead, or rescale the result yourself.

.. warning::

   A trajectory of exactly ``window`` frames with ``gap = 1`` takes the legacy
   single-origin branch: one origin spawns, and the final lag bin
   (``lag = window``) can never be sampled and is written as exactly 0.0. A
   warning is emitted when the analysis is set up. Since the fit always uses trailing points,
   that zero is inside any requested diffusion fit and corrupts it. Use a
   longer trajectory or a smaller ``window``.

References
----------

* [Einstein1905]_ derives the linear growth of the mean square displacement
  that the diffusion fit assumes.
* [Rahman1964]_ is the first molecular-dynamics measurement of single-particle
  displacements and velocity correlations in a liquid.
* [Allen2017]_ and [Frenkel2002]_ cover the multiple-time-origin estimator,
  coordinate unwrapping in a periodic cell and the practical limits of
  extracting :math:`D` from a finite trajectory.
