VACF and Spectra
================

The normalized velocity autocorrelation function describes how rapidly atomic
velocities lose memory of their initial direction [Rahman1964]_, [Allen2017]_:

.. math::

   C_{vv}(t) =
   \left\langle
   \frac{\sum_i \mathbf{v}_i(t_0)\cdot\mathbf{v}_i(t_0+t)}
        {\sum_i \mathbf{v}_i(t_0)\cdot\mathbf{v}_i(t_0)}
   \right\rangle_{t_0}.

The brackets denote an average over admissible time origins, and the sum runs
over the selected atoms. This is the default, legacy-compatible estimator
[thhTools]_ and gives :math:`C_{vv}(0)=1`. The ``fft`` estimator instead
averages the numerator and denominator separately over all available origins
before normalization, evaluating the correlation through the power spectrum as
the Wiener-Khinchin theorem allows [Wiener1930]_, [Khintchine1934]_.

PQAnalysis can transform the correlation to a wavenumber-domain spectrum. That
transform of the velocity autocorrelation function is the vibrational density
of states [Dickey1969]_, [Thomas2013]_. If static or time-dependent partial
charges are supplied, it correlates :math:`q_i\mathbf{v}_i` instead, producing
a charge-flux spectrum that approximates an infrared spectrum [Thomas2013]_.

The correlation written to ``out_file`` is not apodized. When a spectrum is
requested, ``window_function`` multiplies a copy of the correlation before the
cosine transform [Harris1978]_. The optional ``windowed_out_file`` records that
copy.

Correlation and spectrum
------------------------

.. plot:: _plots/vacf.py
   :alt: Analytical normalized VACF, its exponentially windowed copy and the
      resulting spectrum
   :caption: Analytical normalized VACF for two Gaussian-broadened bands
      centered at 300 and 600 cm⁻¹, with dephasing times of 0.22 and 0.12 ps.
      The dashed curve applies an exponential window with a decay coefficient
      of 4 ps⁻¹ before the PQAnalysis cosine transform. Spectrum amplitudes
      are scaled to unit maximum.

Minimal input
-------------

.. code-block:: text

   traj_files = trajectory.vel
   target_selection = all
   out_file = vacf.dat
   time_step = 0.001
   window = 2500
   gap = 5
   spectrum_file = spectrum.dat
   ftsize = 5000
   window_function = exponential
   window_param = 4.0

.. code-block:: console

   $ pqanalysis vacf vacf.in

The time step is specified in ps. ``window`` is the maximum correlation lag in
frames; it is distinct from the apodization selected by ``window_function``.
The example multiplies only the spectrum input by
:math:`\exp[-(4\ \mathrm{ps}^{-1})t]`. ``window_function`` also accepts
``hann``, ``blackman`` and ``none``; ``none`` is the default. The default
sliding-origin method matches the legacy calculation. ``gap`` controls the
spacing between its time origins, not the lag-time spacing in ``out_file``;
``method = fft`` selects a denser-origin Wiener-Khinchin estimator.

Interpretation
--------------

* A rapidly decaying VACF indicates fast velocity decorrelation.
* In liquids, negative regions often indicate backscattering or cage motion;
  in solids, sign oscillations reflect bound vibrational motion.
* Charge-flux spectra require physically meaningful partial charges and should
  not be interpreted as absolute IR intensities without further calibration.

Validity and interpretation
---------------------------

The frequency axis, and what sets it
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two different parameters control the two different properties of the spectrum,
and they are easy to confuse. Write :math:`\Delta t` for ``time_step`` in ps,
:math:`W` for ``window`` in frames and :math:`F` for ``ftsize``.

**The grid spacing is set by** ``ftsize``. PQAnalysis mirrors the padded
correlation into an even extension and labels the transform with

.. math::

   \Delta\tilde\nu = \frac{1}{2(F-1)\,\Delta t\,c}
   \approx \frac{16.68}{(F-1)\,\Delta t[\mathrm{ps}]}\ \mathrm{cm}^{-1}.

The axis starts at :math:`\Delta\tilde\nu` — there is no :math:`\tilde\nu = 0`
point — and runs up to :math:`F\,\Delta\tilde\nu`. For ``time_step = 0.001``
and ``ftsize = 5000`` that is a 3.34 cm⁻¹ grid reaching 16682 cm⁻¹.

**The upper limit is the Nyquist wavenumber**, set by the sampling interval
alone:

.. math::

   \tilde\nu_{\max} = \frac{1}{2\,\Delta t\,c}
   \approx \frac{16.68}{\Delta t[\mathrm{ps}]}\ \mathrm{cm}^{-1},

which is 16678 cm⁻¹ for frames written every 1 fs and 1668 cm⁻¹ for frames
written every 10 fs. What matters is the interval between *written* frames, not
the integration time step of the underlying dynamics: an X-H stretch near
3000 cm⁻¹ needs a trajectory stride below about 5.5 fs to be represented at
all, and motion above the limit is aliased back into the spectrum rather than
discarded. The legacy axis is one grid point longer than Nyquist,
:math:`F/(F-1)\cdot\tilde\nu_{\max}`, and the last spectrum point duplicates
its predecessor; both are deliberate legacy conventions.

.. note::

   The axis carries a deliberate legacy calibration: the spacing is computed
   with a period of :math:`2(F-1)` points although the underlying even
   extension has :math:`2F-1` points. Reported wavenumbers are therefore high
   by a factor :math:`(2F-1)/(2F-2)`. The absolute offset grows linearly with
   wavenumber and reaches about half a grid spacing at the top of the axis, so
   it is at most a sub-bin effect, but it is a systematic one.

Resolution comes from the correlation length, not from padding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The correlation extends to a maximum lag time :math:`T = W\,\Delta t`.
Truncating a correlation at :math:`T` convolves the spectrum with a kernel of
width

.. math::

   \delta\tilde\nu \approx \frac{1}{c\,T}
   \approx \frac{33.4}{T[\mathrm{ps}]}\ \mathrm{cm}^{-1},

and features narrower than that are not resolved no matter how fine the grid
is. Two lines 20 cm⁻¹ apart merge into a single peak for :math:`T = 0.5` ps
(:math:`\delta\tilde\nu = 67` cm⁻¹) and separate cleanly for :math:`T = 2.5` ps
(:math:`\delta\tilde\nu = 13` cm⁻¹). Increasing ``ftsize`` beyond the number of
correlation points only interpolates the same information onto a denser grid;
resolving finer structure requires a longer ``window``, which in turn requires
a longer trajectory.

.. warning::

   ``ftsize`` also truncates. The correlation is zero-padded **or cut** to
   ``ftsize`` points, so with ``ftsize`` smaller than ``window + 1`` everything
   beyond the first ``ftsize`` lags is silently discarded before the transform.
   The default ``ftsize`` is 2000 while the default ``window`` is 1000, so the
   defaults are safe — but raising ``window`` without raising ``ftsize`` throws
   the extra correlation away. Keep ``ftsize`` at least ``window + 1``.

Apodization: leakage against band shape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A correlation that has not decayed to zero at the end of the window produces
sinc-like ringing and leakage in the transform. An apodization window
suppresses the discontinuity, at the cost of widening every band and lowering
every peak [Harris1978]_. Always report which window and which parameters were
used; band widths from differently apodized spectra are not comparable.

The ``exponential`` window multiplies the correlation by
:math:`\exp(-a\,t)` with :math:`a` = ``window_param`` in ps⁻¹. It simply adds
:math:`a` to the decay rate of the correlation, which broadens a Lorentzian
band by roughly

.. math::

   \Delta\tilde\nu_{\mathrm{apod}} \approx \frac{a}{\pi c}
   \approx 10.6\,a[\mathrm{ps}^{-1}]\ \mathrm{cm}^{-1}.

On a test band with a 2 ps dephasing time, the unapodized width is about
9 cm⁻¹, and ``window_param`` values of 2, 4 and 8 ps⁻¹ widen it to about 26, 47
and 89 cm⁻¹ while the peak height drops to 28 %, 16 % and 8 % of its
unapodized value. Apodization strong enough to tame leakage is also strong
enough to dominate the linewidth, so a width read off an apodized spectrum is a
property of the window, not of the dynamics.

.. warning::

   ``hann`` and ``blackman`` do nothing at their default settings. Both are
   built from ``window_start`` (default 0.0 ps) and ``window_stop``
   (default 1000.0 ps), and over a correlation of a few ps the resulting
   factors deviate from unity by less than :math:`3\times10^{-5}` — the
   spectrum is indistinguishable from ``window_function = none``. To use them,
   set ``window_stop`` to the correlation length ``window * time_step``. With
   ``window_stop = 2.5`` on a 2.5 ps correlation, ``hann`` widened the test
   band from 9.2 to 15.0 cm⁻¹ and ``blackman`` to 17.5 cm⁻¹, with peak heights
   falling by 38 % and 45 %. The ``exponential`` window is unaffected by this,
   since it decays from ``window_start`` onwards.

   The legacy ``hann`` and ``blackman`` formulas are additionally non-standard
   — the Hann window is mirrored and the Blackman denominators use the stop
   index rather than the window width. Both quirks are reproduced deliberately;
   see :func:`PQAnalysis.analysis.vacf.spectrum.apodization_window`.

The normalization discards absolute magnitude
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The correlation is normalized to :math:`C(0)=1` before anything else happens:
the ``direct`` estimator divides each time origin by its own aggregate squared
velocity norm, and the ``fft`` estimator divides by its lag-zero value. The
cosine transform is linear, so the spectrum inherits that normalization and its
amplitudes are in arbitrary units. Relative peak *areas* within one spectrum are
meaningful — apodization broadens bands but conserves their area — while peak
*heights* are only comparable between spectra that used the same apodization.
Zero-padding does not rescale amplitudes; it samples the same envelope more
finely, so a coarse grid can under-read the apex of a narrow band by a few per
cent. Absolute intensities are not available at all, for the velocity spectrum
as much as for the charge-flux spectrum.

The sums over atoms are also unweighted — no masses enter anywhere in the VACF
code. Each atom therefore contributes in proportion to its own mean squared
velocity, which at equipartition is proportional to :math:`1/m`, so light atoms
dominate a mixed-element selection. This is not the mass-weighted vibrational
density of states of the textbook definition. For element-resolved band
assignments, run separate analyses with a per-element ``target_selection``.

Classical dynamics limits what a peak position means
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The spectrum describes the classical nuclear motion actually present in the
trajectory. PQAnalysis applies no quantum correction of any kind: no zero-point
energy, no harmonic quantum correction factor to the intensities, no frequency
scaling factor. Band positions therefore carry the classical-nuclei error and
the full anharmonic and thermal shift of the underlying dynamics at the
simulated temperature, and are not directly comparable to harmonic normal-mode
wavenumbers — see :doc:`vibrations` for the harmonic route. Comparisons with
experiment must state the temperature, the level of theory and the fact that
the peak positions are classical.

Estimator choice changes the statistics of the tail
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``method = direct`` (default) spawns an origin every ``gap`` frames while a
  full window still fits. Every origin covers the full window, so every lag
  from 0 to ``window`` is averaged over the same number of origins.
* ``method = fft`` uses every frame as an origin and ignores ``gap``, but
  divides lag :math:`\tau` by its own origin count :math:`N-\tau`. Its origin
  count *does* fall off with lag, so the far tail of the correlation is
  progressively noisier — precisely the part of the correlation that determines
  the low-wavenumber structure of the spectrum. It also holds all velocities in
  memory.

.. warning::

   A velocity trajectory of exactly ``window`` frames with ``gap = 1`` takes
   the legacy single-origin branch: one origin is spawned, the correlation is
   an unaveraged single-origin estimate, and the final lag bin stays exactly
   zero. Use a longer trajectory or a smaller ``window``.

Output and API
--------------

See :ref:`analysis-output-vacf` for correlation and spectrum columns. The
input-file entry point is :func:`PQAnalysis.analysis.vacf.api.vacf`. Direct
calculations use :class:`PQAnalysis.analysis.vacf.vacf.VACF`, while
:func:`PQAnalysis.analysis.vacf.spectrum.vacf_spectrum` performs the spectral
transform.

Discrete line spectra can be broadened independently with
``pqanalysis build_spectrum``; see :ref:`analysis-output-spectrum` for its
output convention.

References
----------

* [Rahman1964]_ introduced the velocity autocorrelation function as a
  molecular-dynamics observable and describes its negative-lobe behavior in a
  liquid.
* [Green1954]_ and [Kubo1957]_ establish the link between equilibrium time
  correlation functions and transport coefficients. The time integral of the
  unnormalized velocity autocorrelation function is the Green-Kubo expression
  for the self-diffusion coefficient.
* [Wiener1930]_ and [Khintchine1934]_ prove that the autocorrelation function
  and the power spectrum of a stationary process are a Fourier pair, which is
  what the ``fft`` estimator exploits.
* [Dickey1969]_ interprets the transformed velocity autocorrelation function
  as a vibrational density of states.
* [Thomas2013]_ covers vibrational density of states, charge-flux and dipole
  routes to infrared spectra, and the effect of correlation depth and
  windowing on band shapes.
* [Harris1978]_ tabulates the apodization windows and their trade-off between
  sidelobe suppression and main-lobe broadening.
* [Allen2017]_ gives the time-origin averaging and sampling requirements for
  correlation functions.
* [thhTools]_ is the legacy program family whose estimator and Fourier
  conventions the default path reproduces.

Full entries are listed in :doc:`../references`.
