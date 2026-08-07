VACF and Spectra
================

The normalized velocity autocorrelation function describes how rapidly atomic
velocities lose memory of their initial direction:

.. math::

   C_{vv}(t) =
   \left\langle
   \frac{\sum_i \mathbf{v}_i(t_0)\cdot\mathbf{v}_i(t_0+t)}
        {\sum_i \mathbf{v}_i(t_0)\cdot\mathbf{v}_i(t_0)}
   \right\rangle_{t_0}.

The brackets denote an average over admissible time origins, and the sum runs
over the selected atoms. This is the default, legacy-compatible estimator and
gives :math:`C_{vv}(0)=1`. The ``fft`` estimator instead averages the numerator
and denominator separately over all available origins before normalization.

PQAnalysis can transform the correlation to a wavenumber-domain spectrum. If
static or time-dependent partial charges are supplied, it correlates
:math:`q_i\mathbf{v}_i` instead, producing a charge-flux spectrum that
approximates an infrared spectrum.

The correlation written to ``out_file`` is not apodized. When a spectrum is
requested, ``window_function`` multiplies a copy of the correlation before the
cosine transform. The optional ``windowed_out_file`` records that copy.

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
* The sampling interval sets the Nyquist limit, and the correlation length sets
  the resolving power. Zero-padding provides a denser frequency grid but does
  not add spectral resolution.
* Apodization reduces endpoint artifacts but changes band widths and
  amplitudes; report the selected function and its parameters.
* Charge-flux spectra require physically meaningful partial charges and should
  not be interpreted as absolute IR intensities without further calibration.

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
