VACF and Spectra
================

The normalized velocity autocorrelation function describes how rapidly atomic
velocities lose memory of their initial direction:

.. math::

   C_v(t) =
   \frac{\left\langle \sum_i \mathbf{v}_i(0)\cdot\mathbf{v}_i(t)\right\rangle}
        {\left\langle \sum_i \mathbf{v}_i(0)\cdot\mathbf{v}_i(0)\right\rangle}.

PQAnalysis can transform the correlation to a wavenumber-domain spectrum. If
static or time-dependent partial charges are supplied, it correlates
:math:`q_i\mathbf{v}_i` instead, producing a charge-flux spectrum that
approximates an infrared spectrum.

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

The time step is specified in ps. ``window_function`` accepts
``exponential``, ``hann`` and ``blackman``. The default sliding-origin method
matches the legacy calculation; ``method = fft`` selects a denser-origin
Wiener-Khinchin estimator.

Interpretation
--------------

* A rapidly decaying VACF indicates fast velocity decorrelation.
* Negative regions indicate backscattering or cage motion.
* The frequency spectrum depends on the sampling interval, correlation length,
  apodization window and zero-padding size.
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
