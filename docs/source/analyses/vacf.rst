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

The cosine transform of the correlation is the vibrational density of states
[Dickey1969]_, [Thomas2013]_. If partial charges are supplied, PQAnalysis
correlates :math:`q_i\mathbf{v}_i` instead and the transform approximates an
infrared spectrum [Thomas2013]_.

.. plot:: _plots/vacf.py
   :alt: Analytical normalized VACF, its exponentially windowed copy and the resulting spectrum
   :caption: Analytical normalized VACF for two Gaussian-broadened bands
      centered at 300 and 600 cm⁻¹, with dephasing times of 0.22 and 0.12 ps.
      The dashed curve applies an exponential window with a decay coefficient
      of 4 ps⁻¹ before the PQAnalysis cosine transform. Spectrum amplitudes
      are scaled to unit maximum.

Input
-----

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

Saved as ``vacf.in``, it runs with:

.. code-block:: console

   $ pqanalysis vacf vacf.in

``time_step`` is in ps. ``window`` is the maximum correlation lag in frames
(the bundled :doc:`../examples` fixture uses ``window = 8``); ``gap`` spaces
the time origins. ``method = fft`` selects the denser-origin Wiener-Khinchin
estimator.

The correlation written to ``out_file`` is never apodized. When a spectrum is
requested, ``window_function`` (``exponential``, ``hann``, ``blackman`` or the
default ``none``) multiplies a copy before the cosine transform
[Harris1978]_; the example applies :math:`\exp[-(4\ \mathrm{ps}^{-1})t]`. The
optional ``windowed_out_file`` records that copy. The full key table is on
:class:`~PQAnalysis.analysis.vacf.vacf_input_file_reader.VACFInputFileReader`.

Output
------

``out_file`` holds lag time and normalized correlation; ``spectrum_file`` holds
wavenumber and relative amplitude (:ref:`analysis-output-vacf`). In liquids,
negative regions of the VACF indicate backscattering or cage motion; in solids,
sign oscillations reflect bound vibrational motion. Charge-flux spectra are only
as good as the partial charges behind them and are not absolute IR intensities.
Discrete line spectra can be broadened separately with
``pqanalysis build_spectrum`` (:ref:`analysis-output-spectrum`).

Python
------

.. code-block:: python

   from PQAnalysis.analysis import VACF
   from PQAnalysis.analysis.vacf import vacf_spectrum
   from PQAnalysis.io import TrajectoryReader

   time, correlation = VACF(
       TrajectoryReader("examples/water/trajectory.vel"),
       time_step=0.001,
       window_size=8,
       gap=2,
   ).run()
   wavenumbers, amplitudes, windowed = vacf_spectrum(
       time,
       correlation,
       ftsize=5000,
       window_function="exponential",
       window_param=4.0,
   )

:func:`~PQAnalysis.analysis.vacf.api.vacf` runs an input file instead.

Before you trust it
-------------------

* The Nyquist limit is set by the interval between *written* frames:
  16678 cm⁻¹ at 1 fs, 1668 cm⁻¹ at 10 fs. Faster motion is aliased, not
  dropped.
* Resolution comes from ``window * time_step``, about
  :math:`33.4 / T[\mathrm{ps}]` cm⁻¹. ``ftsize`` only interpolates, and it
  truncates the correlation if it is smaller than ``window + 1``.
* Apodization widens every band and lowers every peak; report the window and
  its parameter. ``hann`` and ``blackman`` do nothing until ``window_stop`` is
  set to the correlation length.
* Amplitudes are arbitrary units; the sum is not mass-weighted; peak positions
  are classical.

Each point is derived in :doc:`vacf-details`.

.. toctree::
   :hidden:

   vacf-details
