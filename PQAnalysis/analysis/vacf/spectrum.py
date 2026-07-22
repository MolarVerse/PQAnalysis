"""
A module containing the spectrum functions of the VACF analysis.

The functions in this module port the legacy ``ft.f`` Fourier
transformation tool (``ftvac`` ver. 2.0) of the ``thh_tools``
collection. The auto-correlation function is optionally multiplied
with a one-sided apodization window, mirrored into an even extension
and transformed with a discrete cosine sum. All legacy conventions are
replicated exactly - see the notes of :py:func:`vacf_spectrum` for the
historical quirks that are kept for backwards compatibility.
"""

import numpy as np

from beartype.typing import Tuple

from PQAnalysis.types import Np1DNumberArray, PositiveInt, PositiveReal

from .exceptions import VACFError

SPEED_OF_LIGHT_CM_S = 2.99792458e10

#: The window functions supported by :py:func:`apodization_window`.
WINDOW_FUNCTIONS = ("none", "exponential", "hann", "blackman")

# NOTE: the legacy ft.f writes the Blackman coefficients as default-real
#       (single precision) FORTRAN literals. The constants are rounded
#       to single precision here to replicate the legacy binary bit for
#       bit; the plain double precision values would differ by ~1e-8.
BLACKMAN_A0 = float(np.float32(0.42))
BLACKMAN_A1 = 0.5
BLACKMAN_A2 = float(np.float32(0.08))



def _legacy_nint(value: float) -> int:
    """
    Rounds a value to the nearest integer like the FORTRAN ``nint``.

    Half-way cases are rounded away from zero instead of the banker's
    rounding of the built-in ``round``.

    Parameters
    ----------
    value : float
        The value to round.

    Returns
    -------
    int
        The rounded value.
    """
    if value >= 0.0:
        return int(np.floor(value + 0.5))

    return int(np.ceil(value - 0.5))



def apodization_window(
    n_points: PositiveInt,
    time_step: PositiveReal,
    window_function: str = "none",
    window_param: PositiveReal = 4.0,
    window_start: PositiveReal = 0.0,
    window_stop: PositiveReal = 1000.0,
) -> Np1DNumberArray:
    """
    Calculates the one-sided legacy apodization window.

    The window replicates the legacy ``ft.f`` window functions exactly.
    With the start index ``winsind = nint(window_start / time_step)``
    and the end index ``wineind = nint(window_stop / time_step)`` the
    window factor of the one-based point ``i`` is one for
    ``i <= winsind``, zero for ``i > wineind`` and otherwise:

    - ``exponential``: ``exp(-window_param * time_step * (i - 1 - winsind))``
    - ``hann``: ``(1 - cos(pi * (wineind - i) / (wineind - winsind))) / 2``
    - ``blackman``: ``0.42 + 0.5 * cos(pi * (i - 1 - winsind) / wineind)
      + 0.08 * cos(2 * pi * (i - 1 - winsind) / wineind)``

    Note that the legacy ``hann`` and ``blackman`` formulas are
    non-standard: the ``hann`` window is mirrored (it rises from the
    end of the window range) and the ``blackman`` denominators are
    ``wineind`` instead of the window width. Both quirks are replicated
    deliberately.

    Parameters
    ----------
    n_points : PositiveInt
        The number of points of the correlation function.
    time_step : PositiveReal
        The time step between two points of the correlation function.
    window_function : str, optional
        The window function, one of ``none``, ``exponential``, ``hann``
        and ``blackman``, by default ``none`` (all factors are one).
    window_param : PositiveReal, optional
        The exponential decay coefficient ``a`` of the ``exponential``
        window ``exp(-a * t)``, by default 4.0.
    window_start : PositiveReal, optional
        The time at which the window starts to decay, by default 0.0.
    window_stop : PositiveReal, optional
        The time at which the window becomes zero, by default 1000.0.

    Returns
    -------
    Np1DNumberArray
        The window factors for all points.

    Raises
    ------
    VACFError
        If the window function is unknown.
    VACFError
        If the window range is empty or inverted.
    """
    window_function = window_function.lower()

    if window_function not in WINDOW_FUNCTIONS:
        raise VACFError(
            f"Unknown window function '{window_function}'. Possible "
            f"window functions are: {', '.join(WINDOW_FUNCTIONS)}."
        )

    if window_function == "none":
        return np.ones(n_points, dtype=np.float64)

    win_start_index = _legacy_nint(window_start / time_step)
    win_stop_index = _legacy_nint(window_stop / time_step)

    if win_stop_index <= win_start_index or win_stop_index <= 0:
        raise VACFError(
            "The window range is empty: window_stop must resolve to a "
            "larger index than window_start and must be positive."
        )

    index = np.arange(1, n_points + 1, dtype=np.float64)

    if window_function == "exponential":
        factors = np.exp(
            -window_param * time_step * (index - 1 - win_start_index)
        )
    elif window_function == "hann":
        factors = (
            1.0 - np.cos(
                np.pi * (win_stop_index - index) /
                (win_stop_index - win_start_index)
            )
        ) / 2.0
    else:  # blackman
        phase = np.pi * (index - 1 - win_start_index) / win_stop_index
        factors = (
            BLACKMAN_A0 + BLACKMAN_A1 * np.cos(phase) +
            BLACKMAN_A2 * np.cos(2.0 * phase)
        )

    factors = np.where(index <= win_start_index, 1.0, factors)
    factors = np.where(index > win_stop_index, 0.0, factors)

    return factors



def vacf_spectrum(
    time: Np1DNumberArray,
    correlation: Np1DNumberArray,
    ftsize: PositiveInt = 2000,
    window_function: str = "none",
    window_param: PositiveReal = 4.0,
    window_start: PositiveReal = 0.0,
    window_stop: PositiveReal = 1000.0,
) -> Tuple[Np1DNumberArray, Np1DNumberArray, Np1DNumberArray]:
    """
    Calculates the legacy cosine-transform spectrum of a correlation
    function.

    The (optionally apodized) correlation function is zero-padded to
    ``ftsize`` points and mirrored into an even extension of length
    ``2 * ftsize - 1`` centered at index ``ftsize``. The spectrum is
    the discrete cosine sum

    ``spectrum(l) = |0.5 * sum_k cd(k) *
    cos(2 * pi * l * (k - ftsize) / (2 * ftsize - 1))|``

    for ``l = 1..ftsize``, which is evaluated here via the equivalent
    real part of the FFT of the circularly shifted even extension.

    Notes
    -----
    The frequency axis replicates the historical calibration of the
    legacy ``ft.f`` tool: the wavenumber spacing is calculated with a
    period of ``2 * (ftsize - 1)`` points although the underlying even
    extension has ``2 * ftsize - 1`` points. This slight frequency-axis
    mismatch is kept deliberately so that spectra remain comparable
    with the historical results. The time step is taken from the first
    two entries of the time axis and is assumed to be constant, and the
    last spectrum point duplicates its predecessor (a legacy
    consequence of evaluating the cosine sum at ``l = ftsize``).

    Parameters
    ----------
    time : Np1DNumberArray
        The equidistant time axis of the correlation function in ps.
    correlation : Np1DNumberArray
        The correlation function values.
    ftsize : PositiveInt, optional
        The Fourier transform point size, by default 2000. The
        correlation function is zero-padded (or truncated) to this
        size.
    window_function : str, optional
        The apodization window, one of ``none``, ``exponential``,
        ``hann`` and ``blackman``, by default ``none``. See
        :py:func:`apodization_window`.
    window_param : PositiveReal, optional
        The exponential window decay coefficient, by default 4.0.
    window_start : PositiveReal, optional
        The window start time, by default 0.0.
    window_stop : PositiveReal, optional
        The window stop time, by default 1000.0.

    Returns
    -------
    wavenumbers : Np1DNumberArray
        The wavenumbers in cm^-1 for the indices ``1..ftsize``.
    amplitudes : Np1DNumberArray
        The spectrum amplitudes.
    windowed_correlation : Np1DNumberArray
        The correlation function after applying the apodization window
        (equal to the input for the ``none`` window).

    Raises
    ------
    VACFError
        If less than two points are given, the time axis and the
        correlation function have different lengths or ftsize is
        smaller than two.
    """
    time = np.asarray(time, dtype=np.float64)
    correlation = np.asarray(correlation, dtype=np.float64)

    if time.ndim != 1 or correlation.ndim != 1:
        raise VACFError(
            "The time axis and the correlation function have to be "
            "one-dimensional arrays."
        )

    if len(time) != len(correlation):
        raise VACFError(
            "The time axis and the correlation function must have the "
            "same length."
        )

    if len(time) < 2:
        raise VACFError(
            "At least two correlation points are needed to infer the "
            "time step."
        )

    if ftsize < 2:
        raise VACFError("The ftsize must be at least 2.")

    time_step = time[1] - time[0]

    if time_step <= 0.0:
        raise VACFError("The time axis must be strictly increasing.")

    factors = apodization_window(
        len(correlation),
        time_step,
        window_function=window_function,
        window_param=window_param,
        window_start=window_start,
        window_stop=window_stop,
    )

    windowed_correlation = correlation * factors

    padded = np.zeros(ftsize, dtype=np.float64)
    n_points = min(len(windowed_correlation), ftsize)
    padded[:n_points] = windowed_correlation[:n_points]

    # even extension [c_1 .. c_ftsize, c_ftsize .. c_2] of length
    # 2 * ftsize - 1; its DFT is real and equals the legacy cosine sum
    extension = np.concatenate([padded, padded[-1:0:-1]])
    transform = np.fft.rfft(extension).real

    amplitudes = np.abs(0.5 * transform)

    # the legacy loop runs to l = ftsize, which mirrors back onto the
    # rfft bin ftsize - 1 of the odd-length even extension
    amplitudes = np.concatenate([amplitudes[1:], amplitudes[-1:]])

    frequency_spacing = 1.0 / (
        time_step * 1.0e-12 * 2.0 * (ftsize - 1) * SPEED_OF_LIGHT_CM_S
    )
    wavenumbers = np.arange(1, ftsize + 1, dtype=np.float64)
    wavenumbers = wavenumbers * frequency_spacing

    return wavenumbers, amplitudes, windowed_correlation
