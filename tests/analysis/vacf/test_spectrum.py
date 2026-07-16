"""
Tests for the spectrum functions of the VACF analysis.
"""

import numpy as np
import pytest

from beartype.roar import BeartypeCallHintParamViolation

from PQAnalysis.analysis.vacf.exceptions import VACFError
from PQAnalysis.analysis.vacf.spectrum import (
    BLACKMAN_A0,
    BLACKMAN_A1,
    BLACKMAN_A2,
    SPEED_OF_LIGHT_CM_S,
    WINDOW_FUNCTIONS,
    _legacy_nint,
    apodization_window,
    vacf_spectrum,
)

from .. import pytestmark  # pylint: disable=unused-import

# pylint: disable=protected-access



def _direct_cosine_spectrum(correlation, ftsize):
    """
    Literal (slow) transcription of the legacy ft.f cosine sum.
    """
    padded = np.zeros(ftsize)
    n_points = min(len(correlation), ftsize)
    padded[:n_points] = correlation[:n_points]

    # legacy even extension cd(1..2*ftsize-1) centered at ftsize
    cd = np.zeros(2 * ftsize - 1)
    cd[ftsize - 1] = padded[0]
    for i in range(2, ftsize + 1):
        cd[ftsize - (i - 1) - 1] = padded[i - 1]
        cd[ftsize + (i - 1) - 1] = padded[i - 1]

    fu = 2.0 * np.pi / (2 * ftsize - 1)

    amplitudes = np.zeros(ftsize)
    for spectral_index in range(1, ftsize + 1):
        total = 0.0
        for k in range(1, 2 * ftsize):
            total += cd[k - 1] * np.cos(
                fu * spectral_index * (k - ftsize)
            )
        amplitudes[spectral_index - 1] = abs(0.5 * total)

    return amplitudes



class TestLegacyNint:

    """
    Tests for the legacy FORTRAN nint emulation.
    """

    def test_rounds_half_away_from_zero(self):
        """
        Half-way cases are rounded away from zero.
        """
        assert _legacy_nint(2.5) == 3
        assert _legacy_nint(-2.5) == -3
        assert _legacy_nint(0.5) == 1
        assert _legacy_nint(2.4) == 2
        assert _legacy_nint(-2.4) == -2
        assert _legacy_nint(3.0) == 3



class TestApodizationWindow:

    """
    Tests for the legacy apodization windows.
    """

    def test_supported_window_functions(self):
        """
        The supported window functions are none, exponential, hann
        and blackman.
        """
        assert WINDOW_FUNCTIONS == ("none", "exponential", "hann", "blackman")

    def test_none_window(self):
        """
        The none window consists of ones only.
        """
        factors = apodization_window(10, 0.1, window_function="none")

        assert np.all(factors == 1.0)

    def test_exponential_window_shape(self):
        """
        The exponential window is one before the start index, decays
        with exp(-a dt (i - 1 - winsind)) inside the window range and
        is zero after the end index.
        """
        factors = apodization_window(
            20,
            0.1,
            window_function="exponential",
            window_param=2.0,
            window_start=0.5,
            window_stop=1.5,
        )

        # winsind = nint(0.5 / 0.1) = 5, wineind = nint(1.5 / 0.1) = 15
        index = np.arange(1, 21)
        expected = np.exp(-2.0 * 0.1 * (index - 1 - 5))
        expected[index <= 5] = 1.0
        expected[index > 15] = 0.0

        assert np.allclose(factors, expected, atol=1e-14)

    def test_hann_window_shape(self):
        """
        The legacy hann window is mirrored: it rises from the end of
        the window range and reaches zero exactly at the end index.
        """
        factors = apodization_window(
            20,
            0.1,
            window_function="hann",
            window_start=0.5,
            window_stop=1.5,
        )

        index = np.arange(1, 21)
        expected = (1.0 - np.cos(np.pi * (15 - index) / 10)) / 2.0
        expected[index <= 5] = 1.0
        expected[index > 15] = 0.0

        assert np.allclose(factors, expected, atol=1e-14)
        assert factors[14] == 0.0  # i = wineind = 15

    def test_blackman_window_shape(self):
        """
        The legacy blackman window uses the (non-standard) end index
        as denominator and single precision leading coefficients.
        """
        factors = apodization_window(
            20,
            0.1,
            window_function="blackman",
            window_start=0.5,
            window_stop=1.5,
        )

        index = np.arange(1, 21)
        phase = np.pi * (index - 1 - 5) / 15
        expected = (
            BLACKMAN_A0 + BLACKMAN_A1 * np.cos(phase) +
            BLACKMAN_A2 * np.cos(2.0 * phase)
        )
        expected[index <= 5] = 1.0
        expected[index > 15] = 0.0

        assert BLACKMAN_A0 == float(np.float32(0.42))
        assert BLACKMAN_A2 == float(np.float32(0.08))
        assert np.allclose(factors, expected, atol=1e-14)

    def test_unknown_window_function(self):
        """
        An unknown window function raises a VACFError.
        """
        with pytest.raises(VACFError, match="Unknown window function"):
            apodization_window(10, 0.1, window_function="hamming")

    def test_empty_window_range(self):
        """
        An empty or inverted window range raises a VACFError.
        """
        with pytest.raises(VACFError, match="window range is empty"):
            apodization_window(
                10,
                0.1,
                window_function="hann",
                window_start=1.5,
                window_stop=0.5,
            )



class TestVACFSpectrum:

    """
    Tests for the legacy cosine-transform spectrum.
    """

    def test_matches_direct_cosine_sum(self):
        """
        The rfft based implementation is numerically identical to the
        literal legacy cosine sum.
        """
        rng = np.random.default_rng(42)
        correlation = rng.standard_normal(40)
        time = np.arange(40) * 0.1

        _, amplitudes, _ = vacf_spectrum(time, correlation, ftsize=64)

        reference = _direct_cosine_spectrum(correlation, 64)

        assert np.allclose(amplitudes, reference, atol=1e-9)

    def test_matches_direct_cosine_sum_with_truncation(self):
        """
        A correlation function longer than ftsize is truncated like in
        the legacy tool.
        """
        rng = np.random.default_rng(43)
        correlation = rng.standard_normal(100)
        time = np.arange(100) * 0.05

        _, amplitudes, _ = vacf_spectrum(time, correlation, ftsize=32)

        reference = _direct_cosine_spectrum(correlation, 32)

        assert np.allclose(amplitudes, reference, atol=1e-9)

    def test_frequency_axis_replicates_legacy_calibration(self):
        """
        The frequency axis replicates the historical ft.f calibration
        with a period of 2 * (ftsize - 1) points although the even
        extension has 2 * ftsize - 1 points.
        """
        time = np.arange(10) * 0.002
        correlation = np.ones(10)

        wavenumbers, amplitudes, _ = vacf_spectrum(
            time,
            correlation,
            ftsize=16,
        )

        spacing = 1.0 / (0.002e-12 * 2.0 * 15 * SPEED_OF_LIGHT_CM_S)

        assert len(wavenumbers) == 16
        assert np.allclose(
            wavenumbers,
            np.arange(1, 17) * spacing,
            rtol=1e-12,
        )
        # the last point duplicates its predecessor (legacy quirk of
        # evaluating the cosine sum at l = ftsize)
        assert amplitudes[-1] == amplitudes[-2]

    def test_single_cosine_peak_position(self):
        """
        The spectrum of a pure cosine peaks at the corresponding
        wavenumber (within one frequency spacing).
        """
        time_step = 0.002
        frequency = 25.0  # ps^-1
        time = np.arange(101) * time_step
        correlation = np.cos(2.0 * np.pi * frequency * time)

        wavenumbers, amplitudes, _ = vacf_spectrum(
            time,
            correlation,
            ftsize=256,
        )

        expected = frequency * 1.0e12 / SPEED_OF_LIGHT_CM_S
        spacing = wavenumbers[0]

        assert abs(wavenumbers[np.argmax(amplitudes)] - expected) <= spacing

    def test_windowed_correlation_output(self):
        """
        The windowed correlation function equals the input multiplied
        with the apodization window.
        """
        rng = np.random.default_rng(44)
        correlation = rng.standard_normal(30)
        time = np.arange(30) * 0.1

        _, _, windowed = vacf_spectrum(
            time,
            correlation,
            ftsize=32,
            window_function="hann",
            window_start=0.5,
            window_stop=2.0,
        )

        factors = apodization_window(
            30,
            0.1,
            window_function="hann",
            window_start=0.5,
            window_stop=2.0,
        )

        assert np.allclose(windowed, correlation * factors, atol=1e-14)

    def test_length_mismatch(self):
        """
        A length mismatch between the time axis and the correlation
        function raises a VACFError.
        """
        with pytest.raises(VACFError, match="same length"):
            vacf_spectrum(np.arange(5) * 0.1, np.ones(4))

    def test_too_few_points(self):
        """
        Less than two points raise a VACFError.
        """
        with pytest.raises(VACFError, match="At least two"):
            vacf_spectrum(np.zeros(1), np.ones(1))

    def test_too_small_ftsize(self):
        """
        An ftsize smaller than two raises a VACFError.
        """
        with pytest.raises(VACFError, match="ftsize"):
            vacf_spectrum(np.arange(5) * 0.1, np.ones(5), ftsize=1)

    def test_non_increasing_time_axis(self):
        """
        A non-increasing time axis raises a VACFError.
        """
        with pytest.raises(VACFError, match="strictly increasing"):
            vacf_spectrum(np.zeros(5), np.ones(5))

    def test_non_one_dimensional_input(self):
        """
        Multi-dimensional input raises a VACFError (or a beartype
        violation when running with beartype in DEBUG mode).
        """
        with pytest.raises((VACFError, BeartypeCallHintParamViolation)):
            vacf_spectrum(np.zeros((5, 2)), np.ones((5, 2)))
