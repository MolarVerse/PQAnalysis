"""
Tests for the spectrum broadening numerical routines.
"""

import numpy as np
import pytest

from PQAnalysis.analysis.spectrum_broadening import (
    DEFAULT_ALPHA,
    alpha_from_fwhm,
    broaden,
    fwhm_from_alpha,
    read_stick_spectrum,
    wavenumber_grid,
)
from PQAnalysis.analysis.spectrum_broadening import (
    spectrum_broadening as spectrum_broadening_module,
)
from PQAnalysis.analysis.spectrum_broadening.exceptions import (
    SpectrumBroadeningError,
)

from .. import pytestmark  # pylint: disable=unused-import



class TestWidthConversion:

    """
    Tests for the alpha/FWHM conversion helpers.
    """

    def test_alpha_from_fwhm_roundtrip(self):
        """
        The alpha/FWHM conversions are inverse to each other and the
        default alpha corresponds to a FWHM of about 33.3 cm^-1.
        """
        fwhm = fwhm_from_alpha(DEFAULT_ALPHA)
        assert fwhm == pytest.approx(33.302185, rel=1e-6)
        assert alpha_from_fwhm(fwhm) == pytest.approx(
            DEFAULT_ALPHA, rel=1e-12
        )

    def test_non_positive_values(self):
        """
        Non-positive widths are rejected.
        """
        with pytest.raises(SpectrumBroadeningError):
            alpha_from_fwhm(0.0)

        with pytest.raises(SpectrumBroadeningError):
            fwhm_from_alpha(0.0)



class TestWavenumberGrid:

    """
    Tests for the wavenumber grid construction.
    """

    def test_default_grid_matches_legacy_loop(self):
        """
        The default grid reproduces the legacy awk loop
        'for (i=10; i<4000; i+=0.25)' with 15960 points and an
        exclusive upper bound.
        """
        grid = wavenumber_grid()

        assert grid.size == 15960
        assert grid[0] == 10.0
        assert grid[-1] == 3999.75
        assert grid.dtype == np.float64

    def test_custom_grid(self):
        """
        A custom grid honors min, max and step.
        """
        grid = wavenumber_grid(0.0, 1.0, 0.5)
        assert np.array_equal(grid, [0.0, 0.5])

    def test_invalid_grid_settings(self):
        """
        Invalid grid settings are rejected.
        """
        with pytest.raises(SpectrumBroadeningError):
            wavenumber_grid(10.0, 4000.0, 0.0)

        with pytest.raises(SpectrumBroadeningError):
            wavenumber_grid(4000.0, 10.0, 0.25)



class TestBroaden:

    """
    Tests for the broadening kernel evaluation.
    """

    def test_gaussian_peak_height_convention(self):
        """
        A single stick located on a grid point reaches exactly its
        intensity there (no area normalization).
        """
        grid = wavenumber_grid(1990.0, 2010.0, 0.25)
        result = broaden(
            np.array([2000.0]), np.array([3.0]), grid, alpha=0.0025
        )

        assert result[np.where(grid == 2000.0)[0][0]] == pytest.approx(
            3.0, rel=1e-15
        )

    def test_gaussian_matches_direct_formula(self):
        """
        The vectorized broadening matches a direct per-stick scalar
        evaluation of sum_k I_k * exp(-alpha (g - nu_k)^2).
        """
        wavenumbers = np.array([36.7203, 83.2782, 114.5203])
        intensities = np.array([11.4330, 0.8398, 0.9518])
        grid = wavenumber_grid(10.0, 200.0, 0.25)
        alpha = 0.0025

        expected = np.zeros(grid.size)
        for stick, intensity in zip(wavenumbers, intensities):
            for index, grid_point in enumerate(grid):
                expected[index] += intensity * np.exp(
                    -alpha * (grid_point - stick)**2
                )

        result = broaden(wavenumbers, intensities, grid, alpha=alpha)

        assert np.allclose(result, expected, rtol=1e-13)

    def test_gaussian_fwhm(self):
        """
        The broadened profile of a single stick falls to half of its
        peak height at half the FWHM away from the stick.
        """
        fwhm = 20.0
        alpha = alpha_from_fwhm(fwhm)
        grid = np.array([1000.0, 1000.0 + fwhm / 2.0])

        result = broaden(np.array([1000.0]), np.array([2.0]), grid, alpha)

        assert result[0] == pytest.approx(2.0, rel=1e-15)
        assert result[1] == pytest.approx(1.0, rel=1e-12)

    def test_lorentzian_fwhm(self):
        """
        The Lorentzian kernel shares peak height and FWHM with the
        Gaussian kernel for the same alpha.
        """
        fwhm = 20.0
        alpha = alpha_from_fwhm(fwhm)
        grid = np.array([1000.0, 1000.0 + fwhm / 2.0])

        result = broaden(
            np.array([1000.0]),
            np.array([2.0]),
            grid,
            alpha,
            kernel="lorentzian",
        )

        assert result[0] == pytest.approx(2.0, rel=1e-15)
        assert result[1] == pytest.approx(1.0, rel=1e-12)

    def test_chunked_evaluation(self, monkeypatch):
        """
        Chunked evaluation gives the same result as a single chunk.
        """
        rng = np.random.default_rng(42)
        wavenumbers = rng.uniform(10.0, 4000.0, 17)
        intensities = rng.uniform(0.0, 10.0, 17)
        grid = wavenumber_grid(10.0, 500.0, 0.25)

        expected = broaden(wavenumbers, intensities, grid)

        monkeypatch.setattr(
            spectrum_broadening_module, "_GRID_CHUNK_SIZE", 7
        )
        result = broaden(wavenumbers, intensities, grid)

        # bit-exact equality is not guaranteed because BLAS may sum
        # small matrix products in a different order
        assert np.allclose(result, expected, rtol=1e-12, atol=0.0)

    def test_empty_sticks(self):
        """
        An empty stick spectrum broadens to all zeros like the legacy
        awk implementation.
        """
        grid = wavenumber_grid(10.0, 20.0, 1.0)
        result = broaden(np.array([]), np.array([]), grid)

        assert np.array_equal(result, np.zeros(grid.size))

    def test_invalid_inputs(self):
        """
        Invalid alpha, kernel and shape settings are rejected.
        """
        grid = wavenumber_grid(10.0, 20.0, 1.0)

        with pytest.raises(SpectrumBroadeningError):
            broaden(np.array([1.0]), np.array([1.0]), grid, alpha=0.0)

        with pytest.raises(SpectrumBroadeningError):
            broaden(np.array([1.0]), np.array([1.0]), grid, kernel="voigt")

        with pytest.raises(SpectrumBroadeningError):
            broaden(np.array([1.0, 2.0]), np.array([1.0]), grid)



class TestReadStickSpectrum:

    """
    Tests for the stick spectrum file reader.
    """

    def test_read(self, tmpdir):  # pylint: disable=unused-argument
        """
        A two-column file with comments and blank lines is parsed
        into float64 arrays.
        """
        with open("sticks.dat", "w", encoding="utf-8") as file:
            file.write("# comment\n")
            file.write("36.7203    11.4330\n")
            file.write("\n")
            file.write("83.2782    0.8398 extra-column\n")

        wavenumbers, intensities = read_stick_spectrum("sticks.dat")

        assert np.allclose(wavenumbers, [36.7203, 83.2782])
        assert np.allclose(intensities, [11.4330, 0.8398])
        assert wavenumbers.dtype == np.float64
        assert intensities.dtype == np.float64

    def test_missing_file(self):
        """
        A missing file is rejected.
        """
        with pytest.raises(SpectrumBroadeningError):
            read_stick_spectrum("does-not-exist.dat")

    def test_malformed_lines(self, tmpdir):  # pylint: disable=unused-argument
        """
        Files with too few columns or non-numeric data are rejected.
        """
        with open("one_column.dat", "w", encoding="utf-8") as file:
            file.write("36.7203\n")

        with pytest.raises(SpectrumBroadeningError):
            read_stick_spectrum("one_column.dat")

        with open("non_numeric.dat", "w", encoding="utf-8") as file:
            file.write("36.7203 abc\n")

        with pytest.raises(SpectrumBroadeningError):
            read_stick_spectrum("non_numeric.dat")
