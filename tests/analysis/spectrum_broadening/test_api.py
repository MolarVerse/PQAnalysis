"""
Tests for the spectrum broadening API and output writer.
"""

import numpy as np
import pytest

from PQAnalysis.analysis.spectrum_broadening import (
    SpectrumDataWriter,
    build_spectrum,
)
from PQAnalysis.analysis.spectrum_broadening.exceptions import (
    SpectrumBroadeningError,
)

from .. import pytestmark  # pylint: disable=unused-import



class TestBuildSpectrumAPI:

    """
    Tests for the build_spectrum API function.
    """

    @pytest.mark.parametrize(
        "example_dir", ["spectrum_broadening"], indirect=False
    )
    def test_legacy_parity(self, test_with_data_dir):  # pylint: disable=unused-argument
        """
        The broadened spectrum of the legacy 78-line stick spectrum
        matches the output of the legacy build_spectrum.sh gawk
        implementation on the default grid.

        The reference file spectrum.dat was produced by the legacy
        awk loop with BROAD=0.0025 from lines.dat
        (paste frequencies.dat intensity.dat of the legacy test case).
        """
        grid, values = build_spectrum("lines.dat", output="out.dat")

        reference = np.loadtxt("spectrum.dat")

        assert grid.size == 15960
        assert np.array_equal(grid, reference[:, 0])
        assert np.allclose(
            values, reference[:, 1], rtol=1e-10, atol=1e-300
        )

        with open("out.dat", encoding="utf-8") as file:
            lines = file.readlines()

        with open("spectrum.dat", encoding="utf-8") as file:
            reference_lines = file.readlines()

        assert len(lines) == len(reference_lines)

        # the vast majority of the lines must be byte-identical to the
        # legacy output; the remaining ones only differ in the last
        # digit of the mantissa or in the subnormal underflow region
        identical = sum(
            1 for line, reference_line in zip(lines, reference_lines)
            if line == reference_line
        )
        assert identical >= 15900

    def test_stdout_output(self, tmpdir, capsys):  # pylint: disable=unused-argument
        """
        Without an output file the spectrum is printed to stdout in
        the legacy '%8.4f    %16.12e' format.
        """
        with open("sticks.dat", "w", encoding="utf-8") as file:
            file.write("15.0 2.0\n")

        build_spectrum(
            "sticks.dat",
            wavenumber_min=10.0,
            wavenumber_max=11.0,
            wavenumber_step=0.5,
        )

        captured = capsys.readouterr()
        expected_first = 2.0 * np.exp(-0.0025 * (10.0 - 15.0)**2)
        expected_second = 2.0 * np.exp(-0.0025 * (10.5 - 15.0)**2)

        assert captured.out == (
            f" 10.0000    {expected_first:16.12e}\n"
            f" 10.5000    {expected_second:16.12e}\n"
        )

    def test_fwhm_option(self, tmpdir):  # pylint: disable=unused-argument
        """
        The fwhm parameter is an alternative way to set the width.
        """
        with open("sticks.dat", "w", encoding="utf-8") as file:
            file.write("100.0 2.0\n")

        _, values = build_spectrum(
            "sticks.dat",
            output="out.dat",
            fwhm=20.0,
            wavenumber_min=100.0,
            wavenumber_max=120.0,
            wavenumber_step=10.0,
        )

        assert values[0] == pytest.approx(2.0, rel=1e-15)
        assert values[1] == pytest.approx(1.0, rel=1e-12)

    def test_alpha_and_fwhm_are_mutually_exclusive(self, tmpdir):  # pylint: disable=unused-argument
        """
        Specifying both alpha and fwhm is rejected.
        """
        with open("sticks.dat", "w", encoding="utf-8") as file:
            file.write("100.0 2.0\n")

        with pytest.raises(SpectrumBroadeningError):
            build_spectrum("sticks.dat", alpha=0.0025, fwhm=20.0)



class TestSpectrumDataWriter:

    """
    Tests for the SpectrumDataWriter class.
    """

    def test_write(self, tmpdir):  # pylint: disable=unused-argument
        """
        The writer produces the legacy '%8.4f    %16.12e' rows.
        """
        writer = SpectrumDataWriter("spectrum_out.dat")
        writer.write(
            (
                np.array([10.0, 3999.75]),
                np.array([1.918547669503, 0.0]),
            )
        )

        with open("spectrum_out.dat", encoding="utf-8") as file:
            lines = file.readlines()

        assert lines == [
            " 10.0000    1.918547669503e+00\n",
            "3999.7500    0.000000000000e+00\n",
        ]
