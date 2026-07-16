"""
Tests for the build_spectrum CLI.
"""

import sys

from PQAnalysis.cli import build_spectrum as build_spectrum_cli
from PQAnalysis.cli.build_spectrum import BuildSpectrumCLI
from PQAnalysis.io.formats import FileWritingMode



class TestBuildSpectrumCLI:

    """
    Tests for BuildSpectrumCLI.
    """

    def test_program_name(self):
        """
        The program name is build_spectrum.
        """
        assert BuildSpectrumCLI.program_name() == "build_spectrum"

    def test_main_dispatches_defaults(self, monkeypatch):
        """
        main() dispatches the positional input file and the default
        settings to the build_spectrum API function.
        """
        called = []

        monkeypatch.setattr(
            build_spectrum_cli,
            "build_spectrum",
            lambda **kwargs: called.append(kwargs),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            ["build_spectrum", "lines.dat", "--log-file", "off"],
        )

        build_spectrum_cli.main()

        assert called == [
            {
                "input_file": "lines.dat",
                "output": None,
                "alpha": None,
                "fwhm": None,
                "wavenumber_min": 10.0,
                "wavenumber_max": 4000.0,
                "wavenumber_step": 0.25,
                "kernel": "gaussian",
                "mode": FileWritingMode("w"),
            }
        ]

    def test_main_dispatches_options(self, monkeypatch):
        """
        main() dispatches output, width, grid and kernel options.
        """
        called = []

        monkeypatch.setattr(
            build_spectrum_cli,
            "build_spectrum",
            lambda **kwargs: called.append(kwargs),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "build_spectrum",
                "lines.dat",
                "-o",
                "spectrum.dat",
                "--fwhm",
                "20.0",
                "--min",
                "0.0",
                "--max",
                "1000.0",
                "--step",
                "1.0",
                "--lorentzian",
                "--log-file",
                "off",
            ],
        )

        build_spectrum_cli.main()

        assert called == [
            {
                "input_file": "lines.dat",
                "output": "spectrum.dat",
                "alpha": None,
                "fwhm": 20.0,
                "wavenumber_min": 0.0,
                "wavenumber_max": 1000.0,
                "wavenumber_step": 1.0,
                "kernel": "lorentzian",
                "mode": FileWritingMode("w"),
            }
        ]
