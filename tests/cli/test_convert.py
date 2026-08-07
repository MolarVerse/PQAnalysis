"""
Tests for the analysis output converter CLI.
"""

import sys

from PQAnalysis.cli import convert as convert_cli
from PQAnalysis.cli import main as pqanalysis_cli
from PQAnalysis.cli.convert import ConvertCLI
from PQAnalysis.io.formats import FileWritingMode



class TestConvertCLI:

    """Tests for ConvertCLI."""

    def test_program_name(self):
        assert ConvertCLI.program_name() == "convert"

    def test_main_dispatches_repeated_outputs_and_plot_fields(
        self, monkeypatch
    ):
        called = []
        monkeypatch.setattr(
            convert_cli,
            "convert_analysis_output",
            lambda **kwargs: called.append(kwargs),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "convert",
                "rdf.dat",
                "-o",
                "rdf.csv",
                "--output",
                "rdf.xvg",
                "--x",
                "r_i",
                "--y",
                "g_r_i",
                "--y",
                "N_r_i",
                "--mode",
                "o",
                "--log-file",
                "off",
            ],
        )

        convert_cli.main()

        assert called == [
            {
                "input_file": "rdf.dat",
                "output_files": ["rdf.csv", "rdf.xvg"],
                "x_field": "r_i",
                "y_fields": ["g_r_i", "N_r_i"],
                "mode": FileWritingMode.OVERWRITE,
            }
        ]

    def test_convert_is_registered_on_main_cli(self, monkeypatch):
        called = []
        monkeypatch.setattr(
            convert_cli,
            "convert_analysis_output",
            lambda **kwargs: called.append(kwargs),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "pqanalysis",
                "--log-file",
                "off",
                "convert",
                "rdf.dat",
                "-o",
                "rdf.tsv",
            ],
        )

        pqanalysis_cli.main()

        assert called[0]["output_files"] == ["rdf.tsv"]
