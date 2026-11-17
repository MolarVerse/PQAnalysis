"""
Tests for the rdf CLI.
"""

import sys

import pytest

from PQAnalysis.cli import rdf as rdf_cli
from PQAnalysis.cli.rdf import RDFCLI
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.traj import MDEngineFormat



class TestRDFCLI:

    """Tests for RDFCLI."""

    def test_program_name(self):
        assert RDFCLI.program_name() == "rdf"

    def test_main_dispatches_repeated_exports(self, monkeypatch):
        called = []
        monkeypatch.setattr(
            rdf_cli,
            "rdf",
            lambda input_file, engine, **kwargs: called.append(
                (input_file, engine, kwargs)
            ),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "rdf",
                "rdf.in",
                "--export",
                "rdf.csv",
                "--export",
                "rdf.xvg",
                "--log-file",
                "off",
            ],
        )

        rdf_cli.main()

        assert called == [
            (
                "rdf.in",
                MDEngineFormat.PQ,
                {
                    "export_files": ["rdf.csv", "rdf.xvg"],
                    "mode": FileWritingMode("w"),
                },
            )
        ]

    def test_main_dispatches_overwrite_mode(self, monkeypatch):
        called = []
        monkeypatch.setattr(
            rdf_cli,
            "rdf",
            lambda input_file, engine, **kwargs: called.append(kwargs),
        )
        monkeypatch.setattr(
            sys, "argv", ["rdf", "rdf.in", "--mode", "o", "--log-file", "off"]
        )

        rdf_cli.main()

        assert called == [{"mode": FileWritingMode.OVERWRITE}]

    def test_help_explains_mode(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["rdf", "--help"])

        with pytest.raises(SystemExit):
            rdf_cli.main()

        help_text = " ".join(capsys.readouterr().out.split())
        assert "--mode" in help_text
        assert '"o" overwrites' in help_text
