"""
Tests for the rdf CLI.
"""

import sys

from PQAnalysis.cli import rdf as rdf_cli
from PQAnalysis.cli.rdf import RDFCLI
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
                    "export_files": ["rdf.csv", "rdf.xvg"]
                },
            )
        ]
