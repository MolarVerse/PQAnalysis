"""
Tests for the adf CLI.
"""

import sys

from PQAnalysis.cli import adf as adf_cli
from PQAnalysis.cli.adf import ADFCLI
from PQAnalysis.traj import MDEngineFormat



class TestADFCLI:

    """
    Tests for ADFCLI.
    """

    def test_program_name(self):
        assert ADFCLI.program_name() == "adf"

    def test_main_dispatches_input_file(self, monkeypatch):
        called = []

        monkeypatch.setattr(
            adf_cli,
            "adf",
            lambda input_file, engine: called.append((input_file, engine)),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            ["adf", "input.in", "--log-file", "off"],
        )

        adf_cli.main()

        assert called == [("input.in", MDEngineFormat.PQ)]

    def test_main_dispatches_engine(self, monkeypatch):
        called = []

        monkeypatch.setattr(
            adf_cli,
            "adf",
            lambda input_file, engine: called.append((input_file, engine)),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            ["adf", "input.in", "--engine", "qmcfc", "--log-file", "off"],
        )

        adf_cli.main()

        assert called == [("input.in", MDEngineFormat.QMCFC)]

    def test_add_arguments_registers_input_file_and_engine(self):
        from PQAnalysis.cli._argument_parser import _ArgumentParser

        parser = _ArgumentParser()
        ADFCLI.add_arguments(parser)

        args = parser.parse_args(["input.in"])

        assert args.input_file == "input.in"
        assert args.engine == MDEngineFormat.PQ
