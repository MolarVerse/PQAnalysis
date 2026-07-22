"""
Tests for the vacf CLI.
"""

import sys

from PQAnalysis.cli import vacf as vacf_cli
from PQAnalysis.cli.vacf import VACFCLI
from PQAnalysis.traj import MDEngineFormat



class TestVACFCLI:

    """
    Tests for VACFCLI.
    """

    def test_program_name(self):
        """
        The program name is vacf.
        """
        assert VACFCLI.program_name() == "vacf"

    def test_main_dispatches_defaults(self, monkeypatch):
        """
        main() dispatches the positional input file and the default
        engine to the vacf API function.
        """
        called = []

        monkeypatch.setattr(
            vacf_cli,
            "vacf",
            lambda *args: called.append(args),
        )
        monkeypatch.setattr(sys, "argv", ["vacf", "vacf.in"])

        vacf_cli.main()

        assert called == [("vacf.in", MDEngineFormat.PQ)]

    def test_main_dispatches_engine(self, monkeypatch):
        """
        main() dispatches the selected engine to the vacf API function.
        """
        called = []

        monkeypatch.setattr(
            vacf_cli,
            "vacf",
            lambda *args: called.append(args),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            ["vacf", "vacf.in", "--engine", "qmcfc"],
        )

        vacf_cli.main()

        assert called == [("vacf.in", MDEngineFormat.QMCFC)]
