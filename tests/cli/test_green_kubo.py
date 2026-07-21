"""
Tests for the green_kubo CLI.
"""

import sys

from PQAnalysis.cli import green_kubo as green_kubo_cli
from PQAnalysis.cli.green_kubo import GreenKuboCLI
from PQAnalysis.traj import MDEngineFormat



class TestGreenKuboCLI:

    """
    Tests for GreenKuboCLI.
    """

    def test_program_name(self):
        """
        The program name is green_kubo.
        """
        assert GreenKuboCLI.program_name() == "green_kubo"

    def test_main_dispatches_defaults(self, monkeypatch):
        """
        main() dispatches the positional input file and the default
        engine to the green_kubo API function.
        """
        called = []

        monkeypatch.setattr(
            green_kubo_cli,
            "green_kubo",
            lambda *args: called.append(args),
        )
        monkeypatch.setattr(sys, "argv", ["green_kubo", "gk.in"])

        green_kubo_cli.main()

        assert called == [("gk.in", MDEngineFormat.PQ)]

    def test_main_dispatches_engine(self, monkeypatch):
        """
        main() dispatches the selected engine to the green_kubo API
        function.
        """
        called = []

        monkeypatch.setattr(
            green_kubo_cli,
            "green_kubo",
            lambda *args: called.append(args),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            ["green_kubo", "gk.in", "--engine", "qmcfc"],
        )

        green_kubo_cli.main()

        assert called == [("gk.in", MDEngineFormat.QMCFC)]
