"""
Tests for the msd CLI.
"""

import sys

from PQAnalysis.cli import msd as msd_cli
from PQAnalysis.cli.msd import MSDCLI
from PQAnalysis.traj import MDEngineFormat



class TestMSDCLI:

    """
    Tests for MSDCLI.
    """

    def test_program_name(self):
        assert MSDCLI.program_name() == "msd"

    def test_main_dispatches_input_file(self, monkeypatch):
        called = []

        monkeypatch.setattr(
            msd_cli,
            "msd",
            lambda input_file, engine: called.append((input_file, engine)),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            ["msd", "input.in", "--log-file", "off"],
        )

        msd_cli.main()

        assert called == [("input.in", MDEngineFormat.PQ)]

    def test_main_dispatches_engine(self, monkeypatch):
        called = []

        monkeypatch.setattr(
            msd_cli,
            "msd",
            lambda input_file, engine: called.append((input_file, engine)),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            ["msd", "input.in", "--engine", "qmcfc", "--log-file", "off"],
        )

        msd_cli.main()

        assert called == [("input.in", MDEngineFormat.QMCFC)]
