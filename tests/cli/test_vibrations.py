"""
Tests for the vibrations CLI.
"""

import sys

from PQAnalysis.cli import vibrations as vibrations_cli
from PQAnalysis.cli.vibrations import VibrationsCLI



class TestVibrationsCLI:

    """
    Tests for VibrationsCLI.
    """

    def test_program_name(self):
        assert VibrationsCLI.program_name() == "vibrations"

    def test_main_dispatches_input_file(self, monkeypatch):
        called = []

        monkeypatch.setattr(
            vibrations_cli,
            "vibrations",
            lambda input_file: called.append(input_file),
        )
        monkeypatch.setattr(
            sys,
            "argv",
            ["vibrations", "input.in", "--log-file", "off"],
        )

        vibrations_cli.main()

        assert called == ["input.in"]
