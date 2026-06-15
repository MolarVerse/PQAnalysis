"""
Tests for the vibrations CLI.
"""

from PQAnalysis.cli.vibrations import VibrationsCLI



class TestVibrationsCLI:

    """
    Tests for VibrationsCLI.
    """

    def test_program_name(self):
        assert VibrationsCLI.program_name() == "vibrations"
