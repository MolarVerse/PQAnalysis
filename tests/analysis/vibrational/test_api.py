"""
Tests for the vibrational analysis API.
"""

from pathlib import Path

import pytest

from PQAnalysis.analysis.vibrational import vibrations

from .. import pytestmark  # pylint: disable=unused-import



class TestVibrationalAnalysisAPI:

    """
    Tests for the vibrational analysis API.
    """

    @pytest.mark.parametrize("example_dir", ["vibrational"], indirect=False)
    def test_vibrations(self, test_with_data_dir):
        vibrations("input.in")

        output = Path("wavenumbers.dat")
        normal_modes = Path("normal_modes.dat")

        assert output.is_file()
        assert normal_modes.is_file()
        assert output.read_text(
            encoding="utf-8"
        ).startswith("# Wavenumbers (cm-1)  Intensities (km mol-1)")

    @pytest.mark.parametrize("example_dir", ["vibrational"], indirect=False)
    def test_vibrations_with_mode_output(self, test_with_data_dir):
        vibrations("mode_output.in")

        output = Path("wavenumbers.dat")
        normal_modes = Path("normal_modes.dat")
        modes_file = Path("modes.xyz")
        mode_6 = Path("mode-6.xyz")

        assert output.is_file()
        assert normal_modes.is_file()
        assert modes_file.is_file()
        assert mode_6.is_file()
        assert Path("mode-1.xyz").exists() is False
        assert modes_file.read_text(encoding="utf-8").splitlines(
        )[1].startswith("Properties=species:S:1:pos:R:3:mode:R:3 mode=6")
        assert "frame=1/8" in mode_6.read_text(encoding="utf-8")
