"""
Tests for the vibrational analysis input-file reader.
"""

import pytest

from PQAnalysis.analysis.vibrational import VibrationalAnalysisInputFileReader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError

from .. import pytestmark  # pylint: disable=unused-import



class TestVibrationalAnalysisInputFileReader:

    """
    Tests for VibrationalAnalysisInputFileReader.
    """

    @pytest.mark.parametrize("example_dir", ["vibrational"], indirect=False)
    def test_read(self, test_with_data_dir):
        reader = VibrationalAnalysisInputFileReader("input.in")
        reader.read()

        assert reader.structure_file == "h2o.rst"
        assert reader.hessian_file == "hessian.dat"
        assert reader.moldescriptor_file == "moldescriptor.dat"
        assert reader.out_file == "wavenumbers.dat"
        assert reader.unit == "kcal"
        assert reader.hessian_sign == "auto"
        assert reader.normal_modes_file == "normal_modes.dat"
        assert reader.modes_prefix is None
        assert reader.modes_file is None
        assert reader.modes == "all"
        assert reader.modes_frames == 30
        assert reader.modes_amplitude == 0.25
        assert reader.modes_temperature is None
        assert reader.modes_threshold == 1.0e-8

    @pytest.mark.parametrize("example_dir", ["vibrational"], indirect=False)
    def test_read_mode_output(self, test_with_data_dir):
        reader = VibrationalAnalysisInputFileReader("mode_output.in")
        reader.read()

        assert reader.modes_prefix == "mode"
        assert reader.modes_file == "modes.xyz"
        assert reader.modes == [6, 7, 9]
        assert reader.modes_frames == 8
        assert reader.modes_amplitude == 0.3
        assert reader.modes_temperature == 300.0
        assert reader.modes_threshold == 0.5

    @pytest.mark.parametrize("example_dir", ["vibrational"], indirect=False)
    def test_invalid_unit(self, test_with_data_dir):
        reader = VibrationalAnalysisInputFileReader("invalid_unit.in")

        with pytest.raises(InputFileError) as exception:
            reader.read()

        assert str(
            exception.value
        ) == "The unit key must be one of: kcal, hartree, ev."
