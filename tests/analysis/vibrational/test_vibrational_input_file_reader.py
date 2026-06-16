"""
Tests for the vibrational analysis input-file reader.
"""

import pytest

from PQAnalysis.analysis.vibrational import VibrationalAnalysisInputFileReader
from PQAnalysis.analysis.vibrational.vibrational_input_file_reader import _parse_modes
from PQAnalysis.exceptions import PQKeyError
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

    @pytest.mark.parametrize(
        ("contents", "message"),
        [
            (
                "structure_file = h2o.rst\n"
                "hessian_file = hessian.dat\n"
                "out_file = wavenumbers.dat\n"
                "hessian_sign = sideways\n",
                "The hessian_sign key must be one of: auto, positive, negative, 1, -1.",
            ),
            (
                "structure_file = h2o.rst\n"
                "hessian_file = hessian.dat\n"
                "out_file = wavenumbers.dat\n"
                "modes = invalid\n",
                (
                    "The modes key must be one of: all, nonzero, positive, "
                    "an integer, a list of integers or a range."
                ),
            ),
            (
                "structure_file = h2o.rst\n"
                "hessian_file = hessian.dat\n"
                "out_file = wavenumbers.dat\n"
                "modes = [0]\n",
                "Mode numbers must be positive and one-based.",
            ),
            (
                "structure_file = h2o.rst\n"
                "hessian_file = hessian.dat\n"
                "out_file = wavenumbers.dat\n"
                "modes_temperature = 0.0\n",
                "The modes_temperature key must be positive.",
            ),
        ],
    )
    @pytest.mark.parametrize("example_dir", ["vibrational"], indirect=False)
    def test_invalid_mode_options(self, test_with_data_dir, contents, message):
        with open("invalid.in", "w", encoding="utf-8") as file:
            file.write(contents)

        reader = VibrationalAnalysisInputFileReader("invalid.in")

        with pytest.raises(InputFileError) as exception:
            reader.read()

        assert str(exception.value) == message

    @pytest.mark.parametrize("example_dir", ["vibrational"], indirect=False)
    def test_parse_mode_variants(self, test_with_data_dir):
        with open("integer.in", "w", encoding="utf-8") as file:
            file.write(
                "structure_file = h2o.rst\n"
                "hessian_file = hessian.dat\n"
                "out_file = wavenumbers.dat\n"
                "modes = 3\n"
            )

        reader = VibrationalAnalysisInputFileReader("integer.in")
        reader.read()
        assert reader.modes == [3]

        with open("range.in", "w", encoding="utf-8") as file:
            file.write(
                "structure_file = h2o.rst\n"
                "hessian_file = hessian.dat\n"
                "out_file = wavenumbers.dat\n"
                "modes = 2-4\n"
            )

        reader = VibrationalAnalysisInputFileReader("range.in")
        reader.read()
        assert reader.modes == [2, 3]

    def test_parse_modes_helper_variants(self, monkeypatch):

        class MissingKeyDict:

            def __getitem__(self, key):
                raise PQKeyError(key)

        assert _parse_modes(MissingKeyDict(), "modes") is None
        assert _parse_modes({"modes": (None, "None", None)}, "modes") is None
        assert _parse_modes({"modes": (4, "int", None)}, "modes") == [4]
        assert _parse_modes({"modes": ([4, 5], "list(int)", None)},
                            "modes") == [4, 5]
        assert _parse_modes({"modes": (range(2, 5), "range", None)},
                            "modes") == [2, 3, 4]

        with pytest.raises(InputFileError) as exception:
            _parse_modes({"modes": (1.0, "float", None)}, "modes")

        assert str(exception.value) == (
            "The modes key must be one of: all, nonzero, positive, "
            "an integer, a list of integers or a range."
        )

        monkeypatch.setattr(
            VibrationalAnalysisInputFileReader.logger,
            "error",
            lambda *args, **kwargs: None,
        )
        assert _parse_modes({"modes": (1.0, "float", None)}, "modes") is None
