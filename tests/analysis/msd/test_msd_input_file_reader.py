"""
Tests for the MSDInputFileReader class.
"""

from pathlib import Path

import pytest

from PQAnalysis.analysis.msd.msd_input_file_reader import MSDInputFileReader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError, PQFileNotFoundError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging, assert_logging_with_exception



class TestMSDInputFileReader:

    """
    Tests for the MSDInputFileReader class.
    """

    @pytest.mark.parametrize("example_dir", ["msd"], indirect=False)
    def test__init__(self, test_with_data_dir):
        with pytest.raises(PQFileNotFoundError) as exception:
            MSDInputFileReader("not-a-file")
        assert str(exception.value) == "File not-a-file not found."

        filename = "input.in"

        reader = MSDInputFileReader(filename)
        assert reader.filename == filename
        assert reader.parser.filename == filename

    def test__init__type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("filename", 1.0, str),
            exception=PQTypeError,
            function=MSDInputFileReader,
            filename=1.0,
        )

    @pytest.mark.parametrize("example_dir", ["msd"], indirect=False)
    def test_read_defaults(self, test_with_data_dir):
        reader = MSDInputFileReader("input_defaults.in")

        reader.read()

        assert reader.traj_files == ["traj.xyz"]
        assert reader.target_selection == "O"
        assert reader.out_file == "msd.dat"
        assert reader.window is None
        assert reader.gap is None
        assert reader.n_start is None
        assert reader.time_step is None
        assert reader.fit_window is None
        assert reader.log_file is None
        assert reader.use_full_atom_info is None

    @pytest.mark.parametrize("example_dir", ["msd"], indirect=False)
    def test_read_full(self, test_with_data_dir):
        reader = MSDInputFileReader("input_full.in")

        reader.read()

        assert reader.traj_files == ["traj.xyz"]
        assert reader.target_selection == "O"
        assert reader.out_file == "msd.dat"
        assert reader.log_file == "msd.log"
        assert reader.window == 50
        assert reader.gap == 5
        assert reader.n_start == 7
        assert reader.time_step == 0.25
        assert reader.fit_window == 12
        assert reader.use_full_atom_info is True

    @pytest.mark.parametrize("example_dir", ["msd"], indirect=False)
    def test_start_alias(self, test_with_data_dir):
        reader = MSDInputFileReader("input_start_alias.in")

        reader.read()

        assert reader.n_start == 7

    @pytest.mark.parametrize("example_dir", ["msd"], indirect=False)
    def test_alias_conflict(self, test_with_data_dir, caplog):
        reader = MSDInputFileReader("input_alias_conflict.in")

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSDInputFileReader",
            logging_level="ERROR",
            message_to_test=(
                "The keys 'first_frame' and 'start' are aliases and "
                "cannot be used at the same time."
            ),
            exception=InputFileError,
            function=reader.read,
        )

    @pytest.mark.parametrize("example_dir", ["msd"], indirect=False)
    def test_negative_start(self, test_with_data_dir, caplog):
        # the input file grammar parses negative numbers as strings,
        # so a negative start value is already rejected by the
        # generic int parsing of the base reader
        reader = MSDInputFileReader("input_negative_start.in")

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="PQAnalysisInputFileReader",
            logging_level="ERROR",
            message_to_test=(
                'The "start" value has to be of int type - '
                "actually it is parsed as a str"
            ),
            exception=InputFileError,
            function=reader.read,
        )

    @pytest.mark.parametrize("example_dir", ["msd"], indirect=False)
    def test_zero_time_step(self, test_with_data_dir, caplog):
        # the generic positive-real parsing of the base reader
        # accepts 0.0, so the MSD reader has to reject it itself
        Path("input_zero_time_step.in").write_text(
            (
                "traj_files = traj.xyz\n"
                "target_selection = O\n"
                "out_file = msd.dat\n"
                "time_step = 0.0\n"
            ),
            encoding="utf-8",
        )

        reader = MSDInputFileReader("input_zero_time_step.in")

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSDInputFileReader",
            logging_level="ERROR",
            message_to_test=(
                "The 'time_step' value has to be a positive real "
                "number - It actually is 0.0!"
            ),
            exception=InputFileError,
            function=reader.read,
        )

    @pytest.mark.parametrize("example_dir", ["msd"], indirect=False)
    def test_missing_required_key(self, test_with_data_dir, caplog):
        reader = MSDInputFileReader("input_missing_required.in")

        assert_logging_with_exception(
            caplog=caplog,
            logging_name="MSDInputFileReader",
            logging_level="ERROR",
            message_to_test=(
                "Not all required keys were set in the input file! "
                f"The required keys are: {MSDInputFileReader.required_keys}."
            ),
            exception=InputFileError,
            function=reader.read,
        )

    @pytest.mark.parametrize("example_dir", ["msd"], indirect=False)
    def test_unknown_key_warning(self, test_with_data_dir, caplog):
        Path("input_unknown_key.in").write_text(
            (
                "traj_files = traj.xyz\n"
                "target_selection = O\n"
                "out_file = msd.dat\n"
                "r_max = 5.0\n"
            ),
            encoding="utf-8",
        )

        reader = MSDInputFileReader("input_unknown_key.in")

        assert_logging(
            caplog=caplog,
            logging_name="MSDInputFileReader",
            logging_level="WARNING",
            message_to_test=(
                "Unknown keys were set in the input file! The known keys "
                f"are: {MSDInputFileReader.required_keys + MSDInputFileReader.optional_keys}. "
                "They will be ignored!"
            ),
            function=reader.read,
        )
