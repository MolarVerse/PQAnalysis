"""
A module to test the ADF input file reader.
"""

import pytest

from PQAnalysis.analysis.adf.adf_input_file_reader import ADFInputFileReader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError, PQFileNotFoundError

# import topology marker
from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception



class TestADFInputFileReader:

    @pytest.mark.parametrize("example_dir", ["adf"], indirect=False)
    def test__init__(self, test_with_data_dir):
        with pytest.raises(PQFileNotFoundError) as exception:
            ADFInputFileReader("not-a-file")
        assert str(exception.value) == "File not-a-file not found."

        filename = "input.in"

        reader = ADFInputFileReader(filename)
        assert reader.filename == filename
        assert reader.parser.filename == filename

    def test__init__type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("filename", 1.0, str),
            exception=PQTypeError,
            function=ADFInputFileReader,
            filename=1.0,
        )

    @pytest.mark.parametrize("example_dir", ["adf"], indirect=False)
    def test_read_full(self, test_with_data_dir):
        reader = ADFInputFileReader("input.in")
        reader.read()

        assert reader.traj_files == ["traj.xyz"]
        assert reader.reference_selection == "O"
        assert reader.target_selection == "H"
        assert reader.target_selection_2 == "H"
        assert reader.n_angle_bins == 180
        assert reader.delta_angle is None
        assert reader.out_file == "adf.out"

    @pytest.mark.parametrize("example_dir", ["adf"], indirect=False)
    def test_read_minimal_defaults_optional_to_none(self, test_with_data_dir):
        reader = ADFInputFileReader("input_min.in")
        reader.read()

        assert reader.target_selection_2 is None
        assert reader.n_angle_bins is None
        assert reader.delta_angle is None
        assert reader.r_min_1 is None
        assert reader.r_max_1 is None
        assert reader.r_min_2 is None
        assert reader.r_max_2 is None
        assert reader.log_file is None

    @pytest.mark.parametrize("example_dir", ["adf"], indirect=False)
    def test_read_gated(self, test_with_data_dir):
        reader = ADFInputFileReader("input_gated.in")
        reader.read()

        assert reader.delta_angle == 2.0
        assert reader.r_min_1 == 0.8
        assert reader.r_max_1 == 3.0
        assert reader.r_min_2 == 0.8
        assert reader.r_max_2 == 3.0
        assert reader.log_file == "adf.log"

    @pytest.mark.parametrize("example_dir", ["adf"], indirect=False)
    def test_missing_required_key_raises(self, test_with_data_dir):
        reader = ADFInputFileReader("input_missing_required.in")

        with pytest.raises(InputFileError):
            reader.read()
