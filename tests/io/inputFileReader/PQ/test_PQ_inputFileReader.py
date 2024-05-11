import pytest

from filecmp import cmp as filecmp

from ... import pytestmark

from PQAnalysis.io.input_file_reader.pq.pq_input_file_reader import _increase_digit_string, _get_digit_string_from_filename
from PQAnalysis.io.input_file_reader import PQInputFileReader as InputFileReader
from PQAnalysis.io.input_file_reader.formats import InputFileFormat
from PQAnalysis.exceptions import PQValueError



class TestPQ_inputFileReader:

    def test__increase_digit_string(self):
        digit_string = "1"
        assert _increase_digit_string(digit_string) == "2"

        digit_string = "2"
        assert _increase_digit_string(digit_string) == "3"

        digit_string = "9"
        assert _increase_digit_string(digit_string) == "10"

        digit_string = "01"
        assert _increase_digit_string(digit_string) == "02"

        digit_string = "0004"
        assert _increase_digit_string(digit_string) == "0005"

        digit_string = "0009"
        assert _increase_digit_string(digit_string) == "0010"

        digit_string = "00.09"
        with pytest.raises(PQValueError) as exception:
            _increase_digit_string(digit_string)
        assert str(
            exception.value
        ) == "digit_string 00.09 contains non-digit characters."

        digit_string = ""
        assert _increase_digit_string(digit_string) == "1"

    def test__get_digit_string_from_filename(self):
        filename = "input.in"
        with pytest.raises(PQValueError) as exception:
            _get_digit_string_from_filename(filename)
        assert str(
            exception.value
        ) == "Filename input.in does not contain a number to be continued from. It has to be of the form \"...<number>.<extension>\"."

        filename = "input_0001.in"
        assert _get_digit_string_from_filename(filename) == "0001"

        filename = "input_001.in"
        assert _get_digit_string_from_filename(filename) == "001"

        filename = "input_099.in.asdf"
        assert _get_digit_string_from_filename(filename) == "099"

        filename = "input_00.in"
        assert _get_digit_string_from_filename(filename) == "00"

    @pytest.mark.parametrize(
        "example_dir",
        ["inputFileReader/PQ_input/"],
        indirect=False
    )
    def test__init__(self, test_with_data_dir):
        input_file_reader = InputFileReader("run-08.in")

        assert input_file_reader.filename == "run-08.in"
        assert input_file_reader.format == InputFileFormat("PQ")
        assert input_file_reader.parser.filename == "run-08.in"
        assert input_file_reader.parser.input_format == InputFileFormat("PQ")

    @pytest.mark.parametrize(
        "example_dir",
        ["inputFileReader/PQ_input/"],
        indirect=False
    )
    def test_read(self, test_with_data_dir):
        input_file_reader = InputFileReader("run-08.in")

        input_dictionary = input_file_reader.read()

        assert input_file_reader.raw_input_file == open("run-08.in",
            "r").read()
        assert input_file_reader.start_file == "md-07.rst"
        assert input_file_reader.is_rpmd_start_file_defined == False
        assert input_file_reader.restart_file == "md-08.rst"
        assert input_file_reader.trajectory_file == "md-08.xyz"
        assert input_file_reader.velocity_file == "md-08.vel"
        assert input_file_reader.force_file == "md-08.frc"
        assert input_file_reader.charge_file == "md-08.chrg"
        assert input_file_reader.energy_file == "md-08.en"
        assert input_file_reader.output_file == "md-08.out"
        assert input_file_reader.info_file == "md-08.info"
        assert input_file_reader.file_prefix == "md-08"

        input_file_reader = InputFileReader("run-08.rpmd.in")

        input_dictionary = input_file_reader.read()

        assert input_file_reader.raw_input_file == open("run-08.rpmd.in",
            "r").read()
        assert input_file_reader.start_file == "md-07.rst"
        assert input_file_reader.is_rpmd_start_file_defined == True
        assert input_file_reader.rpmd_start_file == "md-07.rpmd.rst"
        assert input_file_reader.rpmd_restart_file == "md-08.rst"
        assert input_file_reader.rpmd_trajectory_file == "md-08.xyz"
        assert input_file_reader.rpmd_velocity_file == "md-08.vel"
        assert input_file_reader.rpmd_force_file == "md-08.frc"
        assert input_file_reader.rpmd_charge_file == "md-08.chrg"
        assert input_file_reader.rpmd_energy_file == "md-08.en"
        assert input_file_reader.output_file == "md-08.out"
        assert input_file_reader.info_file == "md-08.info"

        # testing if ValueError is raised if no start file is defined
        input_file_reader = InputFileReader("no_start_file.in")

        with pytest.raises(PQValueError) as exception:
            input_file_reader.read()
        assert str(
            exception.value
        ) == "No start file defined in input file no_start_file.in."

    @pytest.mark.parametrize(
        "example_dir",
        ["inputFileReader/PQ_input/"],
        indirect=False
    )
    def test__parse_start_n(self, test_with_data_dir):
        input_file_reader = InputFileReader("run-08.in")
        input_file_reader.read()

        assert input_file_reader._parse_start_n() == "07"

        input_file_reader = InputFileReader("run-08.rpmd.in")
        input_file_reader.read()

        assert input_file_reader._parse_start_n() == "07"

        input_file_reader = InputFileReader("n_not_matching.in")
        input_file_reader.read()

        with pytest.raises(PQValueError) as exception:
            input_file_reader._parse_start_n()
        assert str(
            exception.value
        ) == "N from start_file (07) and rpmd_start_file (08) do not match."

    @pytest.mark.parametrize(
        "example_dir",
        ["inputFileReader/PQ_input/"],
        indirect=False
    )
    def test__parse_actual_n(self, test_with_data_dir):
        input_file_reader = InputFileReader("run-08.in")
        input_file_reader.read()

        assert input_file_reader._parse_actual_n() == "08"

        input_file_reader = InputFileReader("run-08.rpmd.in")
        input_file_reader.read()

        assert input_file_reader._parse_actual_n() == "08"

        # testing if ValueError is raised if n not consistent in output files

        input_file_reader = InputFileReader("n_not_matching.in")
        input_file_reader.read()

        with pytest.raises(PQValueError) as exception:
            input_file_reader._parse_actual_n()
        assert str(
            exception.value
        ) == "Actual n in output files is not consistent."

        # testing if ValueError is raised if no output file is defined

        input_file_reader = InputFileReader("no_output_files.in")
        input_file_reader.read()

        with pytest.raises(PQValueError) as exception:
            input_file_reader._parse_actual_n()
        assert str(
            exception.value
        ) == "No output file found to determine actual n."

    @pytest.mark.parametrize(
        "example_dir",
        ["inputFileReader/PQ_input/"],
        indirect=False
    )
    def test_continue_input_file(self, test_with_data_dir):
        # testing if ValueError is raised if actual n and input file n do not match
        input_file_reader = InputFileReader(
            "n_not_matching_input_file_n-08.in"
        )
        input_file_reader.read()

        with pytest.raises(PQValueError) as exception:
            input_file_reader.continue_input_file(2)
        assert str(
            exception.value
        ) == "Actual n (09) and input file n (08) do not match."

        # testing if ValueError is raised if old n is not exactly one less than actual n
        input_file_reader = InputFileReader(
            "old_n_not_less_one_than_actual_n-09.in"
        )
        input_file_reader.read()

        with pytest.raises(PQValueError) as exception:
            input_file_reader.continue_input_file(2)
        assert str(
            exception.value
        ) == "Old n (07) has to be one less than actual n (09)."

        input_file_reader = InputFileReader("run-08.in")
        input_file_reader.read()
        input_file_reader.continue_input_file(2)

        assert filecmp("run-09.in", "run-09.in.ref")
        assert filecmp("run-10.in", "run-10.in.ref")

        input_file_reader = InputFileReader("run-08.rpmd.in")
        input_file_reader.read()
        input_file_reader.continue_input_file(2)

        assert filecmp("run-09.rpmd.in", "run-09.rpmd.in.ref")
        assert filecmp("run-10.rpmd.in", "run-10.rpmd.in.ref")
