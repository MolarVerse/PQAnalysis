import pytest
import numpy as np

from PQAnalysis.io.base import PQFileNotFoundError
from tests.io.inputFileReader.PQ.test_PQ_inputFileReader import InputFileReader

from ... import pytestmark

from PQAnalysis.io.input_file_reader.pq_analysis.pqanalysis_input_file_reader import PQAnalysisInputFileReader
from PQAnalysis.io.input_file_reader.formats import InputFileFormat

from PQAnalysis.io.input_file_reader.exceptions import InputFileError, InputFileWarning, InputFileFormatError
from PQAnalysis.exceptions import PQValueError



class TestPQAnalysisInputFileReader:

    @pytest.mark.parametrize(
        "example_dir", ["inputFileReader/PQAnalysis_input/"], indirect=False
    )
    def test__init__(self, test_with_data_dir):
        pqanalysis_input_file_reader = PQAnalysisInputFileReader("input_PQ.in")

        assert pqanalysis_input_file_reader.filename == "input_PQ.in"
        assert pqanalysis_input_file_reader.format == InputFileFormat.PQANALYSIS
        assert pqanalysis_input_file_reader.parser.filename == "input_PQ.in"
        assert pqanalysis_input_file_reader.parser.input_format == InputFileFormat.PQANALYSIS

        assert pqanalysis_input_file_reader.dictionary is None
        assert pqanalysis_input_file_reader.raw_input_file is None

    def test__init__no_input_file(self):
        with pytest.raises(PQFileNotFoundError) as exception:
            pqanalysis_input_file_reader = PQAnalysisInputFileReader(
                "no_input_file.in"
            )
        assert str(exception.value) == f"File no_input_file.in not found."

    @pytest.mark.parametrize(
        "example_dir", ["inputFileReader/PQAnalysis_input/"], indirect=False
    )
    def test_read(self, test_with_data_dir):
        pqanalysis_input_file_reader = PQAnalysisInputFileReader("input_PQ.in")

        pqanalysis_input_file_dictonary = pqanalysis_input_file_reader.read()

        assert pqanalysis_input_file_reader.raw_input_file == open(
            "input_PQ.in", "r"
        ).read()
        assert pqanalysis_input_file_reader.dictionary == pqanalysis_input_file_reader.parser.parse(
        )

    @pytest.mark.parametrize(
        "example_dir", ["inputFileReader/PQAnalysis_input/"], indirect=False
    )
    def test_check_required_keys(self, test_with_data_dir):
        pqanalysis_input_file_reader = PQAnalysisInputFileReader("input_PQ.in")
        pqanalysis_input_file_reader.read()

        required_keys = ["required_key1", "required_key2"]

        assert pqanalysis_input_file_reader.check_required_keys(
            required_keys
        ) == None

        required_keys = ["required_key1", "required_key2", "required_key3"]
        with pytest.raises(InputFileError) as exception:
            pqanalysis_input_file_reader.check_required_keys(required_keys)
        assert str(exception.value) == (
            "Not all required keys were set in the input file! The required keys are: ['required_key1', 'required_key2', 'required_key3']."
        )

    @pytest.mark.parametrize(
        "example_dir", ["inputFileReader/PQAnalysis_input/"], indirect=False
    )
    def test_check_known_keys(self, test_with_data_dir):
        pqanalysis_input_file_reader = PQAnalysisInputFileReader("input_PQ.in")
        pqanalysis_input_file_reader.read()

        known_keys = [
            "required_key1", "required_key2", "optional_key1", "optional_key2"
        ]

        assert pqanalysis_input_file_reader.check_known_keys(
            known_keys
        ) == None

    @pytest.mark.parametrize(
        "example_dir", ["inputFileReader/PQAnalysis_input/"], indirect=False
    )
    def test_check_known_keys_no_known_keys(self, test_with_data_dir):
        pqanalysis_input_file_reader = PQAnalysisInputFileReader("input_PQ.in")
        pqanalysis_input_file_reader.read()
        known_keys = None
        assert pqanalysis_input_file_reader.check_known_keys(
            known_keys
        ) == None

    @pytest.mark.parametrize(
        "example_dir", ["inputFileReader/PQAnalysis_input/"], indirect=False
    )
    def test_not_defined_optional_keys(self, test_with_data_dir):
        pqanalysis_input_file_reader = PQAnalysisInputFileReader("input_PQ.in")
        pqanalysis_input_file_reader.read()

        optional_keys = ["optional_key1", "optional_key2", "optional_key3"]

        assert pqanalysis_input_file_reader.not_defined_optional_keys(
            optional_keys
        ) == None
        assert pqanalysis_input_file_reader.dictionary["optional_key3"] == (
            None, "None", "None"
        )
