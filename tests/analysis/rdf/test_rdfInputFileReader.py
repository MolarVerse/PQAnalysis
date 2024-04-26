import pytest

from PQAnalysis.analysis.rdf.rdf_input_file_reader import RDFInputFileReader
from PQAnalysis.io.inputFileReader.exceptions import InputFileError


class TestRDFInputFileReader:
    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test__init__(self, test_with_data_dir):
        with pytest.raises(FileNotFoundError) as exception:
            RDFInputFileReader("not-a-file")
        assert str(exception.value) == "File not-a-file not found."

        filename = "input.in"

        reader = RDFInputFileReader(filename)
        assert reader.filename == filename
        assert reader.parser.filename == filename

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_read(self, test_with_data_dir):
        with pytest.raises(InputFileError) as exception:
            reader = RDFInputFileReader("required_keys_not_given.in")
            reader.read()
        assert str(exception.value) == f"Not all required keys were set in the input file! The required keys are: \
{RDFInputFileReader.required_keys}."
