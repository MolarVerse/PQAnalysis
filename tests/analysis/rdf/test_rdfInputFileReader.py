import pytest

from PQAnalysis.analysis.rdf.rdf_input_file_reader import RDFInputFileReader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis.type_checking import get_type_error_message
from PQAnalysis.exceptions import PQTypeError, PQFileNotFoundError

# import topology marker
from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception



class TestRDFInputFileReader:

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test__init__(self, test_with_data_dir):
        with pytest.raises(PQFileNotFoundError) as exception:
            RDFInputFileReader("not-a-file")
        assert str(exception.value) == "File not-a-file not found."

        filename = "input.in"

        reader = RDFInputFileReader(filename)
        assert reader.filename == filename
        assert reader.parser.filename == filename

    def test__init__type_checking(self, caplog):
        assert_logging_with_exception(
            caplog=caplog,
            logging_name="TypeChecking",
            logging_level="ERROR",
            message_to_test=get_type_error_message("filename",
            1.0,
            str),
            exception=PQTypeError,
            function=RDFInputFileReader,
            filename=1.0,
        )

    @pytest.mark.parametrize("example_dir", ["rdf"], indirect=False)
    def test_read(self, test_with_data_dir):
        with pytest.raises(InputFileError) as exception:
            reader = RDFInputFileReader("required_keys_not_given.in")
            reader.read()
        assert str(
            exception.value
        ) == f"Not all required keys were set in the input file! The required keys are: \
{RDFInputFileReader.required_keys}."
