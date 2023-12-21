import pytest

from PQAnalysis.io import InputFileParser
from PQAnalysis.io.inputFileReader.formats import InputFileFormat, InputFileFormatError


class TestInputFileParser:
    @pytest.mark.parametrize("example_dir", ["inputFileReader"], indirect=False)
    def test__init__(self, test_with_data_dir):
        with pytest.raises(FileNotFoundError) as exception:
            InputFileParser("non-existent-file")
        assert str(
            exception.value) == "File non-existent-file not found."

        input_file = "input.in"

        input_file_parser = InputFileParser(input_file)
        assert input_file_parser.filename == input_file
        assert input_file_parser.format == InputFileFormat.PQANALYSIS

        input_file_parser = InputFileParser(input_file, "pqanalysis")
        assert input_file_parser.filename == input_file
        assert input_file_parser.format == InputFileFormat.PQANALYSIS

        input_file_parser = InputFileParser(input_file, "pimd-qmcf")
        assert input_file_parser.filename == input_file
        assert input_file_parser.format == InputFileFormat.PIMD_QMCF

        input_file_parser = InputFileParser(input_file, "qmcfc")
        assert input_file_parser.filename == input_file
        assert input_file_parser.format == InputFileFormat.QMCFC

        with pytest.raises(InputFileFormatError) as exception:
            InputFileParser(input_file, "non-existent-format")
        assert str(exception.value) == """
'non-existent-format' is not a valid InputFileFormat.
Possible values are: InputFileFormat.PQANALYSIS, InputFileFormat.PIMD_QMCF, InputFileFormat.QMCFC
or their case insensitive string representation: PQANALYSIS, PIMD-QMCF, QMCFC"""

    @pytest.mark.parametrize("example_dir", ["inputFileReader"], indirect=False)
    def test_parse(self, test_with_data_dir):
        input_file_parser = InputFileParser("input.in")
        input_dictionary = input_file_parser.parse()

        dict = {}
        dict["a"] = (2.0, "float", "1")
        dict["b"] = ("test", "str", "1")
        dict["aaa"] = (["adf", "3.0", "4.5", "True"], "list(str)", "5")
        dict["key"] = ([2.0, 355.0, 4.1233], "list(float)", "7-13")
        dict["bool"] = (True, "bool", "15")
        dict["myglob"] = (["input.in"], "glob", "17")

        assert input_dictionary.dict == dict

        input_file_parser = InputFileParser("input_pimd-qmcf.in", "pimd-qmcf")
        input_dictionary = input_file_parser.parse()

        dict = {}
        dict["a"] = (2.0, "float", "1")
        dict["b"] = ("test", "str", "1")
        dict["aaa"] = (["adf", "3.0", "4.5", "True"], "list(str)", "5")
        dict["bool"] = (True, "bool", "7")
        dict["myglob"] = (["input.in"], "glob", "9")

        assert input_dictionary.dict == dict
