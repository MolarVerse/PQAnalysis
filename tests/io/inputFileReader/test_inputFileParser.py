import pytest

from .. import pytestmark

from PQAnalysis.io import InputFileParser
from PQAnalysis.io.input_file_reader.formats import InputFileFormat, InputFileFormatError
from PQAnalysis.exceptions import PQFileNotFoundError



class TestInputFileParser:

    @pytest.mark.parametrize(
        "example_dir",
        ["inputFileReader"],
        indirect=False
    )
    def test__init__(self, test_with_data_dir):
        with pytest.raises(PQFileNotFoundError) as exception:
            InputFileParser("non-existent-file")
        assert str(exception.value) == "File non-existent-file not found."

        input_file = "input.in"

        input_file_parser = InputFileParser(input_file)
        assert input_file_parser.filename == input_file
        assert input_file_parser.input_format == InputFileFormat.PQANALYSIS

        input_file_parser = InputFileParser(input_file, "pqanalysis")
        assert input_file_parser.filename == input_file
        assert input_file_parser.input_format == InputFileFormat.PQANALYSIS

        input_file_parser = InputFileParser(input_file, "PQ")
        assert input_file_parser.filename == input_file
        assert input_file_parser.input_format == InputFileFormat.PQ

        input_file_parser = InputFileParser(input_file, "qmcfc")
        assert input_file_parser.filename == input_file
        assert input_file_parser.input_format == InputFileFormat.QMCFC

        with pytest.raises(InputFileFormatError) as exception:
            InputFileParser(input_file, "non-existent-format")
        assert str(exception.value) == (
            "\n"
            "'non-existent-format' is not a valid InputFileFormat.\n"
            "Possible values are: InputFileFormat.PQANALYSIS, InputFileFormat.PQ, "
            "InputFileFormat.QMCFC or their case insensitive string representation: "
            "PQANALYSIS, PQ, QMCFC"
        )

    @pytest.mark.parametrize(
        "example_dir",
        ["inputFileReader"],
        indirect=False
    )
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

        input_file_parser = InputFileParser("input_PQ.in", "PQ")
        input_dictionary = input_file_parser.parse()

        dict = {}
        dict["a"] = (2.0, "float", "1")
        dict["b"] = ("test", "str", "1")
        dict["aaa"] = (["adf", "3.0", "4.5", "True"], "list(str)", "5")
        dict["bool"] = (True, "bool", "7")
        dict["myglob"] = (["input.in"], "glob", "9")

        assert input_dictionary.dict == dict

    def test_parse_boolean_false(self, tmp_path):
        input_file = tmp_path / "input.in"
        input_file.write_text("enabled = True\ndisabled = False\n", encoding="utf-8")

        input_file_parser = InputFileParser(str(input_file))
        input_dictionary = input_file_parser.parse()

        assert input_dictionary["enabled"] == (True, "bool", "1")
        assert input_dictionary["disabled"] == (False, "bool", "2")

    def test_parse_qmcfc_selectors(self, tmp_path):
        input_file = tmp_path / "qmcfc.in"
        input_file.write_text(
            "jobtype = qmcf-md;\n"
            "nstep = 5000; timestep = 0.2; omega = 3E13;\n"
            "solute_charge = +2.0;\n"
            "qm_center = 8:1;\n"
            "qm_blacklist = 1-6, 9-13, 15-16, 18, 20-48;\n"
            "qm_whitelist = 7, 14, 17, 19;\n",
            encoding="utf-8"
        )

        input_dictionary = InputFileParser(str(input_file), "qmcfc").parse()

        assert input_dictionary["jobtype"] == ("qmcf-md", "str", "1")
        assert input_dictionary["nstep"] == (5000, "int", "2")
        assert input_dictionary["timestep"] == (0.2, "float", "2")
        assert input_dictionary["omega"] == (3e13, "float", "2")
        assert input_dictionary["solute_charge"] == (2.0, "float", "3")
        assert input_dictionary["qm_center"] == ("8:1", "str", "4")
        assert input_dictionary["qm_blacklist"] == (
            ["1-6", "9-13", "15-16", "18", "20-48"],
            "list(str)",
            "5"
        )
        assert input_dictionary["qm_whitelist"] == (
            ["7", "14", "17", "19"],
            "list(str)",
            "6"
        )

    def test_parse_qmcfc_latin1_comments(self, tmp_path):
        input_file = tmp_path / "qmcfc.in"
        input_file.write_bytes(b"# force constant in A\xb2*kcal/mol\njobtype = mm-md;\n")

        parser = InputFileParser(str(input_file), "qmcfc")
        input_dictionary = parser.parse()

        assert "A\u00b2*kcal/mol" in parser.raw_input_file
        assert input_dictionary["jobtype"] == ("mm-md", "str", "2")
