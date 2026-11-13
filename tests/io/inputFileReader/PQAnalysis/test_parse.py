import pytest

from ... import pytestmark

import PQAnalysis.io.input_file_reader.pq_analysis._parse as parse

from PQAnalysis.io.input_file_reader import InputDictionary
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis.traj import MDEngineFormat



def test_parse_real():
    input_dict = InputDictionary()
    input_dict["a"] = (2.0, "float", "1")

    assert parse._parse_real(input_dict, "notAKey") == None
    assert parse._parse_real(input_dict, "a") == 2.0

    input_dict["b"] = (2, "int", "1")

    assert parse._parse_real(input_dict, "b") == 2.0

    input_dict["c"] = ("2.0", "str", "1")

    with pytest.raises(InputFileError) as exception:
        parse._parse_real(input_dict, "c")
    assert str(
        exception.value
    ) == "The \"c\" value has to be of float type - actually it is parsed as a str"



def test_parse_positive_real():
    input_dict = InputDictionary()
    input_dict["a"] = (2.0, "float", "1")

    assert parse._parse_positive_real(input_dict, "notAKey") == None
    assert parse._parse_positive_real(input_dict, "a") == 2.0

    input_dict["b"] = (2, "int", "1")

    assert parse._parse_positive_real(input_dict, "b") == 2.0

    input_dict["c"] = ("2.0", "str", "1")

    with pytest.raises(InputFileError) as exception:
        parse._parse_positive_real(input_dict, "c")
    assert str(
        exception.value
    ) == "The \"c\" value has to be of float type - actually it is parsed as a str"

    input_dict["d"] = (-2.0, "float", "1")

    with pytest.raises(InputFileError) as exception:
        parse._parse_positive_real(input_dict, "d")
    assert str(
        exception.value
    ) == "The \"d\" value has to be a positive real number - It actually is -2.0!"



def test_parse_files():
    input_dict = InputDictionary()
    input_dict["a"] = ("file1", "str", "1")

    assert parse._parse_files(input_dict, "notAKey") == None
    assert parse._parse_files(input_dict, "a") == ["file1"]

    input_dict["b"] = (["file1", "file2"], "list(str)", "1")

    assert parse._parse_files(input_dict, "b") == ["file1", "file2"]

    input_dict["c"] = (["file1", "file2"], "glob", "1")

    assert parse._parse_files(input_dict, "c") == ["file1", "file2"]

    input_dict["d"] = ("None", "None", "1")

    assert parse._parse_files(input_dict, "d") == None

    input_dict["e"] = (2, "int", "1")

    with pytest.raises(InputFileError) as exception:
        parse._parse_files(input_dict, "e")
    assert str(
        exception.value
    ) == "The \"e\" value has to be either a string, glob or a list of strings - actually it is parsed as a int"



def test_parse_int():
    input_dict = InputDictionary()
    input_dict["a"] = (2, "int", "1")

    assert parse._parse_int(input_dict, "notAKey") == None
    assert parse._parse_int(input_dict, "a") == 2

    input_dict["d"] = ("2.0", "str", "1")

    with pytest.raises(InputFileError) as exception:
        parse._parse_int(input_dict, "d")
    assert str(
        exception.value
    ) == "The \"d\" value has to be of int type - actually it is parsed as a str"



def test_parse_positive_int():
    input_dict = InputDictionary()
    input_dict["a"] = (2, "int", "1")

    assert parse._parse_positive_int(input_dict, "notAKey") == None
    assert parse._parse_positive_int(input_dict, "a") == 2

    input_dict["b"] = (-2, "int", "1")

    with pytest.raises(InputFileError) as exception:
        parse._parse_positive_int(input_dict, "b")
    assert str(
        exception.value
    ) == "The \"b\" value has to be a positive integer - It actually is -2!"



def test_parse_string():
    input_dict = InputDictionary()
    input_dict["a"] = ("string", "str", "1")

    assert parse._parse_string(input_dict, "notAKey") == None
    assert parse._parse_string(input_dict, "a") == "string"

    input_dict["b"] = (2, "int", "1")

    with pytest.raises(InputFileError) as exception:
        parse._parse_string(input_dict, "b")
    assert str(
        exception.value
    ) == "The \"b\" value has to be of string type - actually it is parsed as a int"



def test_parse_selection():
    input_dict = InputDictionary()
    input_dict["word"] = ("O", "str", "1")
    input_dict["index"] = (72, "int", "1")
    input_dict["indices"] = (range(72, 80), "range", "1")
    input_dict["step"] = (range(0, 9, 3), "range", "1")
    input_dict["float"] = (2.0, "float", "1")

    assert parse._parse_selection(input_dict, "notAKey") is None
    assert parse._parse_selection(input_dict, "word") == "O"
    assert parse._parse_selection(input_dict, "index") == "72"
    assert parse._parse_selection(input_dict, "indices") == "72-80"
    assert parse._parse_selection(input_dict, "step") == "0..3..9"

    with pytest.raises(InputFileError) as exception:
        parse._parse_selection(input_dict, "float")
    assert str(
        exception.value
    ) == "The \"float\" value has to be of string type - actually it is parsed as a float"



def test_parse_bool():
    input_dict = InputDictionary()
    input_dict["a"] = (True, "bool", "1")

    assert parse._parse_bool(input_dict, "notAKey") == None
    assert parse._parse_bool(input_dict, "a") == True

    input_dict["b"] = (False, "bool", "1")

    assert parse._parse_bool(input_dict, "b") == False

    input_dict["c"] = (2, "int", "1")

    with pytest.raises(InputFileError) as exception:
        parse._parse_bool(input_dict, "c")
    assert str(
        exception.value
    ) == "The \"c\" value has to be of bool type - actually it is parsed as a int"
