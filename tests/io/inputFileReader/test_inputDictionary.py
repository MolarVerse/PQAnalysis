import pytest

from .. import pytestmark

from PQAnalysis.io.input_file_reader.input_file_parser import InputDictionary
from PQAnalysis.exceptions import PQKeyError



class TestInputFileDictionary:

    def test__init__(self):
        dictionary = InputDictionary()
        assert dictionary.dict == {}

    def test__setitem__(self):
        dictionary = InputDictionary()
        dictionary["KeY"] = ("value", "type", "line")
        assert dictionary.dict == {"key": ("value", "type", "line")}

        with pytest.raises(PQKeyError) as exception:
            dictionary["KeY"] = ("value", "type", "line")
        assert str(
            exception.value
        ) == "Input file key \"key\" defined multiple times in input file."

    def test__getitem__(self):
        dictionary = InputDictionary()
        dictionary["key"] = ("value", "type", "line")
        assert dictionary["KeY"] == ("value", "type", "line")

        with pytest.raises(PQKeyError) as exception:
            dictionary["non-existent-key"]
        assert str(
            exception.value
        ) == "Input file key \"non-existent-key\" not defined in input file."

    def test_get_value(self):
        dictionary = InputDictionary()
        dictionary["key"] = ("value", "type", "1")
        assert dictionary.get_value("KeY") == "value"

    def test_get_line(self):
        dictionary = InputDictionary()
        dictionary["key"] = ("value", "type", "1")
        assert dictionary.get_line("KeY") == "1"

    def test_get_type(self):
        dictionary = InputDictionary()
        dictionary["key"] = ("value", "type", "1")
        assert dictionary.get_type("KeY") == "type"

    def test_keys(self):
        dictionary = InputDictionary()
        dictionary["key"] = ("value", "type", "1")
        assert dictionary.keys() == ["key"]

        dictionary["key2"] = ("value", "type", "1")
        assert dictionary.keys() == ["key", "key2"]

    def test__eq__(self):
        dictionary = InputDictionary()
        dictionary["key"] = ("value", "type", "1")

        dictionary2 = InputDictionary()
        dictionary2["key"] = ("value", "type", "1")

        assert dictionary == dictionary
        assert dictionary == dictionary2
        assert dictionary != "not a dictionary"
