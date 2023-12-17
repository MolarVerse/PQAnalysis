import pytest

from PQAnalysis.io.inputFileReader.inputFileParser import InputDictionary


class TestInputFileDictionary:
    def test__init__(self):
        dictionary = InputDictionary()
        assert dictionary.dict == {}

    def test__setitem__(self):
        dictionary = InputDictionary()
        dictionary["KeY"] = "value"
        assert dictionary.dict == {"key": "value"}

        with pytest.raises(KeyError) as exception:
            dictionary["KeY"] = "value"
        assert str(
            exception.value) == "\'Input file key \"key\" defined multiple times in input file.\'"

    def test__getitem__(self):
        dictionary = InputDictionary()
        dictionary["key"] = "value"
        assert dictionary["KeY"] == "value"

        with pytest.raises(KeyError) as exception:
            dictionary["non-existent-key"]
        assert str(
            exception.value) == "\'Input file key \"non-existent-key\" not defined in input file.\'"

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
