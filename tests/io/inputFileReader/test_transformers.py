import pytest

from .. import pytestmark

from lark import Token

from PQAnalysis.io.input_file_reader.input_file_parser import (
    PrimitiveTransformer,
    ComposedDatatypesTransformer
)
from PQAnalysis.exceptions import PQTypeError



class TestPrimitiveTransformer:

    def test__init__(self):
        transformer = PrimitiveTransformer()
        assert transformer.__visit_tokens__ == False

        transformer = PrimitiveTransformer(True)
        assert transformer.__visit_tokens__ == True

        transformer = PrimitiveTransformer(False)
        assert transformer.__visit_tokens__ == False

    def test_float(self):
        transformer = PrimitiveTransformer()
        token = Token("FLOAT", "1.0")
        token.end_line = 1

        assert transformer.float([token]) == (1.0, "float", "1")

    def test_int(self):
        transformer = PrimitiveTransformer()
        token = Token("INT", "1")
        token.end_line = 1

        assert transformer.integer([token]) == (1, "int", "1")

    def test_word(self):
        transformer = PrimitiveTransformer()
        token = Token("WORD", "word")
        token.end_line = 1

        assert transformer.word([token]) == ("word", "str", "1")

    def test_bool(self):
        transformer = PrimitiveTransformer()
        token = Token("BOOL", "True")
        token.end_line = 1

        assert transformer.boolean([token]) == (True, "bool", "1")



class TestComposedTransformer:

    def test__init__(self):
        transformer = ComposedDatatypesTransformer()
        assert transformer.__visit_tokens__ == False

        transformer = ComposedDatatypesTransformer(True)
        assert transformer.__visit_tokens__ == True

        transformer = ComposedDatatypesTransformer(False)
        assert transformer.__visit_tokens__ == False

    def test__infer_most_general_type(self):
        transformer = ComposedDatatypesTransformer()
        assert transformer._infer_most_general_type(
            ["int",
            "float"]
        ) == "float"
        assert transformer._infer_most_general_type(["int", "int"]) == "int"
        assert transformer._infer_most_general_type(
            ["float",
            "int"]
        ) == "float"
        assert transformer._infer_most_general_type(
            ["float",
            "float"]
        ) == "float"
        assert transformer._infer_most_general_type(["bool", "bool"]) == "bool"
        assert transformer._infer_most_general_type(["str", "int"]) == "str"
        assert transformer._infer_most_general_type(["str", "float"]) == "str"
        assert transformer._infer_most_general_type(["str", "bool"]) == "str"
        assert transformer._infer_most_general_type(["str", "str"]) == "str"

        with pytest.raises(PQTypeError) as exception:
            transformer._infer_most_general_type(["bool", "int"])
        assert str(
            exception.value
        ) == "Bool cannot be used with other types. Found [\'bool\', \'int\']"

        with pytest.raises(PQTypeError) as exception:
            transformer._infer_most_general_type(["bool", "float"])
        assert str(
            exception.value
        ) == "Bool cannot be used with other types. Found [\'bool\', \'float\']"

    def test_array(self):
        transformer = ComposedDatatypesTransformer()
        token1 = (1, "int", "1")
        token2 = (2, "int", "1")

        assert transformer.array([token1,
            token2]) == ([1,
            2],
            "list(int)",
            "1")

        token1 = (1.0, "float", "1")
        token2 = (2.0, "float", "1")

        assert transformer.array([token1,
            token2]) == ([1.0,
            2.0],
            "list(float)",
            "1")

        token1 = (True, "bool", "1")
        token2 = (False, "bool", "1")

        assert transformer.array([token1,
            token2]) == ([True,
            False],
            "list(bool)",
            "1")

        token1 = ("word1", "str", "1")
        token2 = ("word2", "str", "1")

        string_list_type = "list(str)"

        assert transformer.array([token1,
            token2]) == (["word1",
            "word2"],
            string_list_type,
            "1")

        token1 = (1, "int", "1")
        token2 = (1.0, "float", "1")

        assert transformer.array([token1,
            token2]) == ([1.0,
            1.0],
            "list(float)",
            "1")

        token1 = (1, "int", "1")
        token2 = (True, "bool", "1")

        with pytest.raises(PQTypeError) as exception:
            transformer.array([token1, token2])
        assert str(
            exception.value
        ) == "Bool cannot be used with other types. Found [\'int\', \'bool\']"

        token1 = (1, "int", "1")
        token2 = ("word", "str", "1")

        assert transformer.array([token1,
            token2]) == (["1",
            "word"],
            string_list_type,
            "1")

        token1 = (1.0, "float", "1")
        token2 = ("word", "str", "1")

        assert transformer.array([token1,
            token2]) == (["1.0",
            "word"],
            string_list_type,
            "1")

        token1 = (True, "bool", "1")
        token2 = ("word", "str", "1")

        assert transformer.array([token1,
            token2]) == (["True",
            "word"],
            string_list_type,
            "1")

        token1 = (1, "int", "1")
        token2 = (2, "range", "1")

        with pytest.raises(PQTypeError) as exception:
            transformer.array([token1, token2])
        assert str(
            exception.value
        ) == "Array elements must be primitive types. Found [\'range\']"

    def test_range(self):
        transformer = ComposedDatatypesTransformer()
        token1 = (1, "int", "1")
        token2 = (2, "int", "1")

        assert transformer.range([token1,
            token2]) == (range(1,
            2),
            "range",
            "1")

        token1 = (1, "int", "1")
        token2 = (4, "int", "1")
        token3 = (20, "int", "1")

        assert transformer.range([token1,
            token2,
            token3]) == (range(1,
            20,
            4),
            "range",
            "1")

    def test_glob(self):
        transformer = ComposedDatatypesTransformer()
        token = Token("GLOB", "*")
        token.end_line = 1

        assert len(transformer.glob([token])[0]) > 0
        assert transformer.glob([token])[1] == "glob"
        assert transformer.glob([token])[2] == "1"

    def test_key(self):
        transformer = ComposedDatatypesTransformer()
        token = "key"
        assert transformer.key([token]) == "key"

    def test_value(self):
        transformer = ComposedDatatypesTransformer()
        token = ("value", "type", "1")
        assert transformer.value([token]) == ("value", "type", "1")
