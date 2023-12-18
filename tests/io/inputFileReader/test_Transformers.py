import pytest

from lark import Token

from PQAnalysis.io.inputFileReader.inputFileParser import PrimitiveTransformer, ComposedDatatypesTransformer


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

        assert transformer.int([token]) == (1, "int", "1")

    def test_word(self):
        transformer = PrimitiveTransformer()
        token = Token("WORD", "word")
        token.end_line = 1

        assert transformer.word([token]) == ("word", "str", "1")

    def test_bool(self):
        transformer = PrimitiveTransformer()
        token = Token("BOOL", "True")
        token.end_line = 1

        assert transformer.bool([token]) == (True, "bool", "1")


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
            ["int", "float"]) == "float"
        assert transformer._infer_most_general_type(["int", "int"]) == "int"
        assert transformer._infer_most_general_type(
            ["float", "int"]) == "float"
        assert transformer._infer_most_general_type(
            ["float", "float"]) == "float"
        assert transformer._infer_most_general_type(["bool", "bool"]) == "bool"
        assert transformer._infer_most_general_type(["str", "int"]) == "str"
        assert transformer._infer_most_general_type(["str", "float"]) == "str"
        assert transformer._infer_most_general_type(["str", "bool"]) == "str"
        assert transformer._infer_most_general_type(["str", "str"]) == "str"

        with pytest.raises(TypeError) as exception:
            transformer._infer_most_general_type(["bool", "int"])
        assert str(
            exception.value) == "Bool cannot be used with other types. Found [\'bool\', \'int\']"

        with pytest.raises(TypeError) as exception:
            transformer._infer_most_general_type(["bool", "float"])
        assert str(
            exception.value) == "Bool cannot be used with other types. Found [\'bool\', \'float\']"
