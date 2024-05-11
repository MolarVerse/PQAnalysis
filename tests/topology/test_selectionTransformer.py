from lark import Token

from PQAnalysis.topology.selection import SelectionTransformer

from . import pytestmark



class TestSelectionTransformer:

    def test__init__(self):
        transformer = SelectionTransformer()
        assert transformer.__visit_tokens__ == False

        transformer = SelectionTransformer(visit_tokens=True)
        assert transformer.__visit_tokens__ == True

        transformer = SelectionTransformer(visit_tokens=False)
        assert transformer.__visit_tokens__ == False

    def test_word(self):
        transformer = SelectionTransformer()
        token = [Token("WORD", "word"), Token("WORD", "word")]
        assert transformer.word(token) == "wordword"

    def test_letters(self):
        transformer = SelectionTransformer()
        token = [Token("LETTERs", "letters"), Token("LETTERs", "letters")]
        assert transformer.letters(token) == "lettersletters"

    def test_integer(self):
        transformer = SelectionTransformer()
        token = [Token("INTEGER", "-12")]
        assert transformer.integer(token) == -12

    def test_unsigned_integer(self):
        transformer = SelectionTransformer()
        token = [Token("UNSIGNED_INTEGER", "12")]
        assert transformer.unsigned_integer(token) == 12

    def test_INT(self):
        transformer = SelectionTransformer()
        token = [Token("INT", "-12")]
        assert transformer.INT(token) == -12

    def test_UNSIGNED_INT(self):
        transformer = SelectionTransformer()
        token = [Token("UNSIGNED_INT", "12")]
        assert transformer.UNSIGNED_INT(token) == 12
