from lark import Tree, Token

from .. import pytestmark

from PQAnalysis.io.input_file_reader.input_file_parser import InputFileVisitor, InputDictionary



class TestInputFileVisitor:

    def test__init__(self):
        visitor = InputFileVisitor()
        assert visitor.dict == InputDictionary()

    def test_assign(self):
        visitor = InputFileVisitor()
        tree = Tree("assign", ["key", ("value", "word", "1")])

        return_tree = visitor.assign(tree)

        assert return_tree == tree
        assert visitor.dict.dict == {"key": ("value", "word", "1")}

    def test_multiline_statement(self):
        visitor = InputFileVisitor()
        token = Token("END", "END")
        token.end_line = 4
        key_token = Token("WORD", "key")
        key_token.end_line = 1
        tree = Tree(
            "multiline_statement",
            [key_token,
            ("1",
            "float",
            "1"),
            ("4",
            "int",
            "4"),
            token]
        )

        return_tree = visitor.multiline_statement(tree)

        assert return_tree == tree
        assert visitor.dict.dict == {"key": ([1.0, 4.0], "list(float)", "1-4")}

    def test_visit(self):
        visitor = InputFileVisitor()
        tree = Tree(
            "input_file",
            [Tree("assign",
            ["key",
            ("value",
            "word",
            "1")])]
        )

        return_dict = visitor.visit(tree)

        assert return_dict == visitor.dict
        assert visitor.dict.dict == {"key": ("value", "word", "1")}
