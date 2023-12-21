from lark import Tree

from PQAnalysis.io.inputFileReader.inputFileParser import InputFileVisitor, InputDictionary


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
        tree = Tree("multiline_statement", [
                    "key", ("1", "float", "1"), ("4", "int", "4")])

        return_tree = visitor.multiline_statement(tree)

        assert return_tree == tree
        assert visitor.dict.dict == {"key": ([1.0, 4.0], "list(float)", "1-4")}

    def test_visit(self):
        visitor = InputFileVisitor()
        tree = Tree("input_file", [
                    Tree("assign", ["key", ("value", "word", "1")])])

        return_dict = visitor.visit(tree)

        assert return_dict == visitor.dict
        assert visitor.dict.dict == {"key": ("value", "word", "1")}
