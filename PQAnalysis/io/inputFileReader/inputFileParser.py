from __future__ import annotations

from lark import Visitor, Transformer, Lark
from glob import glob

from .. import BaseReader
from .formats import InputFileFormat


class InputFileParser(BaseReader):
    def __init__(self, filename: str, format: InputFileFormat = InputFileFormat.PQANALYSIS):
        super().__init__(filename)
        self.format = InputFileFormat(format)

    def parse(self) -> InputDictionary:
        if self.format == InputFileFormat.PQANALYSIS:
            grammar_file = open("inputGrammar.lark", "r")
        elif self.format == InputFileFormat.PIMD_QMCF or self.format == InputFileFormat.QMCFC:
            grammar_file = open("PIMD_QMCF_inputGrammar.lark", "r")

        grammar = grammar_file.read()
        parser = Lark(grammar, propagate_positions=True)

        file = open(self.filename, "r")
        self.tree = parser.parse(file.read())

        self.transformed_tree = PrimitiveTransformer(
            visit_tokens=True).transform(self.tree)
        self.transformed_tree = ComposedDatatypesTransformer(
            visit_tokens=True).transform(self.transformed_tree)

        visitor = InputFileVisitor()
        self.input_dictionary = visitor.visit(self.transformed_tree)

        return self.input_dictionary


class InputDictionary:
    def __init__(self):
        self.dict = {}

    def __getitem__(self, key):
        return self.dict[key]

    def __setitem__(self, key, value):
        self.dict[key] = value

    def add(self, key, value):
        self.dict[key] = value

    def get_value(self, key):
        return self.dict[str(key)][0]

    def get_line(self, key):
        return self.dict[key][2]

    def get_type(self, key):
        return self.dict[key][1]


class PrimitiveTransformer(Transformer):
    def __init__(self, visit_tokens=False):
        super().__init__(visit_tokens)

    def float(self, items):
        return float(items[0]), "float", items[0].end_line

    def int(self, items):
        return int(items[0]), "int", items[0].end_line

    def string(self, items):
        return str(items[0]), "str", items[0].end_line

    def bool(self, items):
        return bool(items[0]), "bool", items[0].end_line


class ComposedDatatypesTransformer(Transformer):
    def __init__(self, visit_tokens=False):
        super().__init__(visit_tokens)

    def array(self, items):
        _list = [item[0] for item in items]
        return list(_list), f"list({items[0][1]})", items[0][2]

    def range(self, items):
        return range(int(items[0][0]), int(items[2][0]), int(items[1][0])), "range"

    def glob(self, items):
        return (glob(items[0].strip()), "glob")

    def key(self, items):
        return items[0]

    def value(self, items):
        return items[0]


class InputFileVisitor(Visitor):
    def __init__(self):
        self.dict = InputDictionary()

    def assign(self, items):
        self.dict[str(items.children[0])] = items.children[1]
        return items

    def multiline_statement(self, items):
        _list = [item[0] for item in items.children[1:]]
        self.dict[str(items.children[0])] = (list(_list), f"list({items.children[1][1]})", range(
            items.children[0].line, items.children[-1][2] + 1))
        return items

    def visit(self, tree):
        super().visit(tree)
        return self.dict
