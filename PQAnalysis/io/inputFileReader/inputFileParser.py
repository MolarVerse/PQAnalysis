from __future__ import annotations

from lark import Visitor, Transformer, Lark, Tree
from glob import glob
from pathlib import Path
from beartype.typing import Any, List, Tuple
from numbers import Integral, Real

from ...types import Range, Bool
from .. import BaseReader
from .formats import InputFileFormat


class InputFileParser(BaseReader):
    def __init__(self, filename: str, format: InputFileFormat = InputFileFormat.PQANALYSIS) -> None:
        super().__init__(filename)
        self.format = InputFileFormat(format)

    def parse(self) -> InputDictionary:
        if self.format == InputFileFormat.PQANALYSIS:
            grammar_file = Path(__file__).parent / "inputGrammar.lark"
        elif self.format == InputFileFormat.PIMD_QMCF or self.format == InputFileFormat.QMCFC:
            grammar_file = Path(__file__).parent / \
                "PIMD_QMCF_inputGrammar.lark"

        parser = Lark.open(grammar_file, rel_to=__file__,
                           propagate_positions=True)

        file = open(self.filename, "r")
        self.raw_input_file = file.read()
        self.tree = parser.parse(self.raw_input_file)

        self.transformed_tree = PrimitiveTransformer(
            visit_tokens=True).transform(self.tree)
        self.transformed_tree = ComposedDatatypesTransformer(
            visit_tokens=True).transform(self.transformed_tree)

        visitor = InputFileVisitor()
        self.input_dictionary = visitor.visit(self.transformed_tree)

        return self.input_dictionary


class InputDictionary:
    """
    Input file dictionary.

    This is a dictionary-like object that stores the input file keys and values.
    It is case-insensitive, so that the keys "key", "Key" and "KEY" are all
    equivalent. It also stores the type of the value and the line where the key
    was defined in the input file.
    """

    def __init__(self) -> None:
        """
        Initialize the dictionary.
        """
        self.dict = {}

    def __getitem__(self, key: str) -> Any:
        """
        Get the value of a key.

        Parameters
        ----------
        key : str
            The key to get the value of. It is case-insensitive.

        Returns
        -------
        Any
            The value of the key.

        Raises
        ------
        KeyError
            If the key is not defined in the input file.
        """
        key = key.lower()

        if key not in self.dict.keys():
            raise KeyError(
                f"Input file key \"{key}\" not defined in input file.")

        return self.dict[key]

    def __setitem__(self, key: str, value: Tuple[Any, str, str]) -> None:
        """
        Set the value of a key.

        Parameters
        ----------
        key : str
            The key to set the value of. It is case-insensitive.
        value : Tuple[Any, str, str]
            The value to set. It is a tuple containing the value, the type of
            the value and the line where the key was defined in the input file.

        Raises
        ------
        KeyError
            If the key is already defined in the input file.
        """
        key = key.lower()

        if key in self.dict.keys():
            raise KeyError(
                f"Input file key \"{key}\" defined multiple times in input file.")

        self.dict[key] = value

    def get_value(self, key: str) -> Any:
        """
        Get the value of a key.

        Parameters
        ----------
        key : str
            The key to get the value of. It is case-insensitive.

        Returns
        -------
        Any
            The value of the key.
        """
        return self.dict[key.lower()][0]

    def get_line(self, key: str) -> str:
        """
        Get the line where the key was defined in the input file.

        Parameters
        ----------
        key : str
            The key to get the line of. It is case-insensitive.

        Returns
        -------
        str
            The line where the key was defined in the input file.
        """
        return self.dict[key.lower()][2]

    def get_type(self, key: str) -> str:
        """
        Get the type of the value of a key.

        Parameters
        ----------
        key : str
            The key to get the type of. It is case-insensitive.

        Returns
        -------
        str
            The type of the value of the key.
        """
        return self.dict[key.lower()][1]

    def keys(self) -> List[str]:
        """
        Get the keys of the dictionary.

        Returns
        -------
        List[str]
            The keys of the dictionary in a list.
        """
        return list(self.dict.keys())


class PrimitiveTransformer(Transformer):
    """
    _summary_

    Parameters
    ----------
    Transformer : _type_
        _description_
    """

    def __init__(self, visit_tokens=False):
        super().__init__(visit_tokens)

    def float(self, items) -> Tuple[Real, str, str]:
        return float(items[0]), "float", str(items[0].end_line)

    def int(self, items) -> Tuple[Integral, str, str]:
        return int(items[0]), "int", str(items[0].end_line)

    def word(self, items) -> Tuple[str, str, str]:
        return str(items[0]), "str", str(items[0].end_line)

    def bool(self, items) -> Tuple[Bool, str, str]:
        return bool(items[0]), "bool", str(items[0].end_line)


class ComposedDatatypesTransformer(Transformer):
    def __init__(self, visit_tokens=False):
        super().__init__(visit_tokens)

    def array(self, items) -> Tuple[List[Any], str, str]:
        _list = [item[0] for item in items]
        return list(_list), f"list({items[0][1]})", str(items[0][2])

    def range(self, items) -> Tuple[Range, str, str]:
        return range(int(items[0][0]), int(items[2][0]), int(items[1][0])), "range", str(items[0].end_line)

    def glob(self, items) -> Tuple[List[str], str, str]:
        return glob(items[0].strip()), "glob", str(items[0].end_line)

    def key(self, items) -> str:
        return items[0]

    def value(self, items) -> Tuple[Any, str, str]:
        return items[0]


class InputFileVisitor(Visitor):
    def __init__(self):
        self.dict = InputDictionary()

    def assign(self, items):
        self.dict[str(items.children[0])] = items.children[1]
        return items

    def multiline_statement(self, items):
        _list = [item[0] for item in items.children[1:]]
        self.dict[str(items.children[0])] = (list(
            _list), f"list({items.children[1][1]})", f"{items.children[0].line}-{items.children[-1][2]}")
        return items

    def visit(self, tree: Tree) -> InputDictionary:
        super().visit(tree)
        return self.dict
