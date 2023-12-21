"""
A module containing the input file parser.

...

Classes
-------
InputFileParser
    Class to parse input files.
InputDictionary
    Class to store the input file keys and values.
PrimitiveTransformer
    Transformer for primitive datatypes.
ComposedDatatypesTransformer
    Transformer for composed datatypes.
InputFileVisitor
    Visitor for input files.
"""

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
    """
    Class to parse input files.

    Parameters
    ----------
    BaseReader : BaseReader
        BaseReader class from PQAnalysis.io
    """

    def __init__(self, filename: str, format: InputFileFormat | str = InputFileFormat.PQANALYSIS) -> None:
        """
        Initialize the parser.

        Parameters
        ----------
        filename : str
            The name of the input file.
        format : InputFileFormat | str, optional
            The format of the input file, by default InputFileFormat.PQANALYSIS
        """
        super().__init__(filename)
        self.format = InputFileFormat(format)

    def parse(self) -> InputDictionary:
        """
        Parse the input file.

        It uses the lark parser to parse the input file. It then uses the
        PrimitiveTransformer and ComposedDatatypesTransformer to transform the
        tree. Finally, it uses the InputFileVisitor to visit the tree and
        return the dictionary.

        Returns
        -------
        InputDictionary: InputDictionary
            The parsed input file dictionary.
        """
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

    def __eq__(self, __value: object) -> bool:
        """
        Compare two InputDictionary objects.

        Parameters
        ----------
        __value : object
            The object to compare to.

        Returns
        -------
        bool
            True if the objects are equal, False otherwise.
        """
        if not isinstance(__value, InputDictionary):
            return False

        return self.dict == __value.dict


class PrimitiveTransformer(Transformer):
    """
    Transformer for primitive datatypes.

    Parameters
    ----------
    Transformer : Transformer
        Transformer class from lark.
    """

    def __init__(self, visit_tokens=False):
        """
        initialize the transformer

        Parameters
        ----------
        visit_tokens : bool, optional
            boolean to visit tokens, by default False
        """
        self.__visit_tokens__ = visit_tokens
        super().__init__(self.__visit_tokens__)

    def float(self, items) -> Tuple[Real, str, str]:
        """
        Method to transform float values

        A "float" token is transformed into a float value, the string "float", and the line where the token was defined.

        Parameters
        ----------
        items: List[Token]
            items containing the float value

        Returns
        -------
        Tuple[Real, str, str]
            tuple containing the float value, the string "float", and the line where the token was defined.
        """
        return float(items[0]), "float", str(items[0].end_line)

    def int(self, items) -> Tuple[Integral, str, str]:
        """
        Method to transform int values

        Parameters
        ----------
        items: List[Token]
            items containing the int value  

        Returns
        -------
        Tuple[Integral, str, str]
            tuple containing the int value, the string "int", and the line where the token was defined.
        """
        return int(items[0]), "int", str(items[0].end_line)

    def word(self, items) -> Tuple[str, str, str]:
        """
        Method to transform word values

        Parameters
        ----------
        items : List[Token]
            items containing the word value

        Returns
        -------
        Tuple[str, str, str]
            tuple containing the word value, the string "str", and the line where the token was defined.
        """
        return str(items[0]), "str", str(items[0].end_line)

    def bool(self, items) -> Tuple[Bool, str, str]:
        """
        Method to transform bool values

        Parameters
        ----------
        items : List[Token]
            items containing the bool value

        Returns
        -------
        Tuple[Bool, str, str]
            tuple containing the bool value, the string "bool", and the line where the token was defined.
        """
        return bool(items[0]), "bool", str(items[0].end_line)


class ComposedDatatypesTransformer(Transformer):
    """
    Transformer for composed datatypes.

    Parameters
    ----------
    Transformer : Transformer
        Transformer class from lark.
    """

    primitive_types = ["float", "int", "str", "bool"]

    def __init__(self, visit_tokens=False):
        """
        initialize the transformer

        Parameters
        ----------
        visit_tokens : bool, optional
            boolean to visit tokens, by default False
        """
        self.__visit_tokens__ = visit_tokens
        super().__init__(self.__visit_tokens__)

    def _infer_most_general_type(self, types: List[str]) -> str:
        """
        Method to infer the most general type of a list of items

        Parameters
        ----------
        items : List[Any]
            list of items

        Returns
        -------
        str
            most general type of the list of items

        Raises
        ------
        TypeError
            if the list of items contains a bool and another type
        """
        if "str" in types:
            return "str"

        if "bool" in types and not all([item == "bool" for item in types]):
            raise TypeError(
                f"Bool cannot be used with other types. Found {types}")
        elif "bool" in types:
            return "bool"
        elif "float" in types:
            return "float"
        elif "int" in types:
            return "int"

    def array(self, items) -> Tuple[List[Any], str, str]:
        """
        Method to transform array values

        It transforms a list of items into a list of values, the string "{most_general_type}", and the line where the token was defined.
        The most general type is inferred from the list of items. If the list contains a bool and another type, a TypeError is raised.
        The list of items must contain only primitive types.

        Parameters
        ----------
        items : List[Any]
            items containing the array value

        Returns
        -------
        Tuple[List[Any], str, str]
            tuple containing the array value, the string "{most_general_type}", and the line where the token was defined.

        Raises
        ------
        TypeError
            if the list of items contains a non-primitive type
        """
        not_primitive_types = [
            item[1] for item in items if item[1] not in self.primitive_types]

        if len(not_primitive_types) > 0:
            raise TypeError(
                f"Array elements must be primitive types. Found {not_primitive_types}")

        most_general_type = self._infer_most_general_type(
            [item[1] for item in items])

        if most_general_type == "str":
            items = [(str(item[0]), item[1], item[2]) for item in items]
        elif most_general_type == "bool":
            items = [(bool(item[0]), item[1], item[2]) for item in items]
        elif most_general_type == "float":
            items = [(float(item[0]), item[1], item[2]) for item in items]
        elif most_general_type == "int":
            items = [(int(item[0]), item[1], item[2]) for item in items]

        _list = [item[0] for item in items]
        return list(_list), f"list({most_general_type})", str(items[0][2])

    def range(self, items) -> Tuple[Range, str, str]:
        """
        Method to transform range values

        Parameters
        ----------
        items : List[Any]
            items containing the range value

        Returns
        -------
        Tuple[Range, str, str]
            tuple containing the range value, the string "range", and the line where the token was defined.
        """
        if len(items) == 2:
            return_range = range(items[0][0], items[1][0])
        else:
            return_range = range(items[0][0], items[2][0], items[1][0])

        return return_range, "range", str(items[0][2])

    def glob(self, items) -> Tuple[List[str], str, str]:
        """
        Method to transform glob values

        Parameters
        ----------
        items : List[Any]
            items containing the glob value

        Returns
        -------
        Tuple[List[str], str, str]
            tuple containing the glob value, the string "glob", and the line where the token was defined.
        """
        return glob("".join(items).strip()), "glob", str(items[0].end_line)

    def key(self, items) -> str:
        """
        Method to transform key values

        Parameters
        ----------
        items : List[Token]
            items containing the key value

        Returns
        -------
        str
            key value
        """
        return items[0]

    def value(self, items) -> Tuple[Any, str, str]:
        """
        Method to transform value values

        Parameters
        ----------
        items
            items containing the value value

        Returns
        -------
        Tuple[Any, str, str]
            tuple containing the value value, the string "value", and the line where the token was defined.
        """
        return items[0]


class InputFileVisitor(Visitor):
    """
    Visitor for input files.

    Parameters
    ----------
    Visitor : 
        Visitor class from lark.
    """

    def __init__(self):
        """
        Initialize the visitor.

        It initializes the dictionary and the composed datatypes transformer.
        """
        self.dict = InputDictionary()
        self.composedDatatypeTransformer = ComposedDatatypesTransformer()

    def assign(self, items: Tree) -> Tree:
        """
        Parse an assign statement. The assign statement is of the form:

        key = value; 

        where key is a string, and the value is a tuple containing the value, the type of the value and the line where the key was defined.

        Parameters
        ----------
        items : Tree
            The assign statement.

        Returns
        -------
        Tree
            The assign statement.
        """
        self.dict[str(items.children[0])] = items.children[1]

        return items

    def multiline_statement(self, items: Tree) -> Tree:
        """
        Parse a multiline statement. The multiline statement is of the form:

        key
        value1
        value2
        ...
        END

        where key is a string, and the values are a tuple containing the value, the type of the value and the line where the key was defined.

        Parameters
        ----------
        items : Tree
            The multiline statement.

        Returns
        -------
        Tree
            The multiline statement.
        """
        array = self.composedDatatypeTransformer.array(
            [item for item in items.children[1:-1]])

        self.dict[str(items.children[0])
                  ] = array[0], array[1], f"{items.children[0].end_line}-{items.children[-1].end_line}"

        return items

    def visit(self, tree: Tree) -> InputDictionary:
        """
        Visit the tree and return the dictionary.

        Parameters
        ----------
        tree : Tree
            The tree to visit.

        Returns
        -------
        InputDictionary
            The parsed input file dictionary.
        """
        super().visit(tree)

        return self.dict
