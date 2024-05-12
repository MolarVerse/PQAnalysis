"""
A module containing the Selection class and related functions/classes.
"""
# library imports
import logging
import numpy as np

# 3rd party object imports
from lark import Visitor, Tree, Lark, Transformer, Token
from beartype.typing import List, TypeVar, Any

# local imports
from PQAnalysis import __base_path__, __package_name__
from PQAnalysis.types import Np1DIntArray
from PQAnalysis.core import Atom, Atoms, Element, Elements
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.exceptions import PQValueError
from PQAnalysis.type_checking import runtime_type_checking

from .topology import Topology

#: | A type variable for the Selection class.
#: | It can be used to specify the type of the selection object.
#: | The following types are supported:
#: |     - str: the given string is parsed and the atoms selected by the selection are selected
#: |     - Atom: the given atom is selected
#: |     - Element: all atoms with the given element type are selected
#: |     - Atoms: all atoms in the given list are selected
#: |     - Elements: all atoms with the given element types are selected
#: |     - Np1DIntArray: the atoms with the given indices are selected
#: |     - List[str]: all atoms with the given atom type names are selected
#: |     - Selection: the given selection is copied
#: |     - None: all atoms are selected
SelectionCompatible = TypeVar(
    'SelectionCompatible',
    str,
    'PQAnalysis.core.Atom',
    List['PQAnalysis.core.Atom'],
    'PQAnalysis.core.Element',
    List['PQAnalysis.core.Element'],
    Np1DIntArray,
    List[str],
    'PQAnalysis.topology.Selection',
    None
)



class Selection:

    """
    A class for representing a selection.

    If the selection object is a Selection object, the selection is copied.

    There are several ways to create a selection:
        - None: all atoms are selected
        - Atom: the given atom is selected
        - Element: all atoms with the given element type are selected
        - Atoms: all atoms in the given list are selected
        - Elements: all atoms with the given element types are selected
        - Np1DIntArray: the atoms with the given indices are selected
        - List[str]: all atoms with the given atom type names are selected
        - str: the given string is parsed and the atoms selected by the 
        selection are selected

            This string will be parsed based on a Lark grammar. Which is 
            defined as follows:

                - simple word containing only letters and numbers: the
                atom type with the given name is selected
                - <integer>: the atom with the given index is selected
                - <integer1>..<integer2>: the indices from integer1 to integer2 
                are selected
                - <integer1>-<integer2>: the indices from integer1 to integer2
                are selected
                - <integer1>..<integer2>..<integer3>: the indices from integer1 to
                integer3 with a step size of integer3 are selected
                - atom(<atomtype>, <atomic_number>): the atom with the given atom 
                type and atomic number is selected
                - atom(<atomtype>, <element_symbol>): the atom with the given atom
                type and element symbol is selected
                - elem(<atomic_number>): all atoms with the given element type
                are selected
                - elem(<element_symbol>): all atoms with the given element type
                are selected
                - `*`: all atoms are selected (same as 'all'), useful if only
                few atoms should be excluded

            All of the above statements can be combined with the following operators:

                - ',': the union of the two statements is selected, meaning that
                the atoms selected by the first statement and the atoms selected
                by the second statement are selected
                - '&': the intersection of the two statements is selected, meaning 
                that only the atoms selected by both statements are selected
                - '|': the set difference of the two statements is selected, 
                meaning that only the atoms selected by the first statement are
                selected, which are not selected by the second statement

            The operators are evaluated in the following order: '|' -> '&' -> ','
            This means that the ',' operator has the lowest precedence and the '|' 
            operator has the highest precedence and therefore binds the strongest.

            Additionally, parentheses can be used to group statements and change
            the order of evaluation.


    Notes
    -----
    The atom counting always starts with 0!
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(self, selection_object: SelectionCompatible = None):
        """
        Parameters
        ----------
        selection_object : SelectionCompatible, optional
            The object to create the selection from.
        """
        if isinstance(selection_object, Selection):
            self.selection_object = selection_object.selection_object
        else:
            self.selection_object = selection_object

    @runtime_type_checking
    def select(
        self,
        topology: Topology,
        use_full_atom_info: bool = False
    ) -> Np1DIntArray:
        """
        Returns the indices of the atoms selected by the selection object.

        If the selection_object is None all atoms are selected. With the 
        parameter use_full_atom_info the selection can be performed based
        on the full atom information (i.e. atom type names and element type
        information). If use_full_atom_info is False, the selection is 
        performed based on the element type information only.

        Parameters
        ----------
        topology : Topology
            The topology to get the indices from.
        use_full_atom_info : bool, optional
            Whether to use the full atom information, by default False

        Returns
        -------
        Np1DIntArray
            The indices of the atoms selected by the selection object.
        """
        if self.selection_object is None:
            return np.arange(topology.n_atoms)

        return np.unique(
            _selection(self.selection_object,
            topology,
            use_full_atom_info)
        )

    def __str__(self) -> str:
        """
        Returns the string representation of the selection object.

        Returns
        -------
        str
            The string representation of the selection object.
        """
        return str(self.selection_object)

    def __repr__(self) -> str:
        """
        Returns the string representation of the selection object.

        Returns
        -------
        str
            The string representation of the selection object.
        """
        return str(self.selection_object)



def _selection(
    atoms: SelectionCompatible,
    topology: Topology,
    use_full_atom_info: bool
) -> Np1DIntArray:
    """
    Overloaded function for selecting atoms based
    on a list of atoms/elements or a single atom/element.

    Parameters
    ----------
    atoms : SelectionCompatible
        The selection compatible object to get the indices of.
    topology : Topology
        The topology to get the indices from.
    use_full_atom_info : bool
        Whether to use the full atom information, by default False

    Returns
    -------
    Np1DIntArray
        The indices of the atoms selected by the selection object.
    """

    if isinstance(atoms, (Atom, Element)):
        return _selection_of_atoms(atoms, topology, use_full_atom_info)

    if (isinstance(atoms,
        List) and len(atoms) > 0 and isinstance(atoms[0],
        (Atom,
        Element))):
        return _selection_of_atoms(atoms[0], topology, use_full_atom_info)

    if isinstance(atoms, List):
        return _selection_of_atomtypes(atoms, topology)

    if isinstance(atoms, str):
        return _selection_of_string(atoms, topology, use_full_atom_info)

    return atoms



def _selection_of_atoms(
    atoms: Atoms | Atom | Element | Elements,
    topology: Topology,
    use_full_atom_info: bool
) -> Np1DIntArray:
    """
    Returns the indices of the atoms selected by the selection object.

    Parameters
    ----------
    atoms : Atoms | Atom | Element | Elements
        The atoms/elements to get the indices of.
    topology : Topology
        The topology to get the indices from.
    use_full_atom_info : bool
        Whether to use the full atom information, by default False

    Returns
    -------
    Np1DIntArray
        The indices of the atoms selected by the selection object.

    Raises
    ------
    PQValueError
        If the use_full_atom_info parameter is True and the atoms parameter is an Element object.
    """
    if isinstance(atoms, (Atom, Element)):
        atoms = [atoms]

    if isinstance(atoms[0], Element):
        if use_full_atom_info:
            Selection.logger.error(
                "The use_full_atom_info parameter is not supported for Element objects.",
                exception=PQValueError
            )

        atoms = [Atom(element.symbol, element.symbol) for element in atoms]

    indices = []

    for atom in atoms:
        if use_full_atom_info:
            indices.append(_indices_by_atom(atom, topology))
        else:
            indices.append(_indices_by_element_types(atom, topology))

    return np.sort(np.concatenate(indices))



def _selection_of_atomtypes(
    atomtype_names: List[str],
    topology: Topology
) -> Np1DIntArray:
    """
    Overloaded function for selecting atoms based on a list of atom type names.

    Parameters
    ----------
    atomtype_names : List[str]
        The atom type names to get the indices of.
    topology : Topology
        The topology to get the indices from.

    Returns
    -------
    Np1DIntArray
        The indices of the atoms selected by the selection object.
    """
    indices = []

    for atomtype_name in atomtype_names:
        indices.append(_indices_by_atom_type_name(atomtype_name, topology))

    return np.sort(np.concatenate(indices))



def _selection_of_string(
    string: str,
    topology: Topology,
    use_full_atom_info: bool
) -> Np1DIntArray:
    """
    Overloaded function for selecting atoms based on a string.

    This function implements a Lark parser for the selection grammar.

    Parameters
    ----------
    string : str
        The string to get the indices of.
    topology : Topology
        The topology to get the indices from.
    use_full_atom_info : bool
        Whether to use the full atom information, by default False

    Returns
    -------
    Np1DIntArray
        The indices of the atoms selected by the selection object.
    """
    grammar_file = "selection.lark"
    grammar_path = __base_path__ / "grammar"

    parser = Lark.open(grammar_path / grammar_file, propagate_positions=True)

    tree = parser.parse(string)
    transformed_tree = SelectionTransformer(
        topology=topology,
        visit_tokens=True,
        use_full_atom_info=use_full_atom_info
    ).transform(tree)
    visitor = SelectionVisitor()

    return np.sort(visitor.visit(transformed_tree))



class SelectionTransformer(Transformer):

    """
    A class for transforming a Lark parse tree.

    Parameters
    ----------
    Transformer : Transformer
        The type of the Transformer class to inherit from.
    """

    def __init__(
        self,
        topology=Topology(),
        visit_tokens=False,
        use_full_atom_info=False
    ):
        """
        Parameters
        ----------
        visit_tokens : bool, optional
            Whether to visit the tokens, by default False
        topology : Topology
            The topology to get the indices from.
        use_full_atom_info : bool, optional
            Whether to use the full atom information, by default False
        """
        self.__visit_tokens__ = visit_tokens
        super().__init__(visit_tokens=visit_tokens)

        self.topology = topology
        self.use_full_atom_info = use_full_atom_info

    def word(self, items: List[Token]) -> str:
        """
        Returns the word of the given token.

        Parameters
        ----------
        items : List[Token]
            The token to get the word of.

        Returns
        -------
        str
            The word of the given token.
        """
        return "".join(items)

    def letters(self, items: List[Token]) -> str:
        """
        Returns the letters of the given token.

        Parameters
        ----------
        items : List[Token]
            The token to get the letters of.

        Returns
        -------
        str
            The letters of the given token.
        """
        return "".join(items)

    def unsigned_integer(self, items: List[Token]) -> int:
        """
        Returns the unsigned integer of the given token.

        Parameters
        ----------
        items : List[Token]
            The token to get the unsigned integer of.

        Returns
        -------
        int
            The unsigned integer of the given token.
        """
        return int(items[0])

    def integer(self, items: List[Token]) -> int:
        """
        Returns the integer of the given token.

        Parameters
        ----------
        items : List[Token]
            The token to get the integer of.

        Returns
        -------
        int
            The integer of the given token.
        """
        return int(items[0])

    def UNSIGNED_INT(self, items: Token | List[Token]) -> int:  # pylint: disable=invalid-name
        """
        Returns the unsigned integer of the given token.

        Parameters
        ----------
        items : List[Token] | Token
            The token to get the unsigned integer of.

        Returns
        -------
        int
            The unsigned integer of the given token.
        """
        return int(items[0])

    def INT(self, items: List[Token] | Token) -> int:  # pylint: disable=invalid-name
        """
        Returns the integer of the given token.

        Parameters
        ----------
        items : List[Token] | Token
            The token to get the integer of.

        Returns
        -------
        int
            The integer of the given token.
        """
        return int(items[0])

    def atomtype(self, items: Any) -> Np1DIntArray:
        """
        Returns the indices of the atoms with the given atom type name.

        Parameters
        ----------
        items : Any
            The atom type name to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the atoms with the given atom type name.
        """
        atomtype = items[0]

        return _indices_by_atom_type_name(atomtype, self.topology)

    def atom(self, items: Any) -> Np1DIntArray:
        """
        Returns the indices of the given atom.

        Parameters
        ----------
        items : Any
            The atom to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the given atom.
        """

        atom = Atom(items[0], items[1])

        if self.use_full_atom_info:
            return _indices_by_atom(atom, self.topology)

        return _indices_by_element_types(atom, self.topology)

    def element(self, items: Any) -> Np1DIntArray:
        """
        Returns the indices of the given element type.

        Parameters
        ----------
        items : Any
            The element type to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the given element type.
        """

        element = Element(items[0])

        return _indices_by_element_types(element, self.topology)

    def residue(self, items: Any) -> Np1DIntArray:
        """
        Returns the indices of the given residue name.

        Parameters
        ----------
        items : Any
            The residue name to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the given residue name.
        """
        return np.array(
            self.topology.get_atom_indices_from_residue_names(items[0])
        )

    def residue_number(self, items: Any) -> Np1DIntArray:
        """
        Returns the indices of the given residue number.

        Parameters
        ----------
        items : Any
            The residue number to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the given residue number.
        """

        return np.array(
            self.topology.get_atom_indices_from_residue_numbers(items[0])
        )

    def index(self, items: List[int]) -> Np1DIntArray:
        """
        Returns the given index as an array.

        Parameters
        ----------
        items: List[int]
            The index to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the index Token.
        """
        index = items[0]

        return np.array([index])

    def indices(self, items: List[int]) -> Np1DIntArray:
        """
        Returns a range of indices based on the given indices.

        if two indices are given, a range from the first index to
        the second index is returned.
        if three indices are given, a range from the first index to
        the third index with a step size of the second index is returned.

        Parameters
        ----------
        items : List[int]
            The indices to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the indices Token.
            
        Raises
        ------
        PQValueError
            If the indices Token has more than 3 items. Should never be reached.
        """

        if len(items) == 2:
            return np.arange(items[0], items[1] + 1)

        if len(items) == 3:
            return np.arange(items[0], items[2] + 1, items[1])

        # should never be reached
        Selection.logger.error(
            "The indices Token must have 2 or 3 items.",
            exception=PQValueError
        )

    def all(self, _) -> Np1DIntArray:
        """
        Returns all indices.

        Returns
        -------
        Np1DIntArray
            The indices of the all Token.
        """

        return np.arange(self.topology.n_atoms)

    def without_statement(self, items: Any) -> Np1DIntArray:
        """
        Returns the indices of the the first item
        without the indices of all other items.

        Parameters
        ----------
        items : List[Np1DIntArray]
            The indices to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the given without statement Token.
        """

        difference = items[0]

        for item in items[1:]:
            difference = np.setdiff1d(difference, item)

        return difference

    def and_statement(self, items: Any) -> Np1DIntArray:
        """
        Returns the indices of the and_statement Token.

        Parameters
        ----------
        items : List[Np1DIntArray]
            The indices to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the and_statement Token.
        """

        intersection = items[0]

        for item in items:
            intersection = np.intersect1d(intersection, item)

        return intersection

    def or_statement(self, items: Any) -> Np1DIntArray:
        """
        Returns the indices of the or_statement Token.

        Parameters
        ----------
        items : List[Np1DIntArray]
            The indices to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the or_statement Token.
        """

        union = items[0]

        for item in items:
            union = np.union1d(union, item)

        return union



class SelectionVisitor(Visitor):

    """
    A class for visiting a Lark parse tree and returning the indices of parsed selection.

    Parameters
    ----------
    Visitor : Visitor
        The type of the Visitor class to inherit from.
    """

    def __init__(self):
        """
        Parameters
        ----------
        topology : Topology
            The topology to get the indices from.
        use_full_atom_info : bool, optional
            Whether to use the full atom information, by default False
        """
        self.selection = []
        super().__init__()

    def expression(self, tree: Tree) -> Np1DIntArray:
        """
        Visits the given tree and returns the indices of the parsed selection.

        Parameters
        ----------
        tree : Tree
            The tree to visit.

        Returns
        -------
        Np1DIntArray
            The indices of the parsed selection.
        """
        self.selection = tree.children[0]

        return np.array(self.selection)

    def visit(self, tree: Tree) -> Np1DIntArray:
        """
        Visits the given tree and returns the indices of the parsed selection.

        Parameters
        ----------
        tree : Tree
            The tree to visit.

        Returns
        -------
        Np1DIntArray
            The indices of the parsed selection.
        """
        super().visit(tree)

        return np.array(self.selection)



def _indices_by_atom_type_name(name: str, topology: Topology) -> Np1DIntArray:
    """
    Returns the indices of the atoms with the given atom type name.

    Parameters
    ----------
    name : str
        The name of the atom type to get the indices of.
    topology : Topology
        The topology to get the indices from.

    Returns
    -------
    Np1DIntArray
        The indices of the atoms with the given atom type name.
    """

    bool_array = np.array(
        [topology_name == name for topology_name in topology.atomtype_names]
    )

    indices = np.argwhere(bool_array).flatten()

    return indices



def _indices_by_atom(atom: Atom, topology: Topology) -> Np1DIntArray:
    """
    Returns the indices of the given atom.

    Parameters
    ----------
    atom : Atom
        The atom to get the indices of.
    topology : Topology
        The topology to get the indices from.

    Returns
    -------
    Np1DIntArray
        The indices of the given atom.
    """

    indices = np.argwhere(np.array(topology.atoms) == atom).flatten()

    return indices



def _indices_by_element_types(
    atom: Atom | Element,
    topology: Topology
) -> Np1DIntArray:
    """
    Returns the indices of the given element type.

    It can be use for both Atom and Element objects.

    Parameters
    ----------
    atom: Atom | Element
        The atom or element to get the indices of.
    topology : Topology
        The topology to get the indices from.

    Returns
    -------
    Np1DIntArray
        The indices of the given element type.
    """

    if isinstance(atom, Atom):
        element = atom.element
    else:
        element = atom

    bool_indices = np.array(
        [topology_atom.element == element for topology_atom in topology.atoms]
    )

    indices = np.argwhere(bool_indices).flatten()

    return indices
