"""
A module containing the Selection class and related functions/classes.

...

Classes
-------
Selection
    A class for representing a selection.
SelectionTransformer
    A class for transforming a Lark parse tree.
CrudeSelectionVisitor
    A super class for visiting a Lark parse tree.
SelectionVisitor
    A class for visiting a Lark parse tree and returning the indices of parsed selection.
    
Functions
---------
_selection
    An overloaded function for selecting atoms based on a selection compatible object.
_indices_by_atom_type_name
    Returns the indices of the atoms with the given atom type name.
_indices_by_atom
    Returns the indices of the given atom.
_indices_by_element_types
    Returns the indices of the given element type.
"""

from __future__ import annotations

import numpy as np

from lark import Visitor, Tree, Lark, Transformer, Token
from beartype.typing import List
from multimethod import overload

from .topology import Topology
from .. import __base_path__
from ..types import Np1DIntArray
from ..core import Atom, Atoms, Element

"""
A type hint for a selection compatible object.

A selection compatible object can be:
    - a string
    - an Atoms object (i.e. a list of Atom objects)
    - an Atom object
    - a numpy.ndarray with dtype=int
    - a list of strings
    - a Selection object
    - None
"""
SelectionCompatible = str | Atoms | Atom | Np1DIntArray | List[str] | 'Selection' | None


class Selection:
    def __init__(self, selection_object: str | Atoms | Atom | Np1DIntArray | List[str] | Selection | None):
        """
        A class for representing a selection.

        If the selection object is a Selection object, the selection is copied.

        Parameters
        ----------
        selection_object : str | Atoms | Atom | Np1DIntArray | List[str] | Selection | None
            The object to create the selection from.
        """
        if isinstance(selection_object, Selection):
            self.selection_object = selection_object.selection_object
        else:
            self.selection_object = selection_object

    def select(self, topology: Topology, use_full_atom_info: bool = False) -> Np1DIntArray:
        """
        Returns the indices of the atoms selected by the selection object.

        If the selection_object is None all atoms are selected.
        With the parameter use_full_atom_info the selection can be performed 
        based on the full atom information (i.e. atom type names and element type information).
        If use_full_atom_info is False, the selection is performed based on the element type information only.

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

        return _selection(self.selection_object, topology, use_full_atom_info)


@overload
def _selection(atoms: Atoms | Atom, topology: Topology, use_full_atom_info: bool) -> Np1DIntArray:
    """
    Overloaded function for selecting atoms based on a list of atoms or a single atom.

    Parameters
    ----------
    atoms : Atoms | Atom
        The atoms to get the indices of.
    topology : Topology
        The topology to get the indices from.
    use_full_atom_info : bool
        Whether to use the full atom information, by default False

    Returns
    -------
    Np1DIntArray
        The indices of the atoms selected by the selection object.
    """
    if isinstance(atoms, Atom):
        atoms = [atoms]

    indices = []

    for atom in atoms:
        if use_full_atom_info:
            indices.append(_indices_by_atom(atom, topology))
        else:
            indices.append(_indices_by_element_types(atom, topology))

    return np.sort(np.concatenate(indices))


@_selection.register
def _selection(atomtype_names: List[str], topology: Topology, *_) -> Np1DIntArray:
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


@_selection.register
def _selection(indices: Np1DIntArray, *_) -> Np1DIntArray:
    """
    Overloaded function for selecting atoms based on a list of indices.

    Parameters
    ----------
    indices : Np1DIntArray
        The indices to get the indices of.

    Returns
    -------
    Np1DIntArray
        The indices of the atoms selected by the selection object.
    """
    return indices


@_selection.register
def _selection(string: str, topology: Topology, use_full_atom_info: bool) -> Np1DIntArray:
    """
    Overloaded function for selecting atoms based on a string.

    This function implements a Lark parser for the selection grammar.
    Here the selection grammar is defined as follows:
        TODO: define selection grammar

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

    parser = Lark.open(grammar_path / grammar_file,
                       propagate_positions=True)

    tree = parser.parse(string)
    transformed_tree = SelectionTransformer(visit_tokens=True).transform(tree)
    visitor = SelectionVisitor(topology, use_full_atom_info)

    return np.sort(visitor.visit(transformed_tree))


class SelectionTransformer(Transformer):
    """
    A class for transforming a Lark parse tree.

    Parameters
    ----------
    Transformer : Transformer
        The type of the Transformer class to inherit from.
    """

    def __init__(self, visit_tokens=False):
        """
        Initializes the SelectionTransformer with the given parameters.

        Parameters
        ----------
        visit_tokens : bool, optional
            Whether to visit the tokens, by default False
        """
        self.__visit_tokens__ = visit_tokens
        super().__init__(visit_tokens=visit_tokens)

    def word(self, items: Token) -> str:
        """
        Returns the word of the given token.

        Parameters
        ----------
        items : Token
            The token to get the word of.

        Returns
        -------
        str
            The word of the given token.
        """
        return "".join(items)

    def letters(self, items: Token) -> str:
        """
        Returns the letters of the given token.

        Parameters
        ----------
        items : Token
            The token to get the letters of.

        Returns
        -------
        str
            The letters of the given token.
        """
        return "".join(items)

    def unsigned_int(self, items: Token) -> int:
        """
        Returns the unsigned integer of the given token.

        Parameters
        ----------
        items : Token
            The token to get the unsigned integer of.

        Returns
        -------
        int
            The unsigned integer of the given token.
        """
        return int("".join(items))

    def int(self, items: Token) -> int:
        """
        Returns the integer of the given token.

        Parameters
        ----------
        items : Token
            The token to get the integer of.

        Returns
        -------
        int
            The integer of the given token.
        """
        return int("".join(items))


class CrudeSelectionVisitor(Visitor):
    """
    A super class for visiting a Lark parse tree.

    It is used to define the visit methods for the SelectionVisitor class.

    Parameters
    ----------
    Visitor : Visitor
        The type of the Visitor class to inherit from.
    """

    def atom(self, items: Tree) -> Atom:
        """
        Returns the atom of the given items.

        Parameters
        ----------
        items : Tree
            The items to get the atom of.

        Returns
        -------
        Atom
            The atom of the given items.
        """
        return Atom(items.children[0], items.children[1])

    def element(self, items: Tree) -> Element:
        """
        Returns the element of the given items.

        Parameters
        ----------
        items : Tree
            The items to get the element of.

        Returns
        -------
        Element
            The element of the given items.
        """
        return Element(items.children[0])

    def atom_type(self, items: Tree) -> str:
        """
        Returns the atom type of the given items.

        Parameters
        ----------
        items : Tree
            The items to get the atom type of.

        Returns
        -------
        str
            The atom type of the given items.
        """
        return items.children[0]

    def index(self, items: Tree) -> int:
        """
        Returns the index of the given items.

        Parameters
        ----------
        items : Tree
            The items to get the index of.

        Returns
        -------
        int
            The index of the given items.
        """
        return items.children[0]

    def indices(self, items: Tree) -> Np1DIntArray:
        """
        Returns the indices of the given items.

        Parameters
        ----------
        items : Tree
            The items to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the given items.
        """
        if len(items.children) == 2:
            return np.arange(items.children[0], items.children[1]+1)
        elif len(items.children) == 3:
            return np.arange(items.children[0], items.children[2]+1, items.children[1])

    def residue(self, items: Tree) -> int | str:
        """
        Returns the residue of the given items.

        Parameters
        ----------
        items : Tree
            The items to get the residue of.

        Returns
        -------
        int | str
            The residue of the given items.
        """
        return items.children[0]

    def visit(self, tree: Tree) -> None:
        """
        Visits the given tree.

        Parameters
        ----------
        tree : Tree
            The tree to visit.
        """
        super().visit(tree)


class SelectionVisitor(CrudeSelectionVisitor):
    """
    A class for visiting a Lark parse tree and returning the indices of parsed selection.

    Parameters
    ----------
    Visitor : Visitor
        The type of the Visitor class to inherit from.
    """

    def __init__(self, topology: Topology, use_full_atom_info: bool = False):
        """
        Initializes the SelectionVisitor with the given parameters.

        Parameters
        ----------
        topology : Topology
            The topology to get the indices from.
        use_full_atom_info : bool, optional
            Whether to use the full atom information, by default False
        """
        self.topology = topology
        self.use_full_atom_info = use_full_atom_info
        self.selection = []
        self.super().__init__()

    def atomtype(self, items) -> Np1DIntArray:
        """
        Returns the indices of the atoms with the given atom type name.

        Parameters
        ----------
        items : List[str]
            The atom type name to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the atoms with the given atom type name.
        """
        atomtype = super().atomtype(items)

        for index in _indices_by_atom_type_name(atomtype, self.topology):
            self.selection.append(index)

        return np.array(self.selection)

    def atom(self, items) -> Np1DIntArray:
        """
        Returns the indices of the given atom.

        Parameters
        ----------
        items : 
            The atom to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the given atom.
        """

        atom = super().atom(items)

        if self.use_full_atom_info:
            for index in _indices_by_atom(atom, self.topology):
                self.selection.append(index)
        else:
            for index in _indices_by_element_types(atom, self.topology):
                self.selection.append(index)

        return np.array(self.selection)

    def element(self, items) -> Np1DIntArray:
        """
        Returns the indices of the given element type.

        Parameters
        ----------
        items : 
            The element type to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the given element type.
        """

        element = super().element(items)

        for index in _indices_by_element_types(element, self.topology):
            self.selection.append(index)

    def residue(self, items) -> Np1DIntArray:
        raise NotImplementedError("Residue selection is not implemented yet.")

    def index(self, items) -> Np1DIntArray:
        """
        Returns the indices of the given atom.

        Parameters
        ----------
        items : List[int]
            The index to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the given atom.
        """
        index = super().index(items)

        self.selection.append(index)

        return np.array(self.selection)

    def indices(self, items) -> Np1DIntArray:
        """
        Returns the indices of the given atoms.

        Parameters
        ----------
        items : List[Np1DIntArray]
            The indices to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the given atoms.
        """

        indices = super().indices(items)

        for index in indices:
            self.selection.append(index)

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
        [topology_name == name for topology_name in topology.atomtype_names])

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


def _indices_by_element_types(atom: Atom | Element, topology: Topology) -> Np1DIntArray:
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

    bool_indices = np.array(
        [topology_atom.element == element for topology_atom in topology.atoms])

    indices = np.argwhere(bool_indices).flatten()

    return indices
