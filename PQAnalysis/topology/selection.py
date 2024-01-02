"""
A module containing the Selection class and related functions/classes.

...

Classes
-------
Selection
    A class for representing a selection.
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

from lark import Visitor, Tree, Lark
from beartype.typing import List
from multimethod import overload

from .topology import Topology
from .. import __base_path__
from ..types import Np1DIntArray
from ..core import Atom, Atoms

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
    visitor = SelectionVisitor(topology, use_full_atom_info)

    return np.sort(visitor.visit(tree))


class SelectionVisitor(Visitor):
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
        for index in _indices_by_atom_type_name(items[0], self.topology):
            self.selection.append(index)

        return np.array(self.selection)

    def atom(self, items) -> Np1DIntArray:
        """
        Returns the indices of the given atom.

        Parameters
        ----------
        items : List[Atom]
            The atom to get the indices of.

        Returns
        -------
        Np1DIntArray
            The indices of the given atom.
        """
        if self.use_full_atom_info:
            for index in _indices_by_atom(items[0], self.topology):
                self.selection.append(index)
        else:
            for index in _indices_by_element_types(items[0], self.topology):
                self.selection.append(index)

        return np.array(self.selection)

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
        self.selection.append(items[0])

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
        if len(items) == 2:
            for index in range(items[0], items[1]):
                self.selection.append(index)
        elif len(items) == 3:
            for index in range(items[0], items[2], items[1]):
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


def _indices_by_element_types(atom: Atom, topology: Topology) -> Np1DIntArray:
    """
    Returns the indices of the given element type.

    Parameters
    ----------
    atom: Atom
        The atom type to get the indices of.
    topology : Topology
        The topology to get the indices from.

    Returns
    -------
    Np1DIntArray
        The indices of the given element type.
    """
    bool_indices = np.array(
        [topology_atom.element == atom.element for topology_atom in topology.atoms])

    indices = np.argwhere(bool_indices).flatten()

    return indices
