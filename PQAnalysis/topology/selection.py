from __future__ import annotations

import numpy as np

from lark import Visitor, Tree, Lark
from beartype.typing import List
from multimethod import overload

from .. import __base_path__
from ..types import Np1DIntArray
from ..core.atom import Atom, Atoms, is_same_element_type
from .topology import Topology

SelectionCompatible = str | Atoms | Atom | Np1DIntArray | List[str] | 'Selection' | None


class Selection:
    def __init__(self, selection_object: str | Atoms | Atom | Np1DIntArray | List[str] | Selection | None):

        if isinstance(selection_object, Selection):
            self.selection_object = selection_object.selection_object
        else:
            self.selection_object = selection_object

    def select(self, topology: Topology, use_full_atom_info: bool = False) -> Np1DIntArray:

        if self.selection_object is None:
            return np.arange(topology.n_atoms)

        return _selection(self.selection_object, topology, use_full_atom_info)


@overload
def _selection(atoms: Atoms | Atom, topology: Topology, use_full_atom_info: bool) -> Np1DIntArray:

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

    indices = []

    for atomtype_name in atomtype_names:
        indices.append(_indices_by_atom_type_name(atomtype_name, topology))

    return np.sort(np.concatenate(indices))


@_selection.register
def _selection(indices: Np1DIntArray, *_) -> Np1DIntArray:
    return indices


@_selection.register
def _selection(string: str, topology: Topology, use_full_atom_info: bool) -> Np1DIntArray:
    grammar_file = "selection.lark"
    grammar_path = __base_path__ / "grammar"

    parser = Lark.open(grammar_path / grammar_file,
                       propagate_positions=True)

    tree = parser.parse(string)
    visitor = SelectionVisitor(topology, use_full_atom_info)

    return np.sort(visitor.visit(tree))


class SelectionVisitor(Visitor):
    def __init__(self, topology: Topology, use_full_atom_info: bool = False):
        self.topology = topology
        self.use_full_atom_info = use_full_atom_info
        self.selection = []
        self.super().__init__()

    def atomtype(self, items) -> Np1DIntArray:
        for index in _indices_by_atom_type_name(items[0], self.topology):
            self.selection.append(index)

        return np.array(self.selection)

    def atom(self, items) -> Np1DIntArray:
        if self.use_full_atom_info:
            for index in _indices_by_atom(items[0], self.topology):
                self.selection.append(index)
        else:
            for index in _indices_by_element_types(items[0], self.topology):
                self.selection.append(index)

        return np.array(self.selection)

    def index(self, items) -> Np1DIntArray:
        self.selection.append(items[0])

        return np.array(self.selection)

    def indices(self, items) -> Np1DIntArray:
        if len(items) == 2:
            for index in range(items[0], items[1]):
                self.selection.append(index)
        elif len(items) == 3:
            for index in range(items[0], items[2], items[1]):
                self.selection.append(index)

        return np.array(self.selection)

    def visit(self, tree: Tree) -> Np1DIntArray:
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


def _indices_by_element_types(element: Atom, topology: Topology) -> Np1DIntArray:
    """
    Returns the indices of the given element type.

    Parameters
    ----------
    element: Atom
        The element type to get the indices of.
    topology : Topology
        The topology to get the indices from.

    Returns
    -------
    Np1DIntArray
        The indices of the given element type.
    """
    bool_indices = np.array(
        [is_same_element_type(atom, element) for atom in topology.atoms])

    indices = np.argwhere(bool_indices).flatten()

    return indices
