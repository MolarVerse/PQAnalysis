import numpy as np

from lark import Visitor, Tree, Lark

from ... import __base_path__
from ..types import Np1DIntArray
from ..topology import Topology
from ..core.atom import Atom, is_same_element_type


class Selection:
    def __init__(self, string: str, topology: Topology):
        self.string = string
        self.topology = topology
        self.selection = []

    @property
    def selection(self) -> Np1DIntArray:
        return selection(self.string, self.topology)


def selection(string: str, topology: Topology) -> Np1DIntArray:
    grammar_file = "selection.lark"
    grammar_path = __base_path__ / "grammar"

    parser = Lark.open(grammar_path / grammar_file,
                       propagate_positions=True)

    tree = parser.parse(string)
    visitor = SelectionVisitor(topology)

    return visitor.visit(tree)


class SelectionVisitor(Visitor):
    def __init__(self, topology: Topology):
        self.topology = topology
        self.selection = []
        self.super().__init__()

    def atomtype(self, items) -> Np1DIntArray:
        for index in self.atomic_system._indices_by_atom_type_names([items[0]]):
            self.selection.append(index)

        return np.array(self.selection)

    def atom(self, items) -> Np1DIntArray:
        for index in self.atomic_system._indices_by_atom([items[0]]):
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


def _indices_by_atom_type_name(self, name: str, topology: Topology) -> Np1DIntArray:
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
        [atom.name == name for atom in topology.atoms])

    indices = np.argwhere(bool_array).flatten()

    return indices


def _indices_by_atom(self, atom: Atom, topology: Topology) -> Np1DIntArray:
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

    indices = np.argwhere(np.array(self.atoms) == atom).flatten()

    return indices


def _indices_by_element_types(self, element: Atom) -> Np1DIntArray:
    """
    Returns the indices of the given element type.

    Parameters
    ----------
    element: Atom
        The element type to get the indices of.

    Returns
    -------
    Np1DIntArray
        The indices of the given element type.
    """
    bool_indices = np.array(
        [is_same_element_type(atom, element) for atom in self.atoms])

    indices = np.argwhere(bool_indices).flatten()

    return indices
