import numpy as np

from typing import List, Union
from collections.abc import Iterable

from PQAnalysis.atom.element import Element
from PQAnalysis.pbc.cell import Cell


class Molecule:
    def __init__(self, atoms: List[Union[str, Element]], xyz=None, name: str = None):

        # setting up the atoms
        if not isinstance(atoms, Iterable):
            raise TypeError('atoms must be iterable.')
        elif all(isinstance(atom, Element) for atom in atoms):
            atoms = atoms
        elif all(isinstance(atom, str) for atom in atoms):
            atoms = [Element(atom) for atom in atoms]
        else:
            raise TypeError('atoms must be either a list of Element or str.')

        # setting up name
        if not isinstance(name, str) and name is not None:
            raise TypeError('name must be a str.')
        elif name is None:
            name = ''.join(atoms)

        self.atoms = atoms
        self.xyz = xyz
        self.name = name

    @property
    def atom_masses(self) -> List[float]:
        return [atom.mass for atom in self.atoms]

    @property
    def mass(self) -> float:
        return sum(self.atom_masses)

    def com(self, cell: Cell = None) -> np.array:
        if self.xyz is None:
            raise ValueError('xyz must be provided when computing com.')

        if cell is not None:
            xyz_imaged = cell.image(self.xyz - self.xyz[0]) + self.xyz[0]
        else:
            xyz_imaged = self.xyz

        com = sum(np.array([self.atom_masses]).T *
                  xyz_imaged) / sum(self.atom_masses)

        if cell is not None:
            com = cell.image(com)

        return com
