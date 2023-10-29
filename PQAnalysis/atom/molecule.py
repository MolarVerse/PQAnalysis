"""
A module containing the Molecule class.

...

Classes
-------
Molecule
    A class representing a molecule
"""

import numpy as np

from collections.abc import Iterable
from typing import List, Union

from PQAnalysis.atom.element import Element
from PQAnalysis.pbc.cell import Cell


class Molecule:
    """
    A class representing a molecule

    ...

    Attributes
    ----------
    atoms : list of Element
        The list of atoms in the molecule
    xyz : np.array
        The xyz coordinates of the molecule
    name : str
        The name of the molecule

    Methods
    -------
    atom_masses
        Returns the masses of the atoms in the molecule.
    mass
        Returns the total mass of the molecule.
    com(cell)
        Returns the center of mass of the molecule.
    """

    def __init__(self, atoms: Union[List[str], List[Element]], xyz: np.array = None, name: str = None):
        """
        Parameters
        ----------
        atoms : list of str or list of Element
            The list of atoms in the molecule. Can be either a list of Element or str.
        xyz : np.array, optional
            The xyz coordinates of the molecule. If not provided, the com method will not work.
        name : str, optional
            The name of the molecule. If not provided, the name will be the concatenation of the atoms.

        Raises
        ------
        TypeError
            If the given atoms is not a list of Element or str.
        TypeError
            If the given name is not a str.
        """
        # setting up the atoms
        if not isinstance(atoms, Iterable) or isinstance(atoms, str):
            raise TypeError('atoms must be iterable or list of str.')
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
            name = self.__atoms_to_name__(atoms)

        self.atoms = atoms
        self.xyz = xyz
        self.name = name

    @property
    def atom_masses(self) -> List[float]:
        """
        Returns the masses of the atoms in the molecule.
        """
        return [atom.mass for atom in self.atoms]

    @property
    def mass(self) -> float:
        """
        Returns the total mass of the molecule.
        """
        return sum(self.atom_masses)

    def com(self, cell: Cell = None) -> np.array:
        """
        Returns the center of mass of the molecule.

        Parameters
        ----------
        cell : Cell, optional
            The cell to use for imaging the molecule. If not provided, no imaging will be done.

        Raises
        ------
        ValueError
            If xyz is not provided.
        """
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

    def __atoms_to_name__(self, atoms: Union[List[Element], List[str]]) -> str:
        """
        Returns the name of the molecule from the given atoms.

        The name is the concatenation of the names of the atoms.

        Parameters
        ----------
        atoms : list of Element or list of str
            The list of atoms in the molecule. Can be either a list of Element or str.
        """
        if all(isinstance(atom, str) for atom in atoms):
            name = ''.join(atoms)
        else:
            name = ''.join([atom.name for atom in atoms])

        return name
