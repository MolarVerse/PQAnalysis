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

from PQAnalysis.atomicUnits.element import Elements
from PQAnalysis.coordinates.coordinates import Coordinates, image
from PQAnalysis.pbc.cell import Cell


class Molecule:
    """
    A class representing a molecule.

    ...

    Attributes
    ----------
    coordinates : Coordinates
        the coordinates of the atoms in the molecule.
    atoms : Elements
        the atoms in the molecule as elements
    name : str
        the name of the molecule
    """

    def __init__(self, coordinates: Coordinates = None, atoms: Elements = None, name: str = None):
        """
        Initializes a Molecule object.

        Parameters
        ----------
        coordinates : Coordinates, optional
            the coordinates of the atoms in the molecule, by default None       
        atoms : Elements, optional
            the atoms in the molecule as elements, by default None
        name : str, optional
            the name of the molecule, by default None

        Raises
        ------
        TypeError
            If atoms is not an Elements object.
        TypeError
            If coordinates is not a Coordinates object.
        """
        if not isinstance(atoms, Elements) and atoms is not None:
            raise TypeError('atoms must be an Elements object.')

        if not isinstance(coordinates, Coordinates) and coordinates is not None:
            raise TypeError('coordinates must be a Coordinates object.')

        self.coordinates = coordinates
        self.atoms = atoms
        self.name = name

    @property
    def name(self) -> str:
        """
        If name is None, the name of the molecule is the concatenation of the names of the atoms 
        otherwise it is the name provided by the user.

        Returns
        -------
        str
            The name of the molecule.
        """
        if self._name is None:
            return ''.join(self.atoms.names)
        else:
            return self._name

    @name.setter
    def name(self, name: str):
        """
        Parameters
        ----------
        name : str
            The name of the molecule.

        Raises
        ------
        TypeError
            If name is not a str.
        """
        if name is not None and not isinstance(name, str):
            raise TypeError('name must be a str.')

        self._name = name

    @property
    def atom_masses(self) -> np.array:
        """
        Returns
        -------
        np.array
            The masses of the atoms in the molecule.
        """
        return self.atoms.masses

    @property
    def mass(self) -> float:
        """
        Returns the total mass of the molecule.
        """
        return sum(self.atoms.masses)

    @property
    def cell(self) -> Cell:
        """
        Returns
        -------
        Cell
            The cell of the molecule.
        """
        return self.coordinates.cell

    @property
    def center_of_mass(self) -> Coordinates:
        """
        Computes the center of mass of the molecule.

        Returns
        -------
        Coordinates
            The center of mass of the molecule in the same cell as the molecule as a Coordinates object.

        Raises
        ------
        ValueError
            If coordinates is None.
        """

        if self.coordinates is None:
            raise ValueError(
                'coordinates must be provided when computing com.')

        xyz_imaged = image(self.coordinates -
                           self.coordinates[0]) + self.coordinates[0]

        com = Coordinates(
            sum(np.array([self.atom_masses]).T * xyz_imaged.xyz) / self.mass, self.cell)

        com = image(com)

        return com
