"""
A module containing the Frame class and its associated methods.

...

Classes
-------
Frame
    A class to represent a frame of a trajectory.
"""

import numpy as np

from PQAnalysis.atomicUnits.molecule import Molecule
from PQAnalysis.pbc.cell import Cell
from PQAnalysis.atomicUnits.element import Elements
from PQAnalysis.utils.exceptions import ElementNotFoundError
from PQAnalysis.coordinates.coordinates import Coordinates


class BaseFrame:
    def __init__(self, coordinates: Coordinates):
        self.coordinates = coordinates


class Frame(BaseFrame):
    def __init__(self, coordinates: Coordinates, atoms: Elements):

        # just to make sure it is an Elements object
        self.atoms = Elements(atoms)
        # just to make sure it is a Coordinates object
        coordinates = Coordinates(coordinates)

        super().__init__(coordinates)

        if self.coordinates.n_atoms != self.atoms.n_atoms:
            raise ValueError(
                'coordinates and atoms must have the same length.')

    @property
    def n_atoms(self):
        try:
            return self.coordinates.n_atoms
        except AttributeError:
            return 0

    @property
    def cell(self):
        return self.coordinates.cell

    @property
    def PBC(self):
        return self.coordinates.PBC

    def __getitem__(self, index):
        frame = Frame(self.coordinates[index], self.atoms[index])

        if frame.n_atoms == 0:
            raise ValueError('Selection is empty.')

        return frame

    def compute_com_frame(self, group=None):
        if group is None:
            group = self.n_atoms

        elif self.n_atoms % group != 0:
            raise ValueError(
                'Number of atoms in selection is not a multiple of group.')

        molecules = []

        j = 0
        for i in range(0, self.n_atoms, group):
            molecule = Molecule(
                atoms=self.atoms[i:i+group], coordinates=self.coordinates[i:i+group])

            molecules.append(molecule)

            j += 1

        return MolecularFrame(molecules=molecules, cell=self.cell)

    def is_combinable(self, other: 'Frame') -> bool:
        """
        Checks if two Frames can be combined.

        Parameters
        ----------
        other : Frame
            The Frame to check if it can be combined with.

        Returns
        -------
        bool
            True if the Frames can be combined, False otherwise.
        """

        if not isinstance(other, Frame):
            return False

        if self.n_atoms != other.n_atoms:
            return False

        if self.atoms != other.atoms:
            return False

        return True

    def __eq__(self, other: 'Frame') -> bool:
        """
        Checks if two Frames are equal.

        Parameters
        ----------
        other : Frame
            The Frame to compare to.

        Returns
        -------
        bool
            True if the Frames are equal, False otherwise.
        """

        if not isinstance(other, Frame):
            return False

        if self.n_atoms != other.n_atoms:
            return False

        if self.atoms != other.atoms:
            return False

        if self.coordinates != other.coordinates:
            return False

        if self.cell != other.cell:
            return False

        return True


class MolecularFrame(BaseFrame):
    def __init__(self, molecules: list, coordinates: Coordinates = None, cell: Cell = None):

        self.molecules = molecules

        if coordinates is None:
            coordinates = self.center_of_masses

        super().__init__(coordinates)

        if self.coordinates.n_atoms != len(molecules):
            raise ValueError(
                'coordinates and atoms must have the same length.')

        if cell is not None:
            self.cell = cell
        else:
            self.cell = self.coordinates.cell

    @property
    def center_of_masses(self):
        coordinates = Coordinates()

        for molecule in self.molecules:
            coordinates.append(molecule.center_of_mass)

        return coordinates

    @property
    def n_molecules(self):
        return len(self.molecules)
