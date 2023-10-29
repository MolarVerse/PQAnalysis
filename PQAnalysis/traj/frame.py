"""
A module containing the Frame class and its associated methods.

...

Classes
-------
Frame
    A class to represent a frame of a trajectory.
"""

import numpy as np

from PQAnalysis.selection.selection import Selection
from PQAnalysis.atom.molecule import Molecule
from PQAnalysis.pbc.cell import Cell


class Frame:
    '''
    A class to represent a frame of a trajectory.

    ...

    Attributes
    ----------
    n_atoms : int
        The number of atoms in the frame.
    xyz : np.array(n_atoms, 3)
        The xyz coordinates of the frame.
    atoms : np.array(n_atoms)
        The atoms in the frame.
    cell : Cell
        The cell of the frame.

    Methods
    -------
    __getindex__(index)
        Returns a new Frame with the given index.
    compute_com(group=None)
        Computes the center of mass of the frame.
    '''

    def __init__(self, n_atoms: int, xyz: np.array, atoms: np.array, cell: Cell = None):
        """
        Initializes the Frame with the given number of atoms, xyz coordinates, atoms, and cell.

        Parameters
        ----------
        n_atoms : int
            The number of atoms in the frame.
        xyz : np.array(n_atoms, 3)
            The xyz coordinates of the frame.
        atoms : np.array(n_atoms)
            The atoms in the frame.
        cell : Cell, optional
            The cell of the frame.
        """

        self.n_atoms = n_atoms
        self.xyz = np.array(xyz)
        self.atoms = np.array(atoms)
        self.cell = cell

    @property
    def PBC(self):
        """
        Returns True if the frame has a cell, False otherwise.

        Returns
        -------
        bool
            True if the frame has a cell, False otherwise.
        """
        if self.cell is not None:
            return True
        else:
            return False

    def __getitem__(self, index):
        """
        Makes the Frame indexable.

        Parameters
        ----------
        index : int or Selection
            The index of the new Frame.

        Raises
        ------
        ValueError
            If the selection is empty.

        Returns
        -------
        Frame
            The new Frame with the given index.
        """
        if isinstance(index, Selection):
            index = index.selection

        atoms = self.atoms[index]
        xyz = self.xyz[index]

        frame = Frame(len(atoms), xyz, atoms, cell=self.cell)

        if frame.n_atoms == 0:
            raise ValueError('Selection is empty.')

        return frame

    def compute_com(self, group=None):
        """
        Computes the center of mass of the frame.

        Divides the frame into groups of atoms and computes the center of mass of each group.
        If group is None, the group is all atoms.

        Parameters
        ----------
        group : int, optional
            The group to compute the center of mass for.

        Raises
        ------
        ValueError
            If the number of atoms in the selection is not a multiple of group.

        Returns
        -------
        Frame
            The new Frame with the center of mass of each group.
        """
        if group is None:
            group = self.n_atoms
        elif self.n_atoms % group != 0:
            raise ValueError(
                'Number of atoms in selection is not a multiple of group.')

        com = np.zeros((self.n_atoms // group, 3))
        molecule_names = np.zeros(self.n_atoms // group, dtype=object)

        for i in range(0, self.n_atoms, group):
            molecule = Molecule(self.atoms[i:i+group], self.xyz[i:i+group])

            com[i] = molecule.com(self.cell)

            molecule_names[i] = molecule.name

        return Frame(self.n_atoms // group, com, molecule_names, cell=self.cell)

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

        if not np.array_equal(self.xyz, other.xyz):
            return False

        if not np.array_equal(self.atoms, other.atoms):
            return False

        if self.cell != other.cell:
            return False

        return True
