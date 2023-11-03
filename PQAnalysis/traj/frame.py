"""
A module containing the Frame class and its associated methods.

...

Classes
-------
Frame
    A class to represent a frame of a trajectory.
"""

import numpy as np

from PQAnalysis.traj.selection import Selection
from PQAnalysis.atom.molecule import Molecule
from PQAnalysis.pbc.cell import Cell
from PQAnalysis.atom.element import Element
from PQAnalysis.utils.exceptions import ElementNotFoundError


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
    '''

    def __init__(self, xyz: np.array, atoms: np.array, cell: Cell = None, atom_types_map: dict[str: Element] = None):
        """
        Initializes the Frame with the given number of atoms, xyz coordinates, atoms, and cell.

        Parameters
        ----------
        n_atoms : int
            The number of atoms in the frame.
        xyz : np.array(n_atoms, 3)
            The xyz coordinates of the frame.
        atoms : np.array(n_atoms)
            The atoms in the frame. Can be either an Iterable of strings
            with the atomic symbols or atom type names or an Iterable of 
            Element objects (recommended). For all post processing involving
            atomic masses and atomic numbers (e.g. computing the center of mass)
            it is required to use atomic symbols or Element objects otherwise a
            dicitonary with types {str, Element} has to be given.
        cell : Cell, optional
            The cell of the frame.
        atom_types_map : dict[str: Element], optional
            A dictionary mapping the atom type names to Element objects.

        Raises
        ------
        ValueError
            If the xyz coordinates and atoms are not of the same length.
        """

        self._xyz = xyz
        self._atoms = atoms

        self._atom_types_map = atom_types_map

        self.cell = cell

        if len(self.xyz) != len(self.atoms):
            raise ValueError('xyz and atoms must have the same length.')

    @property
    def atom_types_names(self):
        """
        Returns the atom type names in the frame.

        Returns
        -------
        np.array(n_atoms)
            The atom type names in the frame.
        """
        return self._atom_type_names

    @atom_types_names.setter
    def atom_types_names(self, atom_type_names):
        self._atom_type_names = np.array(atom_type_names)

    @property
    def atom_types_map(self):
        """
        Returns the atom types map in the frame.

        Returns
        -------
        dict[str: Element]
            The atom types map in the frame.
        """
        return self._atom_types_map

    @atom_types_map.setter
    def atom_types_map(self, atom_types_map):
        if atom_types_map is None:
            self._atom_types_map = None
            return

        if not all(isinstance(atom_type, str) for atom_type in list(atom_types_map.keys())):
            raise TypeError(
                'atom_types_map must be a dictionary with str keys and Element values.')

        if not all(isinstance(atom_type, Element) for atom_type in list(atom_types_map.values())):
            raise TypeError(
                'atom_types_map must be a dictionary with str keys and Element values.')

        self._atom_types_map = atom_types_map

        try:
            self.atoms = np.array([atom_types_map[atom_type_name]
                                  for atom_type_name in self.atom_type_names])
        except ElementNotFoundError:
            raise ElementNotFoundError(
                'atom_types_map does not contain all atom type names.')

    @property
    def n_atoms(self):
        """
        Returns the number of atoms in the frame.

        Returns
        -------
        int
            The number of atoms in the frame.
        """
        return len(self.atoms)

    @property
    def xyz(self):
        """
        Returns the xyz coordinates of the frame.

        Returns
        -------
        np.array(n_atoms, 3)
            The xyz coordinates of the frame.
        """
        return self._xyz

    @xyz.setter
    def xyz(self, xyz):

        xyz = np.array(xyz)

        if len(np.shape(xyz)) != 2 or np.shape(xyz)[1] != 3:
            raise ValueError(
                'xyz must be a iterable with following shape - (n_atoms, 3).')

        self._xyz = xyz

    @property
    def atoms(self):
        """
        Returns the atoms in the frame.

        Returns
        -------
        np.array(n_atoms)
            The atoms in the frame.
        """
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):

        atoms = np.array(atoms)

        if len(np.shape(atoms)) != 1:
            raise ValueError(
                'atoms must be a iterable with following shape - (n_atoms,).')

        if all(isinstance(atom, Element) for atom in atoms):
            self.atom_type_names = np.array(
                [atom.name for atom in atoms])
        elif all(isinstance(atom, str) for atom in atoms):
            self.atom_type_names = np.array(atoms)
        else:
            raise TypeError(
                'atoms must be either an Iterable of Element objects or an Iterable of strings.')

        self._atoms = atoms

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

        if isinstance(index, int):
            atoms = np.array([self.atoms[index]])
            xyz = np.array([self.xyz[index]])
        else:
            atoms = self.atoms[index]
            xyz = self.xyz[index]

        frame = Frame(xyz, atoms, cell=self.cell)

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

        j = 0
        for i in range(0, self.n_atoms, group):
            molecule = Molecule(self.atoms[i:i+group], self.xyz[i:i+group])

            com[j] = molecule.com(self.cell)

            molecule_names[j] = molecule.name

            j += 1

        return Frame(com, molecule_names, cell=self.cell)

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

        if not np.array_equal(self.atoms, other.atoms):
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

        if not np.array_equal(self.xyz, other.xyz):
            return False

        if not np.array_equal(self.atoms, other.atoms):
            return False

        if self.cell != other.cell:
            return False

        return True

    def __setup_atoms__(self):

        if all(isinstance(atom, Element) for atom in self.atoms):
            self.atom_type_names = np.array(
                [atom.name for atom in self.atoms])
            return

        self.atom_type_names = self.atoms
        self.atoms = None

        if self.atom_types_map is None:
            try:
                self.atoms = np.array([Element(atom) for atom in self.atoms])
                self.atom_type_names = np.array(
                    [atom.name for atom in self.atoms])

            except ElementNotFoundError:
                pass

            return

        if not all(isinstance(atom_type, str) for atom_type in list(self.atom_types_map.keys())):
            raise TypeError(
                'atom_types_map must be a dictionary with str keys and Element values.')

        if not all(isinstance(atom_type, Element) for atom_type in list(self.atom_types_map.values())):
            raise TypeError(
                'atom_types_map must be a dictionary with str keys and Element values.')
