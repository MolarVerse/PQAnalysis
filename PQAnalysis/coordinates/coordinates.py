import numpy as np

from PQAnalysis.pbc.cell import Cell
from PQAnalysis.traj.selection import Selection


def image(coordinates: 'Coordinates') -> 'Coordinates':
    '''
    Images a coordinates object in the unit cell.

    This method works for both single position vectors and arrays of position vectors.

    Parameters
    ----------
    pos : Coordinates
        coordinates to image.

    Returns
    -------
    Coordinates
        The imaged coordinates.
    '''

    if coordinates.cell is None:
        return coordinates

    fractional_xyz = [np.linalg.inv(
        coordinates.cell.box_matrix) @ i for i in coordinates.xyz]

    fractional_xyz -= np.round(fractional_xyz)

    xyz = [coordinates.cell.box_matrix @ i for i in fractional_xyz]

    return Coordinates(xyz, coordinates.cell)


class Coordinates():
    def __init__(self, xyz: np.array = None, cell: Cell = None):
        """
        Parameters
        ----------
        xyz : np.array
            The coordinates of the atoms in the frame.
        cell : Cell, optional
            The cell of the frame.
        """

        if xyz is None:
            xyz = np.zeros((0, 3))

        if len(np.shape(xyz)) == 1:
            xyz = np.array([xyz])

        self.xyz = xyz
        self.cell = cell

    @property
    def n_atoms(self):
        """
        Returns
        -------
        int
            The number of atoms in the frame.
        """
        return len(self.xyz)

    @property
    def PBC(self):
        """
        Returns
        -------
        bool
            Whether the frame has periodic boundary conditions.
        """
        return self.cell is not None

    def __add__(self, other):

        return Coordinates(self.xyz + other.xyz, self.cell)

    def __sub__(self, other):

        return Coordinates(self.xyz - other.xyz, self.cell)

    def __mul__(self, other):

        return Coordinates(self.xyz * other, self.cell)

    def __truediv__(self, other):

        return Coordinates(self.xyz / other, self.cell)

    def __getitem__(self, index):
        """
        Makes the Coordinates indexable.

        Parameters
        ----------
        index : int or Selection
            The index of the new Coordinates.

        Raises
        ------
        ValueError
            If the selection is empty.

        Returns
        -------
        Coordinates
            The new Coordinates with the given index.
        """
        if isinstance(index, Selection):
            index = index.selection

        return Coordinates(self.xyz[index], cell=self.cell)

    def append(self, other):
        """
        Appends a Coordinates object to the end of the current Coordinates object.

        Parameters
        ----------
        other : Coordinates
            The Coordinates object to append.

        Raises
        ------
        ValueError
            If the cells of the two Coordinates objects are not the same and if the other Coordinates object any cell defined.
        """

        if self.cell != other.cell and other.cell is not None:
            raise ValueError(
                'The cells of the two Coordinates objects must be the same.')

        self.xyz = np.append(self.xyz, other.xyz, axis=0)

    def __eq__(self, other):
        """
        Checks if two Coordinates objects are equal.

        Parameters
        ----------
        other : Coordinates
            The other Coordinates object to check if it is equal to.

        Returns
        -------
        bool
            Whether the two Coordinates objects are equal.
        """

        if self.cell != other.cell:
            return False

        return np.allclose(self.xyz, other.xyz)

    def copy(self):
        """
        Returns a copy of the Coordinates object.

        Returns
        -------
        Coordinates
            A copy of the Coordinates object.
        """
        return Coordinates(self.xyz.copy(), self.cell)
