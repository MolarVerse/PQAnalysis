"""
A module containing the Cell class.

...

Classes
-------
Cell
    A class for storing unit cell parameters.
"""

import numpy as np


class Cell:
    '''
    Class for storing unit cell parameters.

    ...

    Attributes
    ----------

    x : float
        The length of the first box vector.
    y : float
        The length of the second box vector.
    z : float
        The length of the third box vector.
    alpha : float
        The angle between the second and third box vector.
    beta : float
        The angle between the first and third box vector.
    gamma : float
        The angle between the first and second box vector.
    box_matrix : np.array
        The matrix containing the box vectors as columns.
    '''

    def __init__(self, x: float, y: float, z: float, alpha: float = 90, beta: float = 90, gamma: float = 90):
        """
        Initializes the Cell with the given parameters.

        If no angles are given, the cell is assumed to be orthorhombic.
        The box matrix is calculated from the given parameters.

        Parameters
        ----------
        x : float
            The length of the first box vector.
        y : float
            The length of the second box vector.
        z : float
            The length of the third box vector.
        alpha : float, optional
            The angle between the second and third box vector.
        beta : float, optional
            The angle between the first and third box vector.
        gamma : float, optional
            The angle between the first and second box vector.
        """
        self.x = x
        self.y = y
        self.z = z
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.box_matrix = self.setup_box_matrix()

    def setup_box_matrix(self) -> np.array:
        """
        Calculates the box matrix from the given parameters.

        Returns
        -------
        np.array(3, 3)
            The box matrix.
        """
        matrix = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

        alpha = np.deg2rad(self.alpha)
        beta = np.deg2rad(self.beta)
        gamma = np.deg2rad(self.gamma)

        matrix[0][0] = self.x
        matrix[0][1] = self.y * np.cos(gamma)
        matrix[0][2] = self.z * np.cos(beta)
        matrix[1][1] = self.y * np.sin(gamma)
        matrix[1][2] = self.z * \
            (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
        matrix[2][2] = self.z * np.sqrt(1 - np.cos(beta)**2 - (
            np.cos(alpha) - np.cos(beta) * np.cos(gamma))**2 / np.sin(gamma)**2)

        return matrix

    @property
    def bounding_edges(self) -> np.array:
        """
        calculates the coordinates of the eight corners of the unit cell.

        Returns
        -------
        np.array(8, 3)
            The coordinates of the eight corners of the unit cell.
        """
        edges = np.zeros((8, 3))
        for i, x in enumerate([-0.5, 0.5]):
            for j, y in enumerate([-0.5, 0.5]):
                for k, z in enumerate([-0.5, 0.5]):
                    edges[i*4+j*2+k, :] = self.box_matrix @ np.array([x, y, z])

        return edges

    @property
    def volume(self) -> float:
        """
        Returns the volume of the unit cell.

        Returns
        -------
        float
            The volume of the unit cell.
        """
        return np.linalg.det(self.box_matrix)

    @property
    def box_lengths(self) -> np.array:
        """
        Returns the lengths of the box vectors.

        Returns
        -------
        np.array(3)
            The lengths of the box vectors.
        """
        return np.array([self.x, self.y, self.z])

    @property
    def box_angles(self) -> np.array:
        """
        Returns the angles between the box vectors.

        Returns
        -------
        np.array(3)
            The lengths of the box vectors.
        """
        return np.array([self.alpha, self.beta, self.gamma])

    def image(self, pos: np.array) -> np.array:
        '''
        Images a position vector in the unit cell.

        This method works for both single position vectors and arrays of position vectors.

        Parameters
        ----------
        pos : np.array
            The position vector to image.

        Returns
        -------
        np.array
            The imaged position vector.
        '''
        if np.shape(pos) == (3,):
            pos = np.reshape(pos, (1, 3))

        fractional_pos = [np.linalg.inv(self.box_matrix) @ i for i in pos]

        fractional_pos -= np.round(fractional_pos)

        pos = [self.box_matrix @ i for i in fractional_pos]

        if np.shape(pos) == (1, 3):
            pos = np.reshape(pos, (3,))

        return pos

    def __eq__(self, __value: 'Cell') -> bool:
        """
        Checks if the Cell is equal to another Cell.

        Parameters
        ----------
        __value : Cell
            The Cell to compare with.

        Returns
        -------
        bool
            True if the Cells are equal, False otherwise.
        """

        if not isinstance(__value, Cell):
            return False

        is_equal = True
        is_equal &= self.x == __value.x
        is_equal &= self.y == __value.y
        is_equal &= self.z == __value.z
        is_equal &= self.alpha == __value.alpha
        is_equal &= self.beta == __value.beta
        is_equal &= self.gamma == __value.gamma
        return is_equal
