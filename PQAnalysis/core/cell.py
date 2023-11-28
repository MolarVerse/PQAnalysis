"""
A module containing the Cell class.

...

Classes
-------
Cell
    A class for storing unit cell parameters.
"""

from __future__ import annotations

import numpy as np
import sys

from beartype.typing import Any
from numbers import Real

from ..types import Np3x3NumberArray, Np2DNumberArray, Np1DNumberArray


class Cell:
    '''
    Class for storing unit cell parameters.

    ...

    Attributes
    ----------

    x : Real, optional
        The length of the first box vector. Default is sys.float_info.max.
    y : Real, optional
        The length of the second box vector. Default is sys.float_info.max.
    z : Real, optional
        The length of the third box vector. Default is sys.float_info.max.
    alpha : Real, optional
        The angle between the second and third box vector. Default is 90.
    beta : Real, optional
        The angle between the first and third box vector. Default is 90.
    gamma : Real, optional
        The angle between the first and second box vector. Default is 90.
    box_matrix : np.array
        The matrix containing the box vectors as columns.
    '''

    def __init__(self,
                 x: Real = sys.float_info.max,
                 y: Real = sys.float_info.max,
                 z: Real = sys.float_info.max,
                 alpha: Real = 90,
                 beta: Real = 90,
                 gamma: Real = 90
                 ) -> None:
        """
        Initializes the Cell with the given parameters.

        If no angles are given, the cell is assumed to be orthorhombic.
        The box matrix is calculated from the given parameters.

        Parameters
        ----------
        x : Real
            The length of the first box vector.
        y : Real
            The length of the second box vector.
        z : Real
            The length of the third box vector.
        alpha : Real, optional
            The angle between the second and third box vector.
        beta : Real, optional
            The angle between the first and third box vector.
        gamma : Real, optional
            The angle between the first and second box vector.
        """
        self.x = x
        self.y = y
        self.z = z
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.box_matrix = self.setup_box_matrix()

    def setup_box_matrix(self) -> Np3x3NumberArray:
        """
        Calculates the box matrix from the given parameters.

        Returns
        -------
        box matrix: Np3x3NumberArray
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
    def bounding_edges(self) -> Np2DNumberArray:
        """
        calculates the coordinates of the eight corners of the unit cell.

        Returns
        -------
        edges: Np2DNumberArray of shape (8, 3)
            The coordinates of the eight corners of the unit cell.
        """
        edges = np.zeros((8, 3))
        for i, x in enumerate([-0.5, 0.5]):
            for j, y in enumerate([-0.5, 0.5]):
                for k, z in enumerate([-0.5, 0.5]):
                    edges[i*4+j*2+k, :] = self.box_matrix @ np.array([x, y, z])

        return edges

    @property
    def volume(self) -> Real:
        """
        Returns the volume of the unit cell.

        Returns
        -------
        Real
            The volume of the unit cell.
        """
        return np.linalg.det(self.box_matrix)

    @property
    def box_lengths(self) -> Np1DNumberArray:
        """
        Returns the lengths of the box vectors.

        Returns
        -------
        box_lengths: Np1DNumberArray of shape (3,)
            The lengths of the box vectors.
        """
        return np.array([self.x, self.y, self.z])

    @property
    def box_angles(self) -> Np1DNumberArray:
        """
        Returns the angles between the box vectors.

        Returns
        -------
        box_angles: Np1DNumberArray of shape (3,)
            The lengths of the box vectors.
        """
        return np.array([self.alpha, self.beta, self.gamma])

    def image(self, pos: Np2DNumberArray | Np1DNumberArray) -> Np2DNumberArray | Np1DNumberArray:
        """
        Returns the image of the given position in the unit cell.

        Parameters
        ----------
        pos : Np2DNumberArray, Np1DNumberArray
            The position to get the image of.

        Returns
        -------
        imaged_positions: Np2DNumberArray, Np1DNumberArray
            The image of the position(s) in the unit cell.
        """

        original_shape = np.shape(pos)

        if original_shape == (3,):
            pos = np.reshape(pos, (1, 3))

        fractional_pos = np.array(
            [np.linalg.inv(self.box_matrix) @ pos_i for pos_i in pos])

        fractional_pos -= np.round(fractional_pos)

        pos = np.array(
            [self.box_matrix @ fractional_pos_i for fractional_pos_i in fractional_pos])

        return np.reshape(pos, original_shape)

    def __eq__(self, other: Any) -> bool:
        """
        Checks if the Cell is equal to another Cell.

        Parameters
        ----------
        other : Cell
            The Cell to compare with.

        Returns
        -------
        bool
            True if the Cells are equal, False otherwise.
        """

        if not isinstance(other, Cell):
            return False

        is_equal = True
        is_equal &= np.allclose(self.x, other.x)
        is_equal &= np.allclose(self.y, other.y)
        is_equal &= np.allclose(self.z, other.z)
        is_equal &= np.allclose(self.alpha, other.alpha)
        is_equal &= np.allclose(self.beta, other.beta)
        is_equal &= np.allclose(self.gamma, other.gamma)
        return is_equal
