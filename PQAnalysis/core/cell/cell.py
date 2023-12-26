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

from ...types import Np3x3NumberArray, Np2DNumberArray, Np1DNumberArray
from .standardProperties import _StandardPropertiesMixin


class Cell(_StandardPropertiesMixin):
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
    box_lengths : np.array
        The lengths of the box vectors.
    box_angles : np.array
        The angles between the box vectors.
    box_matrix : np.array
        The matrix containing the box vectors as columns.
    inverse_box_matrix : np.array
        The inverse of the box matrix.
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
        self._box_lengths = np.array([x, y, z])
        self._box_angles = np.array([alpha, beta, gamma])
        self._box_matrix = self.setup_box_matrix()

    def setup_box_matrix(self) -> Np3x3NumberArray:
        """
        Calculates the box matrix from the given parameters.

        Returns
        -------
        box matrix: Np3x3NumberArray
            The box matrix.
        """
        matrix = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

        alpha, beta, gamma = np.deg2rad(self.box_angles)
        cos_alpha, cos_beta, cos_gamma = np.cos([alpha, beta, gamma])
        sin_gamma = np.sin(gamma)
        x, y, z = self.box_lengths

        matrix[0][0] = x
        matrix[0][1] = y * cos_gamma
        matrix[0][2] = z * cos_beta
        matrix[1][1] = y * sin_gamma
        matrix[1][2] = z * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
        matrix[2][2] = z * np.sqrt(1 - cos_beta**2 -
                                   (cos_alpha - cos_beta * cos_gamma)**2 / sin_gamma**2)

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

        if self.alpha == 90 and self.beta == 90 and self.gamma == 90:
            pos = pos - self.box_lengths * np.rint(pos / self.box_lengths)
            return pos

        original_shape = np.shape(pos)

        if original_shape == (3,):
            pos = np.reshape(pos, (1, 3))

        fractional_pos = pos @ self.inverse_box_matrix.T

        fractional_pos -= np.round(fractional_pos)

        pos = fractional_pos @ self.box_matrix.T

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
        is_equal &= np.allclose(self.box_lengths, other.box_lengths)
        is_equal &= np.allclose(self.box_angles, other.box_angles)
        return is_equal
