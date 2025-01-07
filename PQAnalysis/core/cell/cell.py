"""
A module containing the Cell class.
"""

import sys
import warnings

from numbers import Real

import numpy as np

from beartype.typing import Any, NewType, Annotated
from beartype.vale import Is

from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.types import Np3x3NumberArray, Np2DNumberArray, NpnDNumberArray, PositiveReal, Bool
from PQAnalysis.utils.math import allclose_vectorized

from ._standard_properties import _StandardPropertiesMixin



class Cell(_StandardPropertiesMixin):

    """
    Class for storing unit cell parameters.
    """

    @runtime_type_checking
    def __init__(
        self,
        x: Real = sys.float_info.max,
        y: Real = sys.float_info.max,
        z: Real = sys.float_info.max,
        alpha: Real = 90,
        beta: Real = 90,
        gamma: Real = 90
    ) -> None:
        """
        A cell object can be initialized with the following parameters:

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

        Notes
        -----
        A vacuum cell can be created by calling Cell(), which is equivalent to 
        Cell(x=sys.float_info.max, y=sys.float_info.max, z=sys.float_info.max,
        alpha=90, beta=90, gamma=90).
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
        sin_beta = np.sin(beta)
        x, y, z = self.box_lengths

        matrix[0][0] = x
        matrix[0][1] = y * cos_gamma
        matrix[0][2] = z * cos_beta
        matrix[1][1] = y * sin_gamma
        matrix[1][2] = z * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
        matrix[2][2] = z * np.sqrt(
            sin_beta**2 - (cos_alpha - cos_beta * cos_gamma)**2 / sin_gamma**2
        )

        return matrix

    @property
    def bounding_edges(self) -> Np2DNumberArray:
        """Np2DNumberArray: The 8 corners of the bounding box of the unit cell."""
        edges = np.zeros((8, 3))
        for i, x in enumerate([-0.5, 0.5]):
            for j, y in enumerate([-0.5, 0.5]):
                for k, z in enumerate([-0.5, 0.5]):
                    edges[i * 4 + j * 2 +
                          k, :] = self.box_matrix @ np.array([x, y, z])

        return edges

    @property
    def volume(self) -> Real:
        """volume: The volume of the unit cell."""
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="overflow encountered in det"
            )
            volume = np.linalg.det(self.box_matrix)

        return volume

    @property
    def is_vacuum(self) -> bool:
        """bool: Returns whether the unit cell is a vacuum."""
        return bool(self.volume > 1e100)

    @runtime_type_checking
    def image(self, pos: NpnDNumberArray) -> NpnDNumberArray:
        """
        Images the given position(s) into the unit cell.

        This class can be used to image positions of arbitrary shape into the unit cell.
        The shape of the input is preserved. The only requirement is that the last 
        dimension of the input is 3, representing the x, y and z coordinates of the position(s).

        Parameters
        ----------
        pos : NpnDNumberArray
            The position to get the image of.

        Returns
        -------
        imaged_positions: NpnDNumberArray
            The image of the position(s) in the unit cell.
        """

        original_shape = np.shape(pos)
        pos = np.reshape(pos, (-1, 3))

        if self.alpha == 90 and self.beta == 90 and self.gamma == 90:
            pos = pos - self.box_lengths * np.round(pos / self.box_lengths)
        else:

            fractional_pos = pos @ self.inverse_box_matrix.T

            fractional_pos -= np.round(fractional_pos)

            pos = fractional_pos @ self.box_matrix.T

        return np.reshape(pos, original_shape)

    def __eq__(self, other: Any) -> Bool:
        """
        Checks if the Cell is equal to another Cell.

        Parameters
        ----------
        other : Cell
            The Cell to compare with.

        Returns
        -------
        Bool
            True if the Cells are equal, False otherwise.
        """

        return self.isclose(other)

    @runtime_type_checking
    def isclose(
        self,
        other: Any,
        rtol: PositiveReal = 1e-9,
        atol: PositiveReal = 0.0,
    ) -> Bool:
        """
        Checks if the Cell is close to another Cell.

        Parameters
        ----------
        other : Cell
            The Cell to compare with.
        rtol : PositiveReal, optional
            The relative tolerance parameter. Default is 1e-9.
        atol : PositiveReal, optional
            The absolute tolerance parameter. Default is 0.0.

        Returns
        -------
        Bool
            True if the Cells are close, False otherwise.
        """

        if not isinstance(other, Cell):
            return False

        is_equal = allclose_vectorized(
            self.box_lengths,
            other.box_lengths,
            rtol=rtol,
            atol=atol,
        )

        is_equal &= allclose_vectorized(
            self.box_angles,
            other.box_angles,
            rtol=rtol,
            atol=atol,
        )

        return is_equal

    def __str__(self) -> str:
        """
        Returns a string representation of the Cell.

        Returns
        -------
        str
            A string representation of the Cell.
        """
        x = float(self.x)
        y = float(self.y)
        z = float(self.z)
        alpha = float(self.alpha)
        beta = float(self.beta)
        gamma = float(self.gamma)

        if self != Cell():
            return f"Cell(x={x}, y={y}, z={z}, alpha={alpha}, beta={beta}, gamma={gamma})"

        return "Cell()"

    def __repr__(self) -> str:
        """
        Returns a string representation of the Cell.

        Returns
        -------
        str
            A string representation of the Cell.
        """
        return self.__str__()

    @classmethod
    @runtime_type_checking
    def init_from_box_matrix(cls, box_matrix: Np3x3NumberArray) -> "Cell":
        """
        Initializes a Cell object from a box matrix.

        Parameters
        ----------
        box_matrix : Np3x3NumberArray
            The box matrix.

        Returns
        -------
        Cell
            The Cell object.
        """
        x = np.linalg.norm(np.transpose(box_matrix)[0])
        y = np.linalg.norm(np.transpose(box_matrix)[1])
        z = np.linalg.norm(np.transpose(box_matrix)[2])

        gamma = np.arccos(box_matrix[0][1] / y)
        beta = np.arccos(box_matrix[0][2] / z)
        alpha = np.arccos(
            (
                box_matrix[0][1] * box_matrix[0][2] +
                box_matrix[1][1] * box_matrix[1][2]
            ) / (y * z)
        )

        return cls(
            x, y, z, np.rad2deg(alpha), np.rad2deg(beta), np.rad2deg(gamma)
        )



#: A type hint for a list of cells
Cells = NewType(
    "Cells",
    Annotated[list,
              Is[lambda list: all(isinstance(atom, Cell) for atom in list)]]
)
