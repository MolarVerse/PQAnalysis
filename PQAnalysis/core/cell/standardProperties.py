"""
A module containing a Mixin Class with the standard properties of a Cell class (i.e. standard getter and setter methods).
"""
import numpy as np

from numbers import Real

from ...types import Np1DNumberArray, Np3x3NumberArray


class _StandardPropertiesMixin:
    """
    A mixin class containing the standard properties of a Cell class (i.e. standard getter and setter methods).
    """
    @property
    def box_lengths(self) -> Np1DNumberArray:
        """
        The lengths of the box vectors.

        Returns
        -------
        box_lengths: Np1DNumberArray of shape (3,)
            The lengths of the box vectors.
        """
        return self._box_lengths

    @box_lengths.setter
    def box_lengths(self, box_lengths: Np1DNumberArray) -> None:
        """
        Sets the box lengths and calculates the box matrix.

        Parameters
        ----------
        box_lengths : Np1DNumberArray
            The lengths of the box vectors.
        """
        self._box_lengths = box_lengths
        self._box_matrix = self.setup_box_matrix()

    @property
    def box_angles(self) -> Np1DNumberArray:
        """
        The angles between the box vectors.

        Returns
        -------
        box_angles: Np1DNumberArray of shape (3,)
            The angles between the box vectors.
        """
        return self._box_angles

    @box_angles.setter
    def box_angles(self, box_angles: Np1DNumberArray) -> None:
        """
        Sets the box angles and calculates the box matrix.

        Parameters
        ----------
        box_angles : Np1DNumberArray
            The angles between the box vectors.
        """
        self._box_angles = box_angles
        self._box_matrix = self.setup_box_matrix()

    @property
    def x(self) -> Real:
        """
        The length of the first box vector.

        Returns
        -------
        x: Real
            The length of the first box vector.
        """
        return self._box_lengths[0]

    @property
    def y(self) -> Real:
        """
        The length of the second box vector.

        Returns
        -------
        y: Real
            The length of the second box vector.
        """
        return self._box_lengths[1]

    @property
    def z(self) -> Real:
        """
        The length of the third box vector.

        Returns
        -------
        z: Real
            The length of the third box vector.
        """
        return self._box_lengths[2]

    @property
    def alpha(self) -> Real:
        """
        The angle between the second and third box vector.

        Returns
        -------
        alpha: Real
            The angle between the second and third box vector.
        """
        return self._box_angles[0]

    @property
    def beta(self) -> Real:
        """
        The angle between the first and third box vector.

        Returns
        -------
        beta: Real
            The angle between the first and third box vector.
        """
        return self._box_angles[1]

    @property
    def gamma(self) -> Real:
        """
        The angle between the first and second box vector.

        Returns
        -------
        gamma: Real
            The angle between the first and second box vector.
        """
        return self._box_angles[2]

    @property
    def box_matrix(self) -> Np3x3NumberArray:
        """
        The matrix containing the box vectors as columns.

        Returns
        -------
        box_matrix: Np3x3NumberArray
            The box matrix.
        """
        return self._box_matrix

    @property
    def _box_matrix(self) -> Np3x3NumberArray:
        """
        The matrix containing the box vectors as columns.

        Returns
        -------
        __box_matrix: Np3x3NumberArray
            The box matrix.
        """
        return self.__box_matrix

    @_box_matrix.setter
    def _box_matrix(self, box_matrix: Np3x3NumberArray) -> None:
        """
        Sets the box matrix and calculates the inverse box matrix.

        Parameters
        ----------
        box_matrix : Np3x3NumberArray
            The box matrix.
        """
        self.__box_matrix = box_matrix
        self._inverse_box_matrix = np.linalg.inv(box_matrix)

    @property
    def inverse_box_matrix(self) -> Np3x3NumberArray:
        """
        The inverse box matrix.

        Returns
        -------
        inverse_box_matrix: Np3x3NumberArray
            The inverse box matrix.
        """
        return self._inverse_box_matrix
