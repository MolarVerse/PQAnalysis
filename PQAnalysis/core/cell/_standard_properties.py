"""
A module containing a Mixin Class with the standard properties
of a Cell class (i.e. standard getter and setter methods).
"""
from numbers import Real

import numpy as np

from PQAnalysis.types import Np1DNumberArray, Np3x3NumberArray
from PQAnalysis.type_checking import runtime_type_checking_setter



class _StandardPropertiesMixin:

    """
    A mixin class containing the standard properties of a
    Cell class (i.e. standard getter and setter methods).
    """

    @property
    def box_lengths(self) -> Np1DNumberArray:
        """
        Np1DNumberArray: The lengths of the box vectors.

        When setting the box lengths, the box matrix is recalculated.
        """
        return self._box_lengths

    @box_lengths.setter
    @runtime_type_checking_setter
    def box_lengths(self, box_lengths: Np1DNumberArray) -> None:
        self._box_lengths = box_lengths
        self._box_matrix = self.setup_box_matrix()

    @property
    def box_angles(self) -> Np1DNumberArray:
        """
        Np1DNumberArray: The angles between the box vectors.

        When setting the box angles, the box matrix is recalculated.
        """
        return self._box_angles

    @box_angles.setter
    @runtime_type_checking_setter
    def box_angles(self, box_angles: Np1DNumberArray) -> None:
        self._box_angles = box_angles
        self._box_matrix = self.setup_box_matrix()

    @property
    def x(self) -> Real:
        """Real: The length of the first box vector."""
        return self._box_lengths[0]

    @property
    def y(self) -> Real:
        """Real: The length of the second box vector."""
        return self._box_lengths[1]

    @property
    def z(self) -> Real:
        """Real: The length of the third box vector."""
        return self._box_lengths[2]

    @property
    def alpha(self) -> Real:
        """Real: The angle between the second and third box vector."""
        return self._box_angles[0]

    @property
    def beta(self) -> Real:
        """Real: The angle between the first and third box vector."""
        return self._box_angles[1]

    @property
    def gamma(self) -> Real:
        """Real: The angle between the first and second box vector."""
        return self._box_angles[2]

    @property
    def box_matrix(self) -> Np3x3NumberArray:
        """Np3x3NumberArray: The box matrix."""
        return self._box_matrix

    @property
    def _box_matrix(self) -> Np3x3NumberArray:
        """
        Np3x3NumberArray: The box matrix.

        This property is used internally to set the box matrix.
        It is not intended to be used by the user. Use the property 'box_matrix' instead.

        When setting the box matrix, the inverse box matrix is calculated.
        """
        return self.__box_matrix

    @_box_matrix.setter
    @runtime_type_checking_setter
    def _box_matrix(self, box_matrix: Np3x3NumberArray) -> None:
        self.__box_matrix = box_matrix
        self._inverse_box_matrix = np.linalg.inv(box_matrix)

    @property
    def inverse_box_matrix(self) -> Np3x3NumberArray:
        """Np3x3NumberArray: The inverse box matrix."""
        return self._inverse_box_matrix
