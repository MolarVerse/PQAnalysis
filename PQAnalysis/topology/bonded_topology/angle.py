"""
A module containing the Angle class.
"""

from PQAnalysis.types import PositiveInt, PositiveReal
from PQAnalysis.type_checking import runtime_type_checking



class Angle:

    """
    A class to represent an angle in a molecular topology.
    """

    @runtime_type_checking
    def __init__(
        self,
        index1: PositiveInt,
        index2: PositiveInt,
        index3: PositiveInt,
        equilibrium_angle: PositiveReal | None = None,
        angle_type: PositiveInt | None = None,
        is_linker: bool = False,
        comment: str | None = None
    ) -> None:
        """
        Parameters
        ----------
        index1 : PositiveInt
            The index of the first atom in the angle.
        index2 : PositiveInt
            The index of the second atom in the angle.
        index3 : PositiveInt
            The index of the third atom in the angle.
        equilibrium_angle : PositiveReal, optional
            The equilibrium angle of the angle, by default None.
        angle_type : PositiveInt, optional
            The type of the angle, by default None.
        is_linker : bool, optional
            A flag to indicate if the angle is a linker, by default False.
        comment : str, optional
            A comment for the angle, by default None.
        """

        self.index1 = index1
        self.index2 = index2
        self.index3 = index3
        self.equilibrium_angle = equilibrium_angle
        self.angle_type = angle_type
        self.is_linker = is_linker
        self.comment = comment

    def copy(self) -> "Angle":
        """
        A method to create a copy of the angle.

        Returns
        -------
        Angle
            A copy of the angle.
        """
        return Angle(
            index1=self.index1,
            index2=self.index2,
            index3=self.index3,
            equilibrium_angle=self.equilibrium_angle,
            angle_type=self.angle_type,
            is_linker=self.is_linker,
            comment=self.comment
        )

    def __eq__(self, value: object) -> bool:
        """
        Compare the Angle object with another object.

        Parameters
        ----------
        value : object
            The object to compare with the Angle object.

        Returns
        -------
        bool
            True if the objects are equal, False otherwise.
        """

        if not isinstance(value, Angle):
            return False

        return (
            self.index1 == value.index1 and self.index2 == value.index2 and
            self.index3 == value.index3 and
            self.equilibrium_angle == value.equilibrium_angle and
            self.angle_type == value.angle_type and
            self.is_linker == value.is_linker
        )
