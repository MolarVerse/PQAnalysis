"""
A module containing the Angle class.
"""

from PQAnalysis.types import PositiveInt, PositiveReal


class Angle:
    """
    A class to represent an angle in a molecular topology.
    """

    def __init__(self,
                 index1: PositiveInt,
                 index2: PositiveInt,
                 index3: PositiveInt,
                 equilibrium_angle: PositiveReal | None = None,
                 angle_type: PositiveInt | None = None,
                 is_linker: bool = False,
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
        """

        self.index1 = index1
        self.index2 = index2
        self.index3 = index3
        self.equilibrium_angle = equilibrium_angle
        self.angle_type = angle_type
        self.is_linker = is_linker
