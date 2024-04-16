"""
A module containing the Dihedral class.
"""

from PQAnalysis.types import PositiveInt, PositiveReal


class Dihedral:
    """
    A class to represent a dihedral in a molecular topology.
    """

    def __init__(self,
                 index1: PositiveInt,
                 index2: PositiveInt,
                 index3: PositiveInt,
                 index4: PositiveInt,
                 equilibrium_angle: PositiveReal | None = None,
                 dihedral_type: int | None = None,
                 is_linker: bool = False,
                 is_improper: bool = False,
                 ) -> None:
        """
        Parameters
        ----------
        index1 : PositiveInt
            The index of the first atom in the dihedral.
        index2 : PositiveInt
            The index of the second atom in the dihedral.
        index3 : PositiveInt
            The index of the third atom in the dihedral.
        index4 : PositiveInt
            The index of the fourth atom in the dihedral.
        equilibrium_angle : PositiveReal, optional
            The equilibrium angle of the dihedral, by default None.
        dihedral_type : int, optional
            The type of the dihedral, by default None.
        is_linker : bool, optional
            A flag to indicate if the dihedral is a linker, by default False.
        is_improper : bool, optional
            A flag to indicate if the dihedral is an improper dihedral, by default False.
        """

        self.index1 = index1
        self.index2 = index2
        self.index3 = index3
        self.index4 = index4
        self.equilibrium_angle = equilibrium_angle
        self.angle_type = dihedral_type
        self.is_linker = is_linker
        self.is_improper = is_improper
