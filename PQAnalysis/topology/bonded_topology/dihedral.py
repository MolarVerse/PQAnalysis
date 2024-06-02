"""
A module containing the Dihedral class.
"""

from PQAnalysis.types import PositiveInt, PositiveReal
from PQAnalysis.type_checking import runtime_type_checking



class Dihedral:

    """
    A class to represent a dihedral in a molecular topology.
    """

    @runtime_type_checking
    def __init__(
        self,
        index1: PositiveInt,
        index2: PositiveInt,
        index3: PositiveInt,
        index4: PositiveInt,
        equilibrium_angle: PositiveReal | None = None,
        dihedral_type: int | None = None,
        is_linker: bool = False,
        is_improper: bool = False,
        comment: str | None = None
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
        comment : str, optional
            A comment for the dihedral, by default None.
        """

        self.index1 = index1
        self.index2 = index2
        self.index3 = index3
        self.index4 = index4
        self.equilibrium_angle = equilibrium_angle
        self.dihedral_type = dihedral_type
        self.is_linker = is_linker
        self.is_improper = is_improper
        self.comment = comment

    def copy(self) -> "Dihedral":
        """
        A method to create a copy of the dihedral.

        Returns
        -------
        Dihedral
            A copy of the dihedral.
        """
        return Dihedral(
            index1=self.index1,
            index2=self.index2,
            index3=self.index3,
            index4=self.index4,
            equilibrium_angle=self.equilibrium_angle,
            dihedral_type=self.dihedral_type,
            is_linker=self.is_linker,
            is_improper=self.is_improper,
            comment=self.comment
        )

    def __eq__(self, value: object) -> bool:
        """
        Compare the Dihedral object with another object.

        Parameters
        ----------
        value : object
            The object to compare with the Dihedral object.

        Returns
        -------
        bool
            True if the objects are equal, False otherwise.
        """

        if not isinstance(value, Dihedral):
            return False

        return (
            self.index1 == value.index1 and self.index2 == value.index2 and
            self.index3 == value.index3 and self.index4 == value.index4 and
            self.equilibrium_angle == value.equilibrium_angle and
            self.dihedral_type == value.dihedral_type and
            self.is_linker == value.is_linker and
            self.is_improper == value.is_improper
        )
