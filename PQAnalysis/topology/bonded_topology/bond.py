"""
A module containing the Bond class.
"""

from PQAnalysis.types import PositiveInt, PositiveReal
from PQAnalysis.type_checking import runtime_type_checking



class Bond:

    """
    A class to represent a bond in a molecular topology.
    """

    @runtime_type_checking
    def __init__(
        self,
        index1: PositiveInt,
        index2: PositiveInt,
        equilibrium_distance: PositiveReal | None = None,
        bond_type: PositiveInt | None = None,
        is_linker: bool = False,
        is_shake: bool = False,
        comment: str | None = None
    ) -> None:
        """
        Parameters
        ----------
        index1 : PositiveInt
            The index of the first atom in the bond.
        index2 : PositiveInt
            The index of the second atom in the bond.
        equilibrium_distance : PositiveReal, optional
            The equilibrium distance of the bond, by default None.
        bond_type : PositiveInt, optional
            The type of the bond, by default None.
        is_linker : bool, optional
            A flag to indicate if the bond is a linker, by default False.
        is_shake : bool, optional
            A flag to indicate if the bond is a shake bond, by default False.
        comment : str, optional
            A comment for the bond, by default None.
        """

        self.index1 = index1
        self.index2 = index2
        self.equilibrium_distance = equilibrium_distance
        self.bond_type = bond_type
        self.is_linker = is_linker
        self.is_shake = is_shake
        self.comment = comment

    def copy(self) -> "Bond":
        """
        A method to create a copy of the bond.

        Returns
        -------
        Bond
            A copy of the bond.
        """
        return Bond(
            index1=self.index1,
            index2=self.index2,
            equilibrium_distance=self.equilibrium_distance,
            bond_type=self.bond_type,
            is_linker=self.is_linker,
            is_shake=self.is_shake,
            comment=self.comment
        )

    def __eq__(self, value: object) -> bool:
        """
        Compare the bond with another bond.

        Parameters
        ----------
        value : object
            The bond to compare.

        Returns
        -------
        bool
            True if the bonds are equal, False otherwise.
        """

        if not isinstance(value, Bond):
            return False

        return (
            self.index1 == value.index1 and self.index2 == value.index2 and
            self.equilibrium_distance == value.equilibrium_distance and
            self.bond_type == value.bond_type and
            self.is_linker == value.is_linker and
            self.is_shake == value.is_shake
        )
