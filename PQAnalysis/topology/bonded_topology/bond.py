"""
A module containing the Bond class.
"""

from __future__ import annotations

from PQAnalysis.types import PositiveInt, PositiveReal


class Bond:
    """
    A class to represent a bond in a molecular topology.
    """

    def __init__(self,
                 index1: PositiveInt,
                 index2: PositiveInt,
                 equilibrium_distance: PositiveReal | None = None,
                 bond_type: PositiveInt | None = None,
                 is_linker: bool = False,
                 is_shake: bool = False,
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
        """

        self.index1 = index1
        self.index2 = index2
        self.equilibrium_distance = equilibrium_distance
        self.bond_type = bond_type
        self.is_linker = is_linker
        self.is_shake = is_shake

    def copy(self) -> Bond:
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
        )
