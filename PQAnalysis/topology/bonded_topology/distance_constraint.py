"""
A module containing the DistanceConstraint class.
"""

from numbers import Real

from PQAnalysis.types import PositiveInt, PositiveReal
from PQAnalysis.type_checking import runtime_type_checking



class DistanceConstraint:

    """
    A class to represent a distance constraint in a molecular topology.
    """

    # pylint: disable-next=too-many-arguments,too-many-positional-arguments
    @runtime_type_checking
    def __init__(
        self,
        index1: PositiveInt,
        index2: PositiveInt,
        lower_distance: PositiveReal,
        upper_distance: PositiveReal,
        spring_constant: PositiveReal,
        d_spring_constant_dt: Real,
        comment: str | None = None
    ) -> None:
        """
        Parameters
        ----------
        index1 : PositiveInt
            The index of the first atom in the distance constraint.
        index2 : PositiveInt
            The index of the second atom in the distance constraint.
        lower_distance : PositiveReal
            The lower distance of the distance constraint.
        upper_distance : PositiveReal
            The upper distance of the distance constraint.
        spring_constant : PositiveReal
            The spring constant of the distance constraint.
        d_spring_constant_dt : Real
            The time derivative of the spring constant.
        comment : str, optional
            A comment for the distance constraint, by default None.
        """

        self.index1 = index1
        self.index2 = index2
        self.lower_distance = lower_distance
        self.upper_distance = upper_distance
        self.spring_constant = spring_constant
        self.d_spring_constant_dt = d_spring_constant_dt
        self.comment = comment

    def copy(self) -> "DistanceConstraint":
        """
        A method to create a copy of the distance constraint.

        Returns
        -------
        DistanceConstraint
            A copy of the distance constraint.
        """
        return DistanceConstraint(
            index1=self.index1,
            index2=self.index2,
            lower_distance=self.lower_distance,
            upper_distance=self.upper_distance,
            spring_constant=self.spring_constant,
            d_spring_constant_dt=self.d_spring_constant_dt,
            comment=self.comment
        )

    def __eq__(self, value: object) -> bool:
        """
        Compare the DistanceConstraint object with another object.

        Parameters
        ----------
        value : object
            The object to compare with the DistanceConstraint object.

        Returns
        -------
        bool
            True if the objects are equal, False otherwise.
        """

        if not isinstance(value, DistanceConstraint):
            return False

        return (
            self.index1 == value.index1 and self.index2 == value.index2 and
            self.lower_distance == value.lower_distance and
            self.upper_distance == value.upper_distance and
            self.spring_constant == value.spring_constant and
            self.d_spring_constant_dt == value.d_spring_constant_dt
        )
