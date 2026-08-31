"""
A module containing the JCoupling class.
"""

from PQAnalysis.types import PositiveInt
from PQAnalysis.type_checking import runtime_type_checking



class JCoupling:

    """
    A class to represent a j-coupling in a molecular topology.
    """

    # pylint: disable-next=too-many-arguments,too-many-positional-arguments
    @runtime_type_checking
    def __init__(
        self,
        index1: PositiveInt,
        index2: PositiveInt,
        index3: PositiveInt,
        index4: PositiveInt,
        j_coupling_type: PositiveInt | None = None,
        comment: str | None = None
    ) -> None:
        """
        Parameters
        ----------
        index1 : PositiveInt
            The index of the first atom in the j-coupling.
        index2 : PositiveInt
            The index of the second atom in the j-coupling.
        index3 : PositiveInt
            The index of the third atom in the j-coupling.
        index4 : PositiveInt
            The index of the fourth atom in the j-coupling.
        j_coupling_type : PositiveInt, optional
            The type of the j-coupling, by default None.
        comment : str, optional
            A comment for the j-coupling, by default None.
        """

        self.index1 = index1
        self.index2 = index2
        self.index3 = index3
        self.index4 = index4
        self.j_coupling_type = j_coupling_type
        self.comment = comment

    def copy(self) -> "JCoupling":
        """
        A method to create a copy of the j-coupling.

        Returns
        -------
        JCoupling
            A copy of the j-coupling.
        """
        return JCoupling(
            index1=self.index1,
            index2=self.index2,
            index3=self.index3,
            index4=self.index4,
            j_coupling_type=self.j_coupling_type,
            comment=self.comment
        )

    def __eq__(self, value: object) -> bool:
        """
        Compare the JCoupling object with another object.

        Parameters
        ----------
        value : object
            The object to compare with the JCoupling object.

        Returns
        -------
        bool
            True if the objects are equal, False otherwise.
        """

        if not isinstance(value, JCoupling):
            return False

        return (
            self.index1 == value.index1 and self.index2 == value.index2 and
            self.index3 == value.index3 and self.index4 == value.index4 and
            self.j_coupling_type == value.j_coupling_type
        )
