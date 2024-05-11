"""
A module containing the base class for all Format enumerations.
"""

from enum import Enum
from beartype.typing import Any, List



class BaseEnumFormat(Enum):

    """
    An enumeration super class of the various supported trajectory formats.
    """

    @classmethod
    def member_repr(cls) -> str:
        """
        This method returns a string representation of the members of the enumeration.

        Returns
        -------
        str
            A string representation of the members of the enumeration.
        """

        return ', '.join([str(member) for member in cls])

    @classmethod
    def value_repr(cls) -> str:
        """
        This method returns a string representation of
        the values of the members of the enumeration.

        Returns
        -------
        str
            A string representation of the values of the members of the enumeration.
        """

        return ', '.join([str(member.value) for member in cls])

    @classmethod
    def _missing_(cls, value: object, exception: type(Exception)) -> Any:  # pylint: disable=arguments-differ
        """
        This method allows a FileWriteMode to be retrieved from a string.

        Parameters
        ----------
        value : object
            The value to return.
        exception : Exception
            The exception to raise if the value is not found.

        Raises
        ------
        exception
            If the value is not found.
        """
        value = value.lower()

        for member in cls:
            if member.value.lower() == value:
                return member

        raise exception(value, cls)

    @classmethod
    def values(cls) -> List[str]:
        """
        This method returns a list of all values of the enumeration.

        Returns
        -------
        List[str]
            A list of all values of the enumeration.
        """

        return [member.value for member in cls]

    def __eq__(self, other: object) -> bool:
        """
        Checks if the given EnumFormat is equal to this Format or not.
        """

        if not isinstance(other, type(self)) and not isinstance(other, str):
            return False

        other = type(self)(other)

        return self.value == other.value
