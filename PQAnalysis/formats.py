"""
A module containing the base class for all Format enumerations.

...

Classes
-------
Format
    An enumeration super class of the various supported trajectory formats.
"""

from enum import Enum


class Format(Enum):
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
        This method returns a string representation of the values of the members of the enumeration.

        Returns
        -------
        str
            A string representation of the values of the members of the enumeration.
        """

        return ', '.join([str(member.value) for member in cls])
