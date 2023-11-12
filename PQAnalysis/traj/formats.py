"""
A module containing different format types of the trajectory.

...

Classes
-------
Format
    An enumeration super class of the various supported trajectory formats.
TrajectoryFormat
    An enumeration of the supported trajectory formats.
MDEngineFormat
    An enumeration of the supported MD engine formats.
"""

from enum import Enum
from beartype.typing import Any

from ..exceptions import TrajectoryFormatError, MDEngineFormatError


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


class TrajectoryFormat(Format):
    """
    An enumeration of the supported trajectory formats.

    ...

    Attributes
    ----------
    XYZ : str
        The XYZ format.
    VEL : str
        The VEL format.
    FORCE : str
        The FORCE format.
    CHARGE : str
        The CHARGE format.
    """

    XYZ = "XYZ"
    VEL = "VEL"
    FORCE = "FORCE"
    CHARGE = "CHARGE"

    @classmethod
    def _missing_(cls, value: object) -> Any:
        """
        This method allows a TrajectoryFormat to be retrieved from a string.
        """
        value = value.lower()
        for member in cls:
            if member.value.lower() == value:
                return member

        raise TrajectoryFormatError(value, cls)


class MDEngineFormat(Format):
    """
    An enumeration of the supported MD engine formats.

    ...

    Attributes
    ----------
    PIMD-QMCF: str
        The PIMD-QMCF format.
    QMCFC: str
        The QMCFC format.
    """

    PIMD_QMCF = "PIMD-QMCF"
    QMCFC = "QMCFC"

    @classmethod
    def _missing_(cls, value: object) -> Any:
        """
        This method allows an MDEngineFormat format to be retrieved from a string.
        """
        value = value.lower()
        for member in cls:
            if member.value.lower() == value:
                return member

        raise MDEngineFormatError(value, cls)

    @classmethod
    def isQMCFType(cls, format: Any) -> bool:
        """
        This method checks if the given format is a QMCF format.

        Parameters
        ----------
        format : Any
            The format to check.

        Returns
        -------
        bool
            True if the format is a QMCF format, False otherwise.
        """
        return format in [cls.PIMD_QMCF, cls.QMCFC]
