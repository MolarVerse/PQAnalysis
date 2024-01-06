"""
A module containing different format types of the trajectory.

...

Classes
-------
TrajectoryFormat
    An enumeration of the supported trajectory formats.
MDEngineFormat
    An enumeration of the supported MD engine formats.
"""

from beartype.typing import Any

from . import TrajectoryFormatError, MDEngineFormatError
from ..formats import Format


class TrajectoryFormat(Format):
    """
    An enumeration of the supported trajectory formats.
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
