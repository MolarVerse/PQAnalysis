"""
A module containing different format types of the trajectory.
"""

from beartype.typing import Any

from . import TrajectoryFormatError, MDEngineFormatError
from PQAnalysis.formats import BaseEnumFormat


class TrajectoryFormat(BaseEnumFormat):
    """
    An enumeration of the supported trajectory formats.
    """

    XYZ = "XYZ"
    VEL = "VEL"
    FORCE = "FORCE"
    CHARGE = "CHARGE"

    @classmethod
    def _missing_(cls, value: Any) -> Any:
        """
        This method returns the missing value of the enumeration.

        Parameters
        ----------
        value : Any
            The value to return.

        Returns
        -------
        Any
            The value to return.
        """

        return super()._missing_(value, TrajectoryFormatError)


class MDEngineFormat(BaseEnumFormat):
    """
    An enumeration of the supported MD engine formats.
    """

    PIMD_QMCF = "PIMD-QMCF"
    QMCFC = "QMCFC"

    @classmethod
    def _missing_(cls, value: Any) -> Any:
        """
        This method returns the missing value of the enumeration.

        Parameters
        ----------
        value : Any
            The value to return.

        Returns
        -------
        Any
            The value to return.
        """

        return super()._missing_(value, MDEngineFormatError)

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
