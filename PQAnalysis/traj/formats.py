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

from . import TrajectoryFormatError, MDEngineFormatError, BoxFileFormatError
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


class BoxFileFormat(Format):
    """
    An enumeration of the supported box file formats.

    """

    #: The VMD format.
    #: The format looks in general like this:
    #:            8
    #:            Box  1.0 1.0 1.0    90.0 90.0 90.0
    #:            X   0.0 0.0 0.0
    #:            X   1.0 0.0 0.0
    #:            X   0.0 1.0 0.0
    #:            X   1.0 1.0 0.0
    #:            X   0.0 0.0 1.0
    #:            X   1.0 0.0 1.0
    #:            X   0.0 1.0 1.0
    #:            X   1.0 1.0 1.0
    #:            8
    #:            Box  1.0 1.0 1.0    90.0 90.0 90.0
    #:            X   0.0 0.0 0.0
    #:            ...
    #: where all X represent the vertices of the box. The first line contains the number of vertices.
    #: The second line contains the box dimensions and box angles as the comment line for a xyz file.
    VMD = "vmd"

    #: The data file format.
    #: The format looks in general like this:
    #:            1 1.0 1.0 1.0 90.0 90.0 90.0
    #:            2 1.0 1.0 1.0 90.0 90.0 90.0
    #:            ...
    #:            n 1.1 1.1 1.1 90.0 90.0 90.0
    #: where the first column represents the step starting from 1, the second to fourth column
    #: represent the box vectors a, b, c, the fifth to seventh column represent the box angles.
    DATA = "data"

    @classmethod
    def _missing_(cls, value: object) -> Any:
        """
        This method allows a BoxFileFormat to be retrieved from a string.
        """
        value = value.lower()
        for member in cls:
            if member.value.lower() == value:
                return member

        raise BoxFileFormatError(value, cls)
