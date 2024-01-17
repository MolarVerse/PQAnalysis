"""
A module containing the InputFileFormat enumeration.

...

Classes
-------
InputFileFormat
    An enumeration of the supported input formats.
"""

from beartype.typing import Any

from .exceptions import InputFileFormatError
from PQAnalysis.formats import BaseEnumFormat


class InputFileFormat(BaseEnumFormat):
    """
    An enumeration of the supported input formats.
    """

    PQANALYSIS = "PQANALYSIS"
    PIMD_QMCF = "PIMD-QMCF"
    QMCFC = "QMCFC"

    @classmethod
    def _missing_(cls, value: object) -> Any:
        """
        This method allows an InputFileFormat to be retrieved from a string.
        """
        value = value.lower()
        for member in cls:
            if member.value.lower() == value:
                return member

        raise InputFileFormatError(value, cls)

    @classmethod
    def is_qmcf_type(cls, value: object) -> bool:
        """
        Returns True if the given value is a QMCF input file format.
        """
        return value in [cls.PIMD_QMCF, cls.QMCFC]
