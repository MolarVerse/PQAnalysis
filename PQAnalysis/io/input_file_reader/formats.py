"""
A module containing the InputFileFormat enumeration.
"""

from beartype.typing import Any

from PQAnalysis.formats import BaseEnumFormat
from .exceptions import InputFileFormatError



class InputFileFormat(BaseEnumFormat):

    """
    An enumeration of the supported input formats.
    """

    PQANALYSIS = "PQANALYSIS"
    PQ = "PQ"
    QMCFC = "QMCFC"

    @classmethod
    def _missing_(cls, value: object) -> Any:  # pylint: disable=arguments-differ
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
        return value in [cls.PQ, cls.QMCFC]
