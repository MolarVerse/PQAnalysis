from beartype.typing import Any

from ...formats import Format
from .exceptions import InputFileFormatError


class InputFileFormat(Format):
    """
    An enumeration of the supported input formats.

    ...

    Attributes
    ----------
    PQANALYSIS : str
        The PQANALYSIS format.
    PIMD-QMCF : str
        The PIMD-QMCF format.
    QMCFC : str
        The QMCFC format.
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
