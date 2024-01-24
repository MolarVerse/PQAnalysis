"""
A module containing exceptions related to input file errors
"""

from multimethod import multimethod

from PQAnalysis.exceptions import BaseEnumFormatError, PQException, PQWarning


class InputFileFormatError(BaseEnumFormatError):
    """
    Exception raised if the given enum is not valid
    """

    @multimethod
    def __init__(self, value: object, enum: object) -> None:
        super().__init__(value, enum)

    @multimethod
    def __init__(self, message: str) -> None:
        super().__init__(message)


class InputFileError(PQException):
    """
    Exception raised if something is wrong with the input file
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(message)


class InputFileWarning(PQWarning):
    """
    Warning raised if something is wrong with the input file
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(message)
