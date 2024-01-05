"""
A module containing exceptions related to input file errors

Classes
-------
InputFileFormatError
    Exception raised if the given enum is not valid
InputFileError
    Exception raised if something is wrong with the input file
InputFileWarning
    Warning raised if something is wrong with the input file
"""

from multimethod import multimethod

from ...exceptions import FormatEnumError, PQException, PQWarning


class InputFileFormatError(FormatEnumError):
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
