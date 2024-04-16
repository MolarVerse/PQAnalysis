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
        """
        Parameters
        ----------
        value : object
            The value that is not valid.
        enum : object
            The enum that is not valid.
        """
        super().__init__(value, enum)

    @multimethod
    def __init__(self, message: str) -> None:
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        super().__init__(message)


class InputFileError(PQException):
    """
    Exception raised if something is wrong with the input file
    """

    def __init__(self, message: str) -> None:
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        self.message = message
        super().__init__(message)


class InputFileWarning(PQWarning):
    """
    Warning raised if something is wrong with the input file
    """

    def __init__(self, message: str) -> None:
        """
        Parameters
        ----------
        message : str
            The warning message.
        """
        self.message = message
        super().__init__(message)
