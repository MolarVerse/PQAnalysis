"""
A module containing exceptions related to input file errors
"""

from PQAnalysis.exceptions import BaseEnumFormatError, PQException, PQWarning



class InputFileFormatError(BaseEnumFormatError):

    """
    Exception raised if the given enum is not valid
    """

    def __init__(self, value: object, enum: object) -> None:
        """
        Parameters
        ----------
        value : object
            The invalid value.
        enum : object
            The enumeration.
        """
        self.value = value
        self.enum = enum
        super().__init__(value, enum)



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
