"""
A module containing different exceptions which could be useful.
"""

from multimethod import multimethod



class PQException(Exception):

    """
    Base class for exceptions in this package.
    """



class PQWarning(Warning):

    """
    Base class for warnings in this package.
    """



class PQIndexError(PQException):

    """
    Exception raised for errors related to the AtomicSystem class
    """



class PQTypeError(PQException):

    """
    Exception raised for errors related to the AtomicSystem class
    """



class PQNotImplementedError(PQException):

    """
    Exception raised for errors related to the AtomicSystem class
    """



class PQFileNotFoundError(PQException):

    """
    Exception raised for errors related to the AtomicSystem class
    """



class PQKeyError(PQException):

    """
    Exception raised for errors related to the AtomicSystem class
    """



class PQValueError(PQException):

    """
    Exception raised for errors related to the AtomicSystem class
    """



class BaseEnumFormatError(PQException):

    """
    Base class for enum exceptions if the given enum is not valid
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

        self.enum = enum
        self.value = value
        self.message = (
            f"\n"
            f"\'{self.value}\' is not a valid "
            f"{enum.__name__}.\n"
            f"Possible values are: {enum.member_repr()} or their "
            f"case insensitive string representation: {enum.value_repr()}"
        )
        super().__init__(self.message)

    @multimethod
    def __init__(self, message: str) -> None:  # pylint: disable=function-redefined
        super().__init__(message)
