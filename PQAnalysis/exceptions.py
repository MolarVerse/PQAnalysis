"""
A module containing different exceptions which could be useful.
"""

from multimethod import multimethod


class PQException(Exception):
    """
    Base class for exceptions in this package.
    """

    def __init__(self, message: str) -> None:
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        self.message = message
        super().__init__(self.message)


class PQWarning(Warning):
    """
    Base class for warnings in this package.
    """

    def __init__(self, message: str) -> None:
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        self.message = message
        super().__init__(self.message)


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
        self.message = f"""
'{self.value}' is not a valid {enum.__name__}.
Possible values are: {enum.member_repr()}
or their case insensitive string representation: {enum.value_repr()}"""
        super().__init__(self.message)

    @multimethod
    def __init__(self, message: str) -> None:
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        super().__init__(message)
