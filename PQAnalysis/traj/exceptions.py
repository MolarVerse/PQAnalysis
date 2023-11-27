"""
A module containing different exceptions related to the traj subpackage.

...

Classes
-------
FormatEnumError
    Exception raised if the given enum is not valid
TrajectoryFormatError
    Exception raised if the given enum is not valid
MDEngineFormatError
    Exception raised if the given enum is not valid
"""

from multimethod import multimethod

from PQAnalysis.exceptions import PQException


class FormatEnumError(PQException):
    """
    Exception raised if the given enum is not valid
    """

    @multimethod
    def __init__(self, value: object, enum: object) -> None:
        self.enum = enum
        self.value = value
        self.message = f"""
'{self.value}' is not a valid {enum.__name__}.
Possible values are: {enum.member_repr()}
or their case insensitive string representation: {enum.value_repr()}"""
        super().__init__(self.message)

    @multimethod
    def __init__(self, message: str) -> None:
        super().__init__(message)


class TrajectoryFormatError(FormatEnumError):
    """
    Exception raised if the given enum is not valid
    """

    @multimethod
    def __init__(self, value: object, enum: object) -> None:
        super().__init__(value, enum)

    @multimethod
    def __init__(self, message: str) -> None:
        super().__init__(message)


class MDEngineFormatError(FormatEnumError):
    """
    Exception raised if the given enum is not valid
    """

    @multimethod
    def __init__(self, value: object, enum: object) -> None:
        super().__init__(value, enum)

    @multimethod
    def __init__(self, message: str) -> None:
        super().__init__(message)
