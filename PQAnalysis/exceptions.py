"""
A module containing different exceptions which could be useful.

...

Classes
-------
PQException
    Base class for exceptions in this module.
ElementNotFoundError
    Exception raised if the given element id is not valid
TrajectoryFormatError
    Exception raised if the given enum is not valid
"""

from beartype.typing import Any


class PQException(Exception):
    """
    Base class for exceptions in this module.
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)


class ElementNotFoundError(PQException):
    """
    Exception raised if the given element id is not valid
    """

    def __init__(self, id: Any) -> None:
        self.id = id
        self.message = f"""Id {self.id} is not a valid element identifier."""
        super().__init__(self.message)


class FormatEnumError(PQException):
    """
    Exception raised if the given enum is not valid
    """

    def __init__(self, value: object, enum: object) -> None:
        self.enum = enum
        self.value = value
        self.message = f"""
'{self.value}' is not a valid {enum.__name__}.
Possible values are: {enum.member_repr()}
or their case insensitive string representation: {enum.value_repr()}"""
        super().__init__(self.message)


class TrajectoryFormatError(FormatEnumError):
    """
    Exception raised if the given enum is not valid
    """

    def __init__(self, value: object, enum: object) -> None:
        super().__init__(value, enum)


class MDEngineFormatError(FormatEnumError):
    """
    Exception raised if the given enum is not valid
    """

    def __init__(self, value: object, enum: object) -> None:
        super().__init__(value, enum)
