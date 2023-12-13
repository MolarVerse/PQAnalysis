"""
A module containing different exceptions related to the traj subpackage.

...

Classes
-------
TrajectoryFormatError
    Exception raised if the given enum is not valid
MDEngineFormatError
    Exception raised if the given enum is not valid
FrameError
    Exception raised for errors related to the Frame class
"""

from multimethod import multimethod

from ..exceptions import PQException, FormatEnumError


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


class FrameError(PQException):
    """
    Exception raised for errors related to the Frame class
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)
