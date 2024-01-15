"""
A module containing different exceptions related to the io subpackage.
"""

from multimethod import multimethod

from ..exceptions import PQException, FormatEnumError


class BoxWriterError(PQException):
    """
    Exception raised for errors related to the BoxWriter class
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)


class FrameReaderError(PQException):
    """
    Exception raised for errors related to the FrameReader class
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)


class MoldescriptorReaderError(PQException):
    """
    Exception raised for errors related to the MoldescriptorReader class
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)


class RestartFileReaderError(PQException):
    """
    Exception raised for errors related to the RestartFileReader class
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)


class TrajectoryReaderError(PQException):
    """
    Exception raised for errors related to the TrajectoryReader class
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)


class BoxFileFormatError(FormatEnumError):
    """
    Exception raised if the given enum is not valid
    """

    @multimethod
    def __init__(self, value: object, enum: object) -> None:
        super().__init__(value, enum)

    @multimethod
    def __init__(self, message: str) -> None:
        super().__init__(message)
