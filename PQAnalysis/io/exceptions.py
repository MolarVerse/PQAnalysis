"""
A module containing different exceptions related to the io subpackage.
"""

from multimethod import multimethod

from PQAnalysis.exceptions import PQException, BaseEnumFormatError


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


class BoxFileFormatError(BaseEnumFormatError):
    """
    Exception raised if the given enum is not valid
    """

    pass


class FileWritingModeError(BaseEnumFormatError):
    """
    Exception raised if the given enum is not valid
    """

    pass
