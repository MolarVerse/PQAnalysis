"""
A module containing different exceptions related to the io subpackage.
"""

from PQAnalysis.exceptions import PQException, BaseEnumFormatError


class BoxWriterError(PQException):
    """
    Exception raised for errors related to the BoxWriter class
    """

    pass


class FrameReaderError(PQException):
    """
    Exception raised for errors related to the FrameReader class
    """

    pass


class MoldescriptorReaderError(PQException):
    """
    Exception raised for errors related to the MoldescriptorReader class
    """

    pass


class RestartFileReaderError(PQException):
    """
    Exception raised for errors related to the RestartFileReader class
    """

    pass


class TrajectoryReaderError(PQException):
    """
    Exception raised for errors related to the TrajectoryReader class
    """

    pass


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
