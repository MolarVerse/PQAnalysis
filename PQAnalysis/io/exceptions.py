"""
A module containing different exceptions related to the io subpackage.
"""

from PQAnalysis.exceptions import PQException, BaseEnumFormatError


class BoxWriterError(PQException):
    """
    Exception raised for errors related to the BoxWriter class
    """

    pass


class MoldescriptorReaderError(PQException):
    """
    Exception raised for errors related to the MoldescriptorReader class
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


class OutputFileFormatError(BaseEnumFormatError):
    """
    Exception raised if the given enum is not valid
    """

    pass
