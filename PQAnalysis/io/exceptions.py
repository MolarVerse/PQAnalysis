"""
A module containing different exceptions related to the io subpackage.
"""

from PQAnalysis.exceptions import PQException, BaseEnumFormatError



class BoxWriterError(PQException):

    """
    Exception raised for errors related to the BoxWriter class
    """



class MoldescriptorReaderError(PQException):

    """
    Exception raised for errors related to the MoldescriptorReader class
    """



class BoxFileFormatError(BaseEnumFormatError):

    """
    Exception raised if the given enum is not valid
    """



class FileWritingModeError(BaseEnumFormatError):

    """
    Exception raised if the given enum is not valid
    """



class OutputFileFormatError(BaseEnumFormatError):

    """
    Exception raised if the given enum is not valid
    """
