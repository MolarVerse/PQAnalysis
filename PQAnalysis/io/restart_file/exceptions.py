"""
A module containing different exceptions related to the restart_file subpackage.
"""

from PQAnalysis.exceptions import PQException



class RestartFileReaderError(PQException):

    """
    Exception raised for errors related to the RestartFileReader class
    """



class RestartFileWriterError(PQException):

    """
    Exception raised for errors related to the RestartFileWriter class
    """
