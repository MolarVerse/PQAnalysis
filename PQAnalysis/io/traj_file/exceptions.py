"""
A module containing different exceptions related to the traj_file subpackage.
"""

from PQAnalysis.exceptions import PQException



class TrajectoryReaderError(PQException):

    """
    Exception raised for errors related to the TrajectoryReader class
    """



class FrameReaderError(PQException):

    """
    Exception raised for errors related to the FrameReader class
    """
