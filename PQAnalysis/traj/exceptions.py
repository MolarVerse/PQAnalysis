"""
A module containing different exceptions related to the traj subpackage.
"""

from PQAnalysis.exceptions import PQException, BaseEnumFormatError



class TrajectoryFormatError(BaseEnumFormatError):

    """
    Exception raised if the given enum is not valid
    """



class MDEngineFormatError(BaseEnumFormatError):

    """
    Exception raised if the given enum is not valid
    """



class TrajectoryError(PQException):

    """
    Exception raised for errors related to the Trajectory class
    """

    def __init__(self, message: str) -> None:
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        self.message = message
        super().__init__(self.message)
