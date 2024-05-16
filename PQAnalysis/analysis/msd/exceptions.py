"""
A module containing different exceptions and warnings for 
the :py:class:`~PQAnalysis.analysis.msd.msd.MSD` class
"""

from ...exceptions import PQException, PQWarning


class MSDError(PQException):
    """
    Exception raised if something goes wrong during the MSD setup or calculation.
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


class MSDWarning(PQWarning):
    """
    Warning raised if something goes wrong during the MSD setup or calculation.
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
