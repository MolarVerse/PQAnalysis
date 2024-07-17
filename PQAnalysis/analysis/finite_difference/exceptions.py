"""
A module containing different exceptions and warnings for 
the :py:class:`~PQAnalysis.analysis.finite_difference.finite_difference.FiniteDifference` class
"""

from ...exceptions import PQException, PQWarning


class FiniteDifferenceError(PQException):

    """
    Exception raised if something goes wrong during the finite difference setup or calculation.
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


class FiniteDifferenceWarning(PQWarning):

    """
    Warning raised if something goes wrong during the finite difference setup or calculation.
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
