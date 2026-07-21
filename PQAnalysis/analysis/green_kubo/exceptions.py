"""
A module containing different exceptions and warnings for
the :py:class:`~PQAnalysis.analysis.green_kubo.green_kubo.GreenKubo`
class.
"""

from ...exceptions import PQException, PQWarning



class GreenKuboError(PQException):

    """
    Exception raised if something goes wrong during the Green-Kubo
    setup or calculation.
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



class GreenKuboWarning(PQWarning):

    """
    Warning raised if something goes wrong during the Green-Kubo
    setup or calculation.
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
