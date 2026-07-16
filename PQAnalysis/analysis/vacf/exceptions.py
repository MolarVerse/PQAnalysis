"""
A module containing different exceptions and warnings for
the :py:class:`~PQAnalysis.analysis.vacf.vacf.VACF` class.
"""

from ...exceptions import PQException, PQWarning



class VACFError(PQException):

    """
    Exception raised if something goes wrong during the VACF setup or calculation.
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



class VACFWarning(PQWarning):

    """
    Warning raised if something goes wrong during the VACF setup or calculation.
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
