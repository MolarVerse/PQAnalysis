"""
A module containing different exceptions and warnings for
the :py:class:`~PQAnalysis.analysis.adf.adf.ADF` class
"""

from ...exceptions import PQException, PQWarning



class ADFError(PQException):

    """
    Exception raised if something goes wrong during the ADF setup or calculation.
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



class ADFWarning(PQWarning):

    """
    Warning raised if something goes wrong during the ADF setup or calculation.
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
