"""
A module containing different exceptions and warnings for 
the :py:class:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus.BulkModulus` class
"""

from ...exceptions import PQException, PQWarning


class BulkModulusError(PQException):

    """
    Exception raised if something goes wrong during the BulkModulus setup or calculation.
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


class BulkModulusWarning(PQWarning):

    """
    Warning raised if something goes wrong during the BulkModulus setup or calculation.
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
