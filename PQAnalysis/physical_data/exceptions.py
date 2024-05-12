"""
A module containing different exceptions related to the physicalData subpackage.
"""

from PQAnalysis.exceptions import PQException



class EnergyError(PQException):

    """
    Exception raised for errors related to the Energy class
    """

    def __init__(self, message: str):
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        self.message = message
        super().__init__(self.message)
