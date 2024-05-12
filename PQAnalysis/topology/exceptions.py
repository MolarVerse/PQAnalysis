"""
A module containing different exceptions related to the topology subpackage.
"""

from PQAnalysis.exceptions import PQException



class TopologyError(PQException):

    """
    Exception raised for errors related to the Topology class
    """

    def __init__(self, message):
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        self.message = message
        super().__init__(self.message)
