"""
A module containing different exceptions and warnings for 
the :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` class
"""

from ...exceptions import PQException, PQWarning



class RDFError(PQException):

    """
    Exception raised if something goes wrong during the RDF setup or calculation.
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



class RDFWarning(PQWarning):

    """
    Warning raised if something goes wrong during the RDF setup or calculation.
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
