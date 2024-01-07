"""
A module containing different exceptions and warnings for the RDF class.

"""

from ...exceptions import PQException, PQWarning


class RDFError(PQException):
    """
    Exception raised if something goes wrong during the RDF calculation.
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)


class RDFWarning(PQWarning):
    """
    Warning raised if something goes wrong during the RDF calculation.
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)
