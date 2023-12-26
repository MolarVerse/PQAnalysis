"""
A module containing different exceptions for the RDF class.

...

Classes
-------
RDFError
    Exception raised if the given element id is not valid
RDFWarning
    Warning raised if the given element id is not valid
"""

from ...exceptions import PQException, PQWarning


class RDFError(PQException):
    """
    Exception raised if the given element id is not valid
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)


class RDFWarning(PQWarning):
    """
    Warning raised if the given element id is not valid
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)
