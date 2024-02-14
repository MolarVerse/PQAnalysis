"""
A module containing different exceptions related to the core subpackage.
"""

from PQAnalysis.exceptions import PQException


class AtomicSystemPositionsError(PQException):
    """
    Exception raised if atoms is not of the same length as positions
    """

    message = """Atoms and positions must be of the same length."""

    def __init__(self) -> None:
        super().__init__(self.message)


class AtomicSystemMassError(PQException):
    """
    Exception raised if atoms do not contain mass information
    """

    message = """AtomicSystem contains atoms without mass information. Which is required for this operation."""

    def __init__(self) -> None:
        super().__init__(self.message)
