"""
A module containing different exceptions related to the core subpackage.

...

Classes
-------
ElementNotFoundError
    Exception raised if the given element id is not valid
AtomicSystemPositionsError
    Exception raised if atoms is not of the same length as positions
AtomicSystemMassError
    Exception raised if atoms do not contain mass information
AtomicSystemEmptySelectionWarning
    Warning raised if the selection is empty
"""

from PQAnalysis.exceptions import PQException, PQWarning


from beartype.typing import Any


class ElementNotFoundError(PQException):
    """
    Exception raised if the given element id is not valid
    """

    def __init__(self, id: Any) -> None:
        self.id = id
        self.message = f"""Id {self.id} is not a valid element identifier."""
        super().__init__(self.message)


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


class AtomicSystemEmptySelectionWarning(PQWarning):
    """
    Warning raised if the selection is empty
    """

    def __init__(self, message) -> None:
        self.message = message
        super().__init__(self.message)
