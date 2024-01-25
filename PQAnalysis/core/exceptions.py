"""
A module containing different exceptions related to the core subpackage.
"""

from beartype.typing import Any

from PQAnalysis.exceptions import PQException, PQWarning


class ElementNotFoundError(PQException):
    """
    Exception raised if the given element id is not valid
    """

    def __init__(self, id: Any) -> None:
        self.id = id
        self.message = f"""Id {self.id} is not a valid element identifier."""
        super().__init__(self.message)


class ResidueError(PQException):
    """
    Exception raised for errors related to the Residue class
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class ResidueWarning(PQWarning):
    """
    Warning raised for problems related to the Residue class
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
