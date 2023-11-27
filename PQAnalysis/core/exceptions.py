"""
A module containing different exceptions related to the core subpackage.

...

Classes
-------
ElementNotFoundError
    Exception raised if the given element id is not valid
"""

from PQAnalysis.exceptions import PQException


from beartype.typing import Any


class ElementNotFoundError(PQException):
    """
    Exception raised if the given element id is not valid
    """

    def __init__(self, id: Any) -> None:
        self.id = id
        self.message = f"""Id {self.id} is not a valid element identifier."""
        super().__init__(self.message)
