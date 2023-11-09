"""
A module containing different exceptions which could be useful.
"""

from beartype.typing import Any


class ElementNotFoundError(Exception):
    """
    Exception raised if the given element id is not valid
    """

    def __init__(self, id: Any) -> None:
        self.id = id
        self.message = f"""Id {
            self.id} is not a valid element identifier."""
        super().__init__(self.message)
