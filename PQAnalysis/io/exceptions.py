"""
A module containing different exceptions related to the io subpackage.

...

Classes
-------
BoxWriterError
    Exception raised for errors related to the BoxWriter class
"""

from ..exceptions import PQException


class BoxWriterError(PQException):
    """
    Exception raised for errors related to the BoxWriter class
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)
