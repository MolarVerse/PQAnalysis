"""
A module containing different exceptions related to the topology subpackage.

...

Classes
-------
MolTypeError
    Exception raised for errors related to the MolType class
"""

from ..exceptions import PQException


class MolTypeError(PQException):
    """
    Exception raised for errors related to the MolType class
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
