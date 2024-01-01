"""
A module containing different exceptions related to the topology subpackage.

...

Classes
-------
ResidueError
    Exception raised for errors related to the Residue class
"""

from ..exceptions import PQException


class ResidueError(PQException):
    """
    Exception raised for errors related to the Residue class
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
