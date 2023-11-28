"""
A module containing different exceptions which could be useful.

...

Classes
-------
PQException
    Base class for exceptions in this module.
"""


class PQException(Exception):
    """
    Base class for exceptions in this package.
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)
