"""
A module containing different exceptions related to the frame subpackage.
"""

from PQAnalysis.exceptions import PQException, BaseEnumFormatError


class FrameError(PQException):
    """
    Exception raised for errors related to the Frame class
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)
