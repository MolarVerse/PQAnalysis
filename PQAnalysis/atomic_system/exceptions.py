"""
A module containing different exceptions related to the core subpackage.
"""

from PQAnalysis.exceptions import PQException



class AtomicSystemError(PQException):

    """
    Exception raised for errors related to the AtomicSystem class
    """

    def __init__(self, message: str) -> None:
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        self.message = message
        super().__init__(self.message)



class AtomicSystemPositionsError(AtomicSystemError):

    """
    Exception raised if atoms is not of the same length as positions
    """

    message = "Atoms and positions must be of the same length."

    def __init__(self, message: str | None = None) -> None:
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        if message is not None:
            self.message = message

        super().__init__(self.message)



class AtomicSystemMassError(AtomicSystemError):

    """
    Exception raised if atoms do not contain mass information
    """

    message = "AtomicSystem contains atoms without mass information. "
    message += "Which is required for this operation."

    def __init__(self, message: str | None = None) -> None:
        """
        Parameters
        ----------
        message : str
            The error message.
        """

        if message is not None:
            self.message = message

        super().__init__(self.message)
