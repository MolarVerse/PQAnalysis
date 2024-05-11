"""
A module containing different exceptions related to the core subpackage.
"""

from beartype.typing import Any

from PQAnalysis.exceptions import PQException, PQWarning



class ElementNotFoundError(PQException):

    """
    Exception raised if the given element id is not valid
    """

    def __init__(self, element_id: Any) -> None:
        """
        Parameters
        ----------
        element_id : Any
            The id that is not valid.
        """
        self.id = element_id
        self.message = f"""Id {self.id} is not a valid element identifier."""
        super().__init__(self.message)



class ResidueError(PQException):

    """
    Exception raised for errors related to the Residue class
    """

    def __init__(self, message: str):
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        self.message = message
        super().__init__(self.message)



class ResidueWarning(PQWarning):

    """
    Warning raised for problems related to the Residue class
    """

    def __init__(self, message: str):
        """
        Parameters
        ----------
        message : str
            The error message.
        """
        self.message = message
        super().__init__(self.message)



class AtomError(PQException):

    """
    Exception raised for errors related to the Atom class
    """
