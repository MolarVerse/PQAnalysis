from ..exceptions import PQException


class RDFError(PQException):
    """
    Exception raised if the given element id is not valid
    """

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)
