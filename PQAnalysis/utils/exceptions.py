"""
A module containing different exceptions which could be useful.
"""


class ElementNotFoundError(Exception):
    """
    Exception raised if the given element id is not valid
    """

    def __init__(self, id):
        self.id = id
        self.message = f"""Element with id {
            self.id} is not a valid element identifier."""
        super().__init__(self.message)
