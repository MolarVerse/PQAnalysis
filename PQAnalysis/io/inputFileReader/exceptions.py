from multimethod import multimethod

from ...exceptions import FormatEnumError


class InputFileFormatError(FormatEnumError):
    """
    Exception raised if the given enum is not valid
    """

    @multimethod
    def __init__(self, value: object, enum: object) -> None:
        super().__init__(value, enum)

    @multimethod
    def __init__(self, message: str) -> None:
        super().__init__(message)
