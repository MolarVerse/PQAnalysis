import logging

from decorator import decorator
from beartype.door import is_bearable

from PQAnalysis.utils.custom_logging import setup_logger
from .types import (
    Np1DIntArray,
    Np2DIntArray,
    Np1DNumberArray,
    Np2DNumberArray,
    Np3x3NumberArray,
    NpnDNumberArray,
)

__logger_name__ = "PQAnalysis.TypeChecking"

if not logging.getLogger(__logger_name__).handlers:
    logger = setup_logger(logging.getLogger(__logger_name__))
else:
    logger = logging.getLogger(__logger_name__)


@decorator
def runtime_type_checking(func, *args, **kwargs):
    """
    A decorator to check the type of the arguments passed to a function at runtime.
    """

    # Get the type hints of the function
    type_hints = func.__annotations__

    # Check the type of each argument
    for arg_name, arg_value in zip(func.__code__.co_varnames, args):
        if arg_name in type_hints:
            if not is_bearable(arg_value, type_hints[arg_name]):
                logger.error(
                    _get_type_error_message(
                        arg_name,
                        type_hints[arg_name],
                        type(arg_value)
                    ),
                    exception=TypeError,
                )

    # Check the type of each keyword argument
    for kwarg_name, kwarg_value in kwargs.items():
        if kwarg_name in type_hints:
            if not is_bearable(kwarg_value, type_hints[kwarg_name]):
                logger.error(
                    _get_type_error_message(
                        kwarg_name,
                        type_hints[kwarg_name],
                        type(kwarg_value)
                    ),
                    exception=TypeError,
                )

    # Call the function
    return func(*args, **kwargs)


def _get_type_error_message(arg_name, expected_type, actual_type):
    """
    Get the error message for a type error.
    """

    header = (
        f"Argument '{arg_name}' should be of type {expected_type}, "
        f"but got {actual_type}."
    )

    if expected_type is Np1DIntArray:
        header += " Expected a 1D numpy integer array."
    elif expected_type is Np2DIntArray:
        header += " Expected a 2D numpy integer array."
    elif expected_type is Np1DNumberArray:
        header += " Expected a 1D numpy number array."
    elif expected_type is Np2DNumberArray:
        header += " Expected a 2D numpy number array."
    elif expected_type is Np3x3NumberArray:
        header += " Expected a 3x3 numpy number array."
    elif expected_type is NpnDNumberArray:
        header += " Expected an n-dimensional numpy number array."

    return header
