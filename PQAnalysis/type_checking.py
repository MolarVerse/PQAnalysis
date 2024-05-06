"""
A module for type checking of arguments passed to functions at runtime.
"""

import logging

from decorator import decorator
from beartype.door import is_bearable

from PQAnalysis.utils.custom_logging import setup_logger

logger_name = "PQAnalysis.TypeChecking"

if not logging.getLogger(logger_name).handlers:
    logger = setup_logger(logging.getLogger(logger_name))
else:
    logger = logging.getLogger(logger_name)


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
                        arg_value,
                        type_hints[arg_name],
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
                        kwarg_value,
                        type_hints[kwarg_name],
                    ),
                    exception=TypeError,
                )

    # Call the function
    return func(*args, **kwargs)


def _get_type_error_message(arg_name, value, expected_type):
    """
    Get the error message for a type error.
    """

    actual_type = type(value)

    header = (
        f"Argument '{arg_name}' with {value=} should be "
        f"of type {expected_type}, but got {actual_type}."
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
    elif str(expected_type) == '~SelectionCompatible':
        header += " Expected a SelectionCompatible object. For more "
        header += "information, see the documentation for Selections."

    return header
