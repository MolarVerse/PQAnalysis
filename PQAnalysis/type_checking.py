"""
A module for type checking of arguments passed to functions at runtime.
"""

import logging

from decorator import decorator
from beartype.door import is_bearable
from beartype.typing import ForwardRef

from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.exceptions import PQTypeError
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
def runtime_type_checking_setter(func, self, value):
    """
    A decorator to check the type of the arguments passed to a setter function at runtime.
    """

    type_hints = func.__annotations__

    # get var_name and type_hint from func.__annotations__
    var_name = list(type_hints.keys())[0]

    if isinstance(type_hints[var_name], str):
        type_hints[var_name] = ForwardRef(  # pylint: disable=protected-access
            type_hints[var_name]
        )._evaluate(globals(),
            locals(),
            frozenset())

    if not is_bearable(value, type_hints[var_name]):
        logger.error(
            get_type_error_message(
            var_name,
            value,
            type_hints[var_name],
            ),
            exception=PQTypeError,
        )

    # Call the function
    return func(self, value)



@decorator
def runtime_type_checking(func, *args, **kwargs):
    """
    A decorator to check the type of the arguments passed to a function at runtime.
    """

    if "disable_type_checking" in kwargs:
        disable_type_checking = kwargs.pop("disable_type_checking")
        if disable_type_checking:
            return func(*args, **kwargs)

    # Get the type hints of the function
    type_hints = func.__annotations__

    # Check the type of each argument
    for arg_name, arg_value in zip(func.__code__.co_varnames, args):
        if arg_name in type_hints:

            if isinstance(type_hints[arg_name], str):
                type_hints[arg_name] = ForwardRef(  # pylint: disable=protected-access
                    type_hints[arg_name]
                )._evaluate(globals(),
                    locals(),
                    frozenset())

            if not is_bearable(arg_value, type_hints[arg_name]):
                logger.error(
                    get_type_error_message(
                    arg_name,
                    arg_value,
                    type_hints[arg_name],
                    ),
                    exception=PQTypeError,
                )

    # Check the type of each keyword argument
    for kwarg_name, kwarg_value in kwargs.items():
        if kwarg_name in type_hints:
            if not is_bearable(kwarg_value, type_hints[kwarg_name]):
                logger.error(
                    get_type_error_message(
                    kwarg_name,
                    kwarg_value,
                    type_hints[kwarg_name],
                    ),
                    exception=TypeError,
                )

    # Call the function
    return func(*args, **kwargs)



def get_type_error_message(arg_name, value, expected_type):
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
