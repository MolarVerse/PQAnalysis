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
                    f"Argument '{arg_name}' should be of type "
                    f"{type_hints[arg_name]}, but got {type(arg_value)}",
                    exception=TypeError,
                )

    # Check the type of each keyword argument
    for kwarg_name, kwarg_value in kwargs.items():
        if kwarg_name in type_hints:
            if not is_bearable(kwarg_value, type_hints[kwarg_name]):
                logger.error(
                    f"Argument '{kwarg_name}' should be of type "
                    f"{type_hints[kwarg_name]}, but got {type(kwarg_value)}",
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
        f"but got {actual_type}"
    )
