"""
A module containing functions to parse the input file.
"""

import logging

from numbers import Real
from beartype.typing import List

from PQAnalysis.types import PositiveReal, PositiveInt
from PQAnalysis.io.input_file_reader.input_file_parser import InputDictionary
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis.exceptions import PQKeyError
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__

logger = logging.getLogger(__package_name__
                           ).getChild("PQAnalysisInputFileReader")
logger = setup_logger(logger)



def _parse_positive_real(
    input_dict: InputDictionary,
    key: str
) -> PositiveReal | None:
    """
    Gets the value of a key from the input dictionary and checks if it is a positive real number.
    If the key is not in the dictionary, None is returned.

    Parameters
    ----------
    dict : InputDictionary
        the input dictionary
    key : str
        the key to get the value from

    Returns
    -------
    PositiveReal | None
        the value of the key or None if the key is not in the dictionary

    Raises
    ------
    InputFileError
        if the value is not a positive real number
    """
    value = _parse_real(input_dict, key)

    if value is None:
        return None

    if value < 0:
        logger.error(
            f"The \"{key}\" value has to be a positive real number - It actually is {value}!",
            exception=InputFileError
        )

    return value



def _parse_real(input_dict: InputDictionary, key: str) -> Real | None:
    """
    Gets the value of a key from the input dictionary and checks if it is a real number.
    None is returned if the key is not in the dictionary.

    Parameters
    ----------
    dict : InputDictionary
        the input dictionary
    key : str
        the key to get the value from

    Returns
    -------
    PositiveReal | None
        the value of the key or None if the key is not in the dictionary

    Raises
    ------
    InputFileError
        if the value is not a real number
    """
    try:
        data = input_dict[key]
    except PQKeyError:
        return None

    if data[1] not in ["float", "int"]:
        logger.error(
            f"The \"{key}\" value has to be of float type - actually it is parsed as a {data[1]}",
            exception=InputFileError
        )

    return data[0]



def _parse_files(input_dict: InputDictionary, key: str) -> List[str] | None:
    """
    Gets the value of a key from the input dictionary and
    checks if it is a list of strings, a glob or a string.
    If the key is not in the dictionary, None is returned.

    Parameters
    ----------
    dict : InputDictionary
        the input dictionary
    key : str
        the key to get the value from

    Returns
    -------
    List[str] | None
        the value of the key or None if the key is not in the dictionary

    Raises
    ------
    InputFileError
        if the value is not a list of strings, a glob or a string
    """
    try:
        data = input_dict[key]
    except PQKeyError:
        return None

    data_type = data[1]

    if data_type == "str":
        return [data[0]]

    if data_type in {"glob", "list(str)"}:
        return data[0]

    logger.error(
        (
        f"The \"{key}\" value has to be either a "
        "string, glob or a list of strings - actually "
        f"it is parsed as a {data_type}"
        ),
        exception=InputFileError
    )



def _parse_int(input_dict: InputDictionary, key: str) -> int | None:
    """
    Gets the value of a key from the input dictionary and checks if it is an integer.

    Parameters
    ----------
    dict : InputDictionary
        the input dictionary
    key : str
        the key to get the value from

    Returns
    -------
    PositiveInt | None
        the value of the key or None if the key is not in the dictionary

    Raises
    ------
    InputFileError
        if the value is not an integer
    """
    try:
        data = input_dict[key]
    except PQKeyError:
        return None

    if data[1] != "int":
        logger.error(
            f"The \"{key}\" value has to be of int type - actually it is parsed as a {data[1]}",
            exception=InputFileError
        )

    return data[0]



def _parse_positive_int(
    input_dict: InputDictionary,
    key: str
) -> PositiveInt | None:
    """
    Gets the value of a key from the input dictionary and checks if it is a positive integer.

    Parameters
    ----------
    dict : InputDictionary
        the input dictionary
    key : str
        the key to get the value from

    Returns
    -------
    PositiveInt | None
        the value of the key or None if the key is not in the dictionary

    Raises
    ------
    InputFileError
        if the value is not a positive integer
    """
    value = _parse_int(input_dict, key)

    if value is None:
        return None

    if value < 1:
        logger.error(
            f"The \"{key}\" value has to be a positive integer - It actually is {value}!",
            exception=InputFileError
        )

    return value



def _parse_string(input_dict: InputDictionary, key: str) -> str | None:
    """
    Gets the value of a key from the input dictionary and checks if it is a string.

    Parameters
    ----------
    dict : InputDictionary
        the input dictionary
    key : str
        the key to get the value from

    Returns
    -------
    str | None
        the value of the key or None if the key is not in the dictionary

    Raises
    ------
    InputFileError
        if the value is not a string
    """
    try:
        data = input_dict[key]
    except PQKeyError:
        return None

    if data[1] != "str":
        logger.error(
            (
            f"The \"{key}\" value has to be of "
            f"string type - actually it is parsed as a {data[1]}"
            ),
            exception=InputFileError
        )

    return data[0]



def _parse_bool(input_dict: InputDictionary, key: str) -> bool | None:
    """
    Gets the value of a key from the input dictionary and checks if it is a bool.

    Parameters
    ----------
    dict : InputDictionary
        the input dictionary
    key : str
        the key to get the value from

    Returns
    -------
    bool | None
        the value of the key or None if the key is not in the dictionary

    Raises
    ------
    InputFileError
        if the value is not a bool
    """
    try:
        data = input_dict[key]
    except PQKeyError:
        return None

    if data[1] != "bool":
        logger.error(
            f"The \"{key}\" value has to be of bool type - actually it is parsed as a {data[1]}",
            exception=InputFileError
        )

    return data[0]
