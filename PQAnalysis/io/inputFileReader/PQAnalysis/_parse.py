from beartype.typing import List

import PQAnalysis.io.inputFileReader.inputFileParser as inputFileParser
from ....types import PositiveReal, PositiveInt


def _parse_positive_real(dict: inputFileParser.InputDictionary, key: str) -> PositiveReal | None:
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
    value = _parse_real(dict, key)

    if value < 0:
        raise inputFileParser.inputFileParser.InputFileError(
            "The \"{key}\" value has to be a positive real number - It actually is {value}!")

    return value


def _parse_real(dict: inputFileParser.InputDictionary, key: str) -> PositiveReal | None:
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
        data = dict[key]
    except KeyError:
        return None

    data_type = data[1]

    if data_type != "float" or data_type != "int":
        raise inputFileParser.inputFileParser.InputFileError(
            f"The \"{key}\" value has to be of float type - actual it is parsed as a {data_type}")

    return data[0]


def _parse_files(dict: inputFileParser.InputDictionary, key: str) -> List[str] | None:
    """
    Gets the value of a key from the input dictionary and checks if it is a list of strings, a glob or a string.
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
        data = dict[key]
    except KeyError:
        return None

    data_type = data[1]

    if data_type == "str":
        return [data[0]]
    elif data_type == "glob" or data_type == "list(str)":
        return data[0]
    else:
        raise inputFileParser.InputFileError(
            f"The \"{key}\" value has to be either a string, glob or a list of strings - actual it is parsed as a {data_type}")


def _parse_int(dict: inputFileParser.InputDictionary, key: str) -> PositiveInt | None:
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
        data = dict[key]
    except KeyError:
        return None

    data_type = data[1]

    if data_type != "int":
        raise inputFileParser.InputFileError(
            f"The \"{key}\" value has to be of int type - actual it is parsed as a {data_type}")

    return data[0]


def _parse_positive_int(dict: inputFileParser.InputDictionary, key: str) -> PositiveInt | None:
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
    value = _parse_int(dict, key)

    if value < 1:
        raise inputFileParser.InputFileError(
            "The \"{key}\" value has to be a positive integer - It actually is {value}!")

    return value


def _parse_string(dict: inputFileParser.InputDictionary, key: str) -> str | None:
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
        data = dict[key]
    except KeyError:
        return None

    data_type = data[1]

    if data_type != "str":
        raise inputFileParser.InputFileError(
            f"The \"{key}\" value has to be of string type - actual it is parsed as a {data_type}")

    return data[0]


def _parse_bool(dict: inputFileParser.InputDictionary, key: str) -> bool | None:
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
        data = dict[key]
    except KeyError:
        return None

    data_type = data[1]

    if data_type != "bool":
        raise inputFileParser.InputFileError(
            f"The \"{key}\" value has to be of bool type - actual it is parsed as a {data_type}")

    return data[0]
