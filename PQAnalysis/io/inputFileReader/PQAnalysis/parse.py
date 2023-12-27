from beartype.typing import List

from ....types import PositiveInt, PositiveReal
from ..exceptions import InputFileError
from ..inputFileParser import InputDictionary


def _parse_positive_real(dict: InputDictionary, key: str) -> PositiveReal | None:
    try:
        data = dict[key]
    except KeyError:
        return None

    type = data[1]

    if type != "float":
        raise InputFileError(
            f"The \"{key}\" value has to be of float type - actual it is parsed as a {type}")

    if data[0] < 1:
        raise InputFileError(
            f"The \"{key}\" value has to be a positive floating point number!")

    return data[0]


def parse_traj_files_key(dict: InputDictionary) -> List[str] | None:
    try:
        data = dict["traj_files"]
    except KeyError:
        return None

    type = data[1]

    if type == "str":
        return [data[0]]
    elif type == "glob" or type == "list(str)":
        return data[0]
    else:
        raise InputFileError(
            f"The \"traj_files\" value has to be either a string, glob or a list of strings - actual it is parsed as a {type}")


def parse_n_bins_key(dict: InputDictionary) -> PositiveInt | None:
    try:
        data = dict["traj_files"]
    except KeyError:
        return None

    type = data[1]

    if type != "int":
        raise InputFileError(
            f"The \"n_bins\" value has to be of int type - actual it is parsed as a {type}")

    if data[0] < 1:
        raise InputFileError(
            "The \"n_bins\" value has to be a positive integer!")

    return data[0]


def parse_r_max_key(dict: InputDictionary) -> PositiveReal | None:
    return _parse_positive_real(dict, "r_max")


def parse_r_min_key(dict: InputDictionary) -> PositiveReal | None:
    return _parse_positive_real(dict, "r_min")
