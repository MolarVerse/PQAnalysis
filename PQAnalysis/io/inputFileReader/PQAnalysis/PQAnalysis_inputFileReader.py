import warnings

from beartype.typing import List

from ....types import PositiveInt, PositiveReal
from ..exceptions import InputFileError, InputFileWarning
from ..inputFileParser import InputDictionary, InputFileParser
from ..formats import InputFileFormat


class PQAnalysisInputFileReader:

    known_keys = []

    # fmt: off
    r_max_key = "r_max"; known_keys.append(r_max_key)
    r_min_key = "r_min"; known_keys.append(r_min_key)
    delta_r_key = "delta_r"; known_keys.append(delta_r_key)
    n_bins_key = "n_bins"; known_keys.append(n_bins_key)
    traj_files_key = "traj_files"; known_keys.append(traj_files_key)
    selection_key = "selection"; known_keys.append(selection_key)
    reference_selection_key = "reference_selection"; known_keys.append(reference_selection_key)
    target_selection_key = "target_selection"; known_keys.append(target_selection_key)
    use_full_atom_info_key = "use_full_atom_info"; known_keys.append(use_full_atom_info_key)
    out_file_key = "out_file"; known_keys.append(out_file_key)
    log_file_key = "log_file"; known_keys.append(log_file_key)
    #fmt: on

    def __init__(self, filename: str) -> None:
        self.format = InputFileFormat.PQANALYSIS
        self.filename = filename
        self.parser = InputFileParser(filename)

    def read(self):
        self.dictionary = self.parser.parse()
        self.raw_input_file = self.parser.raw_input_file

    def check_required_keys(self, required_keys: List[str]):
        if not all([key in self.dictionary.keys() for key in required_keys]):
            raise InputFileError(
                f"Not all required keys were set in the input file! The required keys are: {required_keys}.")

    def check_known_keys(self, known_keys: List[str] | None = None):
        if known_keys is None:
            known_keys = self.known_keys

        if not all([key in known_keys for key in self.dictionary.keys()]):
            warnings.warn(
                f"Unknown keys were set in the input file! The known keys are: {known_keys}. They will be ignored!", InputFileWarning)

    @property
    def r_max(self) -> PositiveReal | None:
        return _parse_positive_real(self.dictionary, self.r_max_key)

    @property
    def r_min(self) -> PositiveReal | None:
        return _parse_positive_real(self.dictionary, self.r_min_key)

    @property
    def delta_r(self) -> PositiveReal | None:
        return _parse_positive_real(self.dictionary, self.delta_r_key)

    @property
    def traj_files(self) -> List[str] | None:
        return _parse_files(self.dictionary, self.traj_files_key)

    @property
    def out_file(self) -> str | None:
        return _parse_string(self.dictionary, self.out_file_key)

    @property
    def log_file(self) -> str | None:
        return _parse_string(self.dictionary, self.log_file_key)

    @property
    def n_bins(self) -> PositiveInt | None:
        return _parse_positive_int(self.dictionary, self.n_bins_key)

    @property
    def selection(self) -> str | None:
        return _parse_string(self.dictionary, self.selection_key)

    @property
    def reference_selection(self) -> str | None:
        return _parse_string(self.dictionary, self.reference_selection_key)

    @property
    def target_selection(self) -> str | None:
        return _parse_string(self.dictionary, self.target_selection_key)

    @property
    def use_full_atom_info_for_selection(self) -> bool | None:
        return _parse_bool(self.dictionary, self.use_full_atom_info_key)


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


def _parse_files(dict: InputDictionary, key: str) -> List[str] | None:
    try:
        data = dict[key]
    except KeyError:
        return None

    type = data[1]

    if type == "str":
        return [data[0]]
    elif type == "glob" or type == "list(str)":
        return data[0]
    else:
        raise InputFileError(
            f"The \"{key}\" value has to be either a string, glob or a list of strings - actual it is parsed as a {type}")


def _parse_int(dict: InputDictionary, key: str) -> PositiveInt | None:
    try:
        data = dict[key]
    except KeyError:
        return None

    type = data[1]

    if type != "int":
        raise InputFileError(
            f"The \"{key}\" value has to be of int type - actual it is parsed as a {type}")

    return data[0]


def _parse_positive_int(dict: InputDictionary, key: str) -> PositiveInt | None:

    value = _parse_int(dict, key)

    if value < 1:
        raise InputFileError(
            "The \"{key}\" value has to be a positive integer - It actually is {value}!")

    return value


def _parse_string(dict: InputDictionary, key: str) -> str | None:
    try:
        data = dict[key]
    except KeyError:
        return None

    type = data[1]

    if type != "str":
        raise InputFileError(
            f"The \"{key}\" value has to be of string type - actual it is parsed as a {type}")

    return data[0]


def _parse_bool(dict: InputDictionary, key: str) -> bool | None:
    try:
        data = dict[key]
    except KeyError:
        return None

    type = data[1]

    if type != "bool":
        raise InputFileError(
            f"The \"{key}\" value has to be of bool type - actual it is parsed as a {type}")

    return data[0]
