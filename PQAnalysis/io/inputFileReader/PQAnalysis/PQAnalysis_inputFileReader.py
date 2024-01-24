"""
A module for reading the input file for the PQAnalysis.
"""
import warnings

from beartype.typing import List

from ._fileMixin import _FileMixin
from ._selectionMixin import _SelectionMixin
from ._positionsMixin import _PositionsMixin
from ..exceptions import InputFileError, InputFileWarning
from ..inputFileParser import InputFileParser
from ..formats import InputFileFormat


class PQAnalysisInputFileReader(_FileMixin, _SelectionMixin, _PositionsMixin):
    """
    A class to read the input file for the PQAnalysis. It inherits from _FileMixin, _SelectionMixin and _PositionsMixin.

    The main idea of this class is work as a super class for all command line tools of PQAnalysis, but it can also be used as a standalone class.
    But if for example an analysis tool needs only few input file keywords, it is recommended to inherit from this class.

    Attention:
        As this class is not intended to be used as a standalone class, it does not check if all required keys are set in the input file,
        thus also no information/documentation of the values of the keywords is given, so that they can be used accordingly in the implementing
        subclasses.
    """
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
    no_intra_molecular_key = "no_intra_molecular"; known_keys.append(no_intra_molecular_key)
    moldescriptor_file_key = "moldescriptor_file"; known_keys.append(moldescriptor_file_key)
    restart_file_key = "restart_file"; known_keys.append(restart_file_key)
    #fmt: on

    def __init__(self, filename: str) -> None:
        """
        It sets the format to InputFileFormat.PQANALYSIS and the filename to the given filename. It also creates a InputFileParser with the given filename.

        Parameters
        ----------
        filename : str
            the filename of the input file
        """
        self.format = InputFileFormat.PQANALYSIS
        self.filename = filename
        self.parser = InputFileParser(filename)

    def read(self):
        """
        Reads the input file and parses it.
        It also sets the raw_input_file and the dictionary.
        """
        self.dictionary = self.parser.parse()
        self.raw_input_file = self.parser.raw_input_file

    def check_required_keys(self, required_keys: List[str]):
        """
        Checks if all required keys are set in the input file.

        Parameters
        ----------
        required_keys : List[str]
            the required keys

        Raises
        ------
        InputFileError
            if not all required keys are set in the input file
        """
        if not all([key in self.dictionary.keys() for key in required_keys]):
            raise InputFileError(
                f"Not all required keys were set in the input file! The required keys are: {required_keys}.")

    def check_known_keys(self, known_keys: List[str] | None = None):
        """
        Checks if all keys in the input file are known.

        Parameters
        ----------
        known_keys : List[str] | None, optional
            the known keys, by default None

        Warns
        -----
        InputFileWarning
            if unknown keys were set in the input file
        """
        if known_keys is None:
            known_keys = self.known_keys

        if not all([key in known_keys for key in self.dictionary.keys()]):
            warnings.warn(
                f"Unknown keys were set in the input file! The known keys are: {known_keys}. They will be ignored!", InputFileWarning)
