"""
A module for reading the input file for the PQAnalysis.
"""
import logging

from beartype.typing import List

from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis.io.input_file_reader.input_file_parser import InputFileParser
from PQAnalysis.io.input_file_reader.formats import InputFileFormat
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__

from ._file_mixin import _FileMixin
from ._selection_mixin import _SelectionMixin
from ._positions_mixin import _PositionsMixin



class PQAnalysisInputFileReader(_FileMixin, _SelectionMixin, _PositionsMixin):

    """
    A class to read the input file for the PQAnalysis. 
    It inherits from _FileMixin, _SelectionMixin and _PositionsMixin.

    The main idea of this class is work as a super class for all
    command line tools of PQAnalysis, but it can also be used 
    as a standalone class. But if for example an analysis tool
    needs only few input file keywords, it is recommended to 
    inherit from this class.

    Attention:
        As this class is not intended to be used as a standalone class,
        it does not check if all required keys are set in the input file,
        thus also no information/documentation of the values of the 
        keywords is given, so that they can be used accordingly in
        the implementing subclasses.
    """
    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    known_keys = []

    # fmt: off
    r_max_key = "r_max"; known_keys.append(r_max_key) # pylint: disable=multiple-statements
    r_min_key = "r_min"; known_keys.append(r_min_key) # pylint: disable=multiple-statements
    delta_r_key = "delta_r"; known_keys.append(delta_r_key) # pylint: disable=multiple-statements
    n_bins_key = "n_bins"; known_keys.append(n_bins_key) # pylint: disable=multiple-statements
    traj_files_key = "traj_files"; known_keys.append(traj_files_key) # pylint: disable=multiple-statements
    selection_key = "selection"; known_keys.append(selection_key) # pylint: disable=multiple-statements
    reference_selection_key = "reference_selection"; known_keys.append(reference_selection_key) # pylint: disable=multiple-statements
    target_selection_key = "target_selection"; known_keys.append(target_selection_key) # pylint: disable=multiple-statements
    use_full_atom_info_key = "use_full_atom_info"; known_keys.append(use_full_atom_info_key) # pylint: disable=multiple-statements
    out_file_key = "out_file"; known_keys.append(out_file_key) # pylint: disable=multiple-statements
    log_file_key = "log_file"; known_keys.append(log_file_key) # pylint: disable=multiple-statements
    no_intra_molecular_key = "no_intra_molecular"; known_keys.append(no_intra_molecular_key) # pylint: disable=multiple-statements
    moldescriptor_file_key = "moldescriptor_file"; known_keys.append(moldescriptor_file_key) # pylint: disable=multiple-statements
    restart_file_key = "restart_file"; known_keys.append(restart_file_key) # pylint: disable=multiple-statements
    #fmt: on

    def __init__(self, filename: str) -> None:
        """
        It sets the format to InputFileFormat.PQANALYSIS and the
        filename to the given filename. It also creates a 
        InputFileParser with the given filename.

        Parameters
        ----------
        filename : str
            the filename of the input file
        """
        self.format = InputFileFormat.PQANALYSIS
        self.filename = filename
        self.parser = InputFileParser(filename)

        self.dictionary = None
        self.raw_input_file = None

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
        if not all(key in self.dictionary.keys() for key in required_keys):
            self.logger.error(
                (
                "Not all required keys were set in "
                f"the input file! The required keys are: {required_keys}."
                ),
                exception=InputFileError
            )

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

        if not all(key in known_keys for key in self.dictionary.keys()):
            self.logger.warning(
                "Unknown keys were set in the input file! "
                f"The known keys are: {known_keys}. They will be ignored!"
            )
