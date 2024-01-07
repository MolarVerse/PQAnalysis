"""
A module to read the input file for the radial distribution function.
"""

from ...io import PQAnalysisInputFileReader as Reader


class RDFInputFileReader(Reader):
    """
    A class to read the input file for the radial distribution function.

    The following keywords can be included: #TODO: make these keywords dependent on actual variables
        - traj_files
        - reference_selection
        - target_selection
        - out_file
        - r_max
        - r_min
        - delta_r
        - n_bins
        - use_full_atom_info
        - log_file
    """
    required_keys = [
        Reader.traj_files_key,
        Reader.reference_selection_key,
        Reader.target_selection_key,
        Reader.out_file_key
    ]

    known_keys = required_keys + [
        Reader.r_max_key,
        Reader.r_min_key,
        Reader.delta_r_key,
        Reader.n_bins_key,
        Reader.use_full_atom_info_key,
        Reader.log_file_key
    ]

    def __init__(self, filename: str) -> None:
        """
        Initialize the RDFInputFileReader.

        Parameters
        ----------
        filename : str
            the filename of the input file
        """
        self.filename = filename
        super().__init__(filename)

    def read(self):
        """
        Reads the input file and parses it.
        It also sets the raw_input_file and the dictionary.
        It checks if all required keys are set and if all keys are known.

        Raises
        ------
        InputFileError
            if not all required keys are set in the input file
        InputFileWarning
            if unknown keys are set in the input file 
        """
        super().read()
        super().check_required_keys(self.required_keys)
        super().check_known_keys(self.known_keys)

        if self.use_full_atom_info is None:
            self.use_full_atom_info = False
