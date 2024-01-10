"""
A module containing a class to read input files to setup the :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` class.
"""

from ...io import PQAnalysisInputFileReader as Reader
from ...utils.docs import extend_documentation


def fun():
    return """
hello
"""


class RDFInputFileReader(Reader):
    """
    A class to read input files to setup the :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` class.
    """

    #: List[str]: The required keys of the input file
    required_keys = [
        Reader.traj_files_key,
        Reader.reference_selection_key,
        Reader.target_selection_key,
        Reader.out_file_key,
        Reader.restart_file_key,
        Reader.moldescriptor_file_key,
    ]

    #: List[str]: The optional keys of the input file
    optional_keys = required_keys + [
        Reader.r_max_key,
        Reader.r_min_key,
        Reader.delta_r_key,
        Reader.n_bins_key,
        Reader.use_full_atom_info_key,
        Reader.no_intra_molecular_key,
        Reader.log_file_key
    ]

    def __init__(self, filename: str) -> None:
        """
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
        super().check_known_keys(self.required_keys + self.optional_keys)
