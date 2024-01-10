"""
A module containing a class to read input files to setup the :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` class.
"""

from __future__ import annotations

from ...io import PQAnalysisInputFileReader as Reader


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
    ]

    #: List[str]: The optional keys of the input file
    optional_keys = required_keys + [
        Reader.r_max_key,
        Reader.r_min_key,
        Reader.delta_r_key,
        Reader.n_bins_key,
        Reader.log_file_key,
        Reader.use_full_atom_info_key,
        Reader.no_intra_molecular_key,
        Reader.restart_file_key,
        Reader.moldescriptor_file_key,
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


input_keys_documentation = f"""

For the RDF analysis input file several keys are available of which some are required and some are optional.

.. list-table:: Required keys
    :header-rows: 1

    * - Key
        - Value
    * - {Reader.traj_files_key}
        - The trajectory files to read.
    * - {Reader.reference_selection_key}
        - The selection string to select the reference atoms.
    * - {Reader.target_selection_key}
        - The selection string to select the target atoms.
    * - {Reader.out_file_key}
        - The output file to write the RDF to.

.. list-table:: Optional keys
    :header-rows: 1

    * - Key
        - Value
    * - {Reader.r_max_key}
        - The maximum radius to calculate the RDF for.
    * - {Reader.r_min_key}
        - The minimum radius to calculate the RDF for.
    * - {Reader.delta_r_key}
        - The delta radius to calculate the RDF for.
    * - {Reader.n_bins_key}
        - The number of bins to calculate the RDF for.
    * - {Reader.log_file_key}
        - The log file to write the log to.
    * - {Reader.use_full_atom_info_key}
        - Whether to use full atom information.
    * - {Reader.no_intra_molecular_key}
        - Whether to exclude intra molecular pairs.
    * - {Reader.restart_file_key}
        - The restart file to write the restart information to.
    * - {Reader.moldescriptor_file_key}
        - The moldescriptor file to write the moldescriptor to.

"""

RDFInputFileReader.__doc__ += input_keys_documentation
