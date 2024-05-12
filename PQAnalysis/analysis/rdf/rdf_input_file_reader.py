"""
A module containing a class to read input files to setup the 
:py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` class.
"""
import logging

# local imports
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import PQAnalysisInputFileReader as Reader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking



class RDFInputFileReader(Reader):

    """
    A class to read input files to setup the :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` class.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

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

    @runtime_type_checking
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
        InputFileError
            if the no_intra_molecular key is used without a restart file
        InputFileError
            if the no_intra_molecular key is used without a moldescriptor file
        InputFileError
            if the moldescriptor_file key is used without a restart file
        """
        super().read()
        super().check_required_keys(self.required_keys)
        super().check_known_keys(self.required_keys + self.optional_keys)

        if (self.no_intra_molecular is not None and
            (self.restart_file is None or self.moldescriptor_file is None)):
            self.logger.error(
                (
                "The no_intra_molecular key can only be used "
                "if both a restart file and a moldescriptor file are given."
                ),
                exception=InputFileError,
            )

        if self.moldescriptor_file is not None and self.restart_file is None:
            self.logger.error(
                (
                "The moldescriptor_file key can only be "
                "used in a meaningful way if a restart file is given."
                ),
                exception=InputFileError,
            )



input_keys_documentation = f"""

For the RDF analysis input file several keys are available of which some are required and some are optional. For more details on the grammar and syntax of the input file see :ref:`inputFile`.

.. list-table:: Required keys
    :header-rows: 1

    * - Key
        - Value
    * - {Reader.traj_files_key}
        - The trajectory files to read. This can be a single file or a list of files.
    * - {Reader.reference_selection_key}
        - The selection string to select the reference atoms. For more details see :py:class:`~PQAnalysis.topology.selection.Selection`.
    * - {Reader.target_selection_key}
        - The selection string to select the target atoms. For more details see :py:class:`~PQAnalysis.topology.selection.Selection`.
    * - {Reader.out_file_key}
        - The output file to write the RDF data to.

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
        - The log file to write the log information to.
    * - {Reader.use_full_atom_info_key}
        - Whether to use full atom information for the selections.
    * - {Reader.no_intra_molecular_key}
        - Whether to exclude intra molecular pairs.
    * - {Reader.restart_file_key}
        - The restart file to read the topology from.
    * - {Reader.moldescriptor_file_key}
        - The moldescriptor file to read the reference residues from.

Note
----
Optional keys does not mean that they are optional for the analysis.
They are optional in the input file, but they might be required for
the analysis. This means that if an optional keyword is specified
other keywords might be required. For example:

- If the :code:`{Reader.no_intra_molecular_key}` key is specified,
the :code:`{Reader.restart_file_key}` and
:code:`{Reader.moldescriptor_file_key}` keys
are required in order to exclude intra molecular pairs.
- If the :code:`{Reader.moldescriptor_file_key}` key is specified,
the :code:`{Reader.restart_file_key}` key is required in order
use the reference residues in any meaningful way.
- In general, the :code:`{Reader.r_max_key}`,
:code:`{Reader.n_bins_key}` and :code:`{Reader.delta_r_key}`
are mutual exclusive, meaning that they can't be specified at
the same time. Furthermore, at least one of
:code:`{Reader.n_bins}' or :code:`{Reader.delta_r}` is required
(for more information see :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF`).

"""

RDFInputFileReader.__doc__ += input_keys_documentation
