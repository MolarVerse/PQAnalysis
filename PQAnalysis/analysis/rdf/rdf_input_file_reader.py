"""
A module containing a class to read input files to setup the 
:py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` class.
"""
import logging

# local imports
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import PQAnalysisInputFileReader as Reader
from PQAnalysis.io.input_file_reader.pq_analysis._topology_files_mixin import (
    _TopologyFilesMixin,
)
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking



class RDFInputFileReader(_TopologyFilesMixin, Reader):

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
        Reader.reference_com_key,
        Reader.target_com_key,
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
            if the reference_com or target_com key is used without a restart
            file and a moldescriptor file
        InputFileError
            if the moldescriptor_file key is used without a restart file
        """
        super().read()
        super().check_required_keys(self.required_keys)
        super().check_known_keys(self.required_keys + self.optional_keys)
        super().not_defined_optional_keys(self.optional_keys)

        self._set_default_no_intra_molecular()

        if self.no_intra_molecular is True:
            self._require_topology_files(self.no_intra_molecular_key)

        if self.reference_com is True or self.target_com is True:
            self._require_topology_files(
                f"{self.reference_com_key} or {self.target_com_key}"
            )

        self._check_moldescriptor_needs_restart()

    def _set_default_no_intra_molecular(self):
        """
        Enables no_intra_molecular if both topology files are provided
        and the keyword was omitted.
        """
        if self.no_intra_molecular is not None:
            return

        if self.restart_file is None or self.moldescriptor_file is None:
            return

        self._set_inferred_key(self.no_intra_molecular_key, True, "bool")
        self.logger.info(
            (
                "Enabled no_intra_molecular because both restart_file "
                "and moldescriptor_file are set."
            )
        )




input_keys_documentation = f"""

For the RDF analysis input file several keys are available of which some are required and some are optional. For more details on the grammar and syntax of the input file see :ref:`inputFile`.

The following keys are required:

.. list-table::
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
        - The output file to write the RDF data to; see :ref:`RDF output files <analysis-output-rdf>` for its columns, units and normalization formulas.

The following keys are optional:

.. list-table::
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
    * - {Reader.reference_com_key}
        - Whether the reference entities are the centers of mass of the molecules instead of atoms; the reference selection then picks the atoms that define each center. Needs the restart file and the moldescriptor file.
    * - {Reader.target_com_key}
        - Whether the target entities are the centers of mass of the molecules instead of atoms, see :code:`{Reader.reference_com_key}`.
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

- Basic RDF calculations can be run without a restart file. In that
  case, intra molecular pairs are included.
- If both :code:`{Reader.restart_file_key}` and
  :code:`{Reader.moldescriptor_file_key}` are given and
  :code:`{Reader.no_intra_molecular_key}` is omitted,
  :code:`{Reader.no_intra_molecular_key}` defaults to :code:`True`.
- If the :code:`{Reader.no_intra_molecular_key}` key is enabled,
  the :code:`{Reader.restart_file_key}` and
  :code:`{Reader.moldescriptor_file_key}` keys
  are required in order to exclude intra molecular pairs. Missing
  files are inferred from the first trajectory file if possible:
  :code:`trajectory.xyz` maps to :code:`trajectory.rst`, and
  :code:`moldescriptor.dat` is searched in the trajectory directory.
- If :code:`{Reader.reference_com_key}` or :code:`{Reader.target_com_key}`
  is enabled, both files are required as well and are inferred in the same
  way. With centers on both sides, :code:`{Reader.no_intra_molecular_key}`
  has no effect. With centers on one side it removes the pairs of an atom
  with the center of its own molecule and, as for atoms, defaults to
  :code:`True` when both files are given.
- If the :code:`{Reader.moldescriptor_file_key}` key is specified,
  the :code:`{Reader.restart_file_key}` key is required in order
  to use the reference residues in any meaningful way.
- Explicitly provided :code:`{Reader.restart_file_key}` and
  :code:`{Reader.moldescriptor_file_key}` values are used as-is.
- Inferred and defaulted values are written to the normal PQAnalysis
  log output.
- In general, the :code:`{Reader.r_max_key}`,
  :code:`{Reader.n_bins_key}` and :code:`{Reader.delta_r_key}`
  are mutual exclusive, meaning that they can't be specified at
  the same time. Furthermore, at least one of
  :code:`{Reader.n_bins_key}` or :code:`{Reader.delta_r_key}` is required
  (for more information see :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF`).

"""

RDFInputFileReader.__doc__ += input_keys_documentation
