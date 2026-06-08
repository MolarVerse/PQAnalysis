"""
A module containing a class to read input files to setup the 
:py:class:`~PQAnalysis.analysis.rdf.rdf.RDF` class.
"""
import logging
from pathlib import Path

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

    _moldescriptor_file_default = "moldescriptor.dat"

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
        super().not_defined_optional_keys(self.optional_keys)

        self._set_default_no_intra_molecular()

        if (
            self.no_intra_molecular is True and
            not self._infer_required_topology_files()
        ):
            self.logger.error(
                self._missing_topology_files_error_message(),
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

    def _infer_required_topology_files(self) -> bool:
        """
        Infers missing topology files required by no_intra_molecular.
        """
        if self.restart_file is None:
            restart_file = self._infer_restart_file()

            if restart_file is not None:
                self._set_inferred_key(
                    self.restart_file_key,
                    restart_file,
                    "str",
                )
                self.logger.info(
                    f"Inferred restart_file='{restart_file}' from traj_files."
                )

        if self.moldescriptor_file is None:
            moldescriptor_file = self._infer_moldescriptor_file()

            if moldescriptor_file is not None:
                self._set_inferred_key(
                    self.moldescriptor_file_key,
                    moldescriptor_file,
                    "str",
                )
                self.logger.info(
                    (
                        "Inferred moldescriptor_file="
                        f"'{moldescriptor_file}' from traj_files."
                    )
                )

        return (
            self.restart_file is not None and
            self.moldescriptor_file is not None
        )

    def _infer_restart_file(self) -> str | None:
        """
        Infers the restart file from the first trajectory filename.
        """
        restart_file = self._restart_file_candidate()

        if restart_file is None:
            return None

        return restart_file if Path(restart_file).is_file() else None

    def _infer_moldescriptor_file(self) -> str | None:
        """
        Infers the moldescriptor file from the first trajectory directory.
        """
        moldescriptor_file = self._moldescriptor_file_candidate()

        if moldescriptor_file is None:
            return None

        return moldescriptor_file if Path(moldescriptor_file).is_file() else None

    def _restart_file_candidate(self) -> str | None:
        """
        Gets the restart file candidate inferred from the first trajectory.
        """
        traj_file = self._first_traj_file()

        if traj_file is None:
            return None

        return str(Path(traj_file).with_suffix(".rst"))

    def _moldescriptor_file_candidate(self) -> str | None:
        """
        Gets the moldescriptor file candidate inferred from the first trajectory.
        """
        traj_file = self._first_traj_file()

        if traj_file is None:
            return None

        return str(Path(traj_file).parent / self._moldescriptor_file_default)

    def _first_traj_file(self) -> str | None:
        """
        Gets the first trajectory file from the input file.
        """
        traj_files = self.traj_files

        if not traj_files:
            return None

        return traj_files[0]

    def _missing_topology_files_error_message(self) -> str:
        """
        Gets the error message for missing topology files.
        """
        traj_file = self._first_traj_file()

        if traj_file is None:
            return (
                "The no_intra_molecular key requires both a restart file "
                "and a moldescriptor file. Could not infer missing files "
                "because no trajectory files were available."
            )

        return (
            "The no_intra_molecular key requires both a restart file "
            "and a moldescriptor file. Could not infer missing files "
            f"from trajectory file '{traj_file}'. "
            f"Tried restart_file='{self._restart_file_candidate()}' "
            "and moldescriptor_file="
            f"'{self._moldescriptor_file_candidate()}'."
        )

    def _set_inferred_key(self, key: str, value, value_type: str):
        """
        Sets a value in the input dictionary after it was inferred.
        """
        self.dictionary.dict[key] = (value, value_type, "inferred")



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
