"""
A module containing a class to read input files to setup the
:py:class:`~PQAnalysis.analysis.adf.adf.ADF` class.
"""
import logging

# local imports
from PQAnalysis.types import PositiveInt, PositiveReal
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import PQAnalysisInputFileReader as Reader
from PQAnalysis.io.input_file_reader.pq_analysis._parse import (
    _parse_string,
    _parse_positive_int,
    _parse_positive_real,
)
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking



class ADFInputFileReader(Reader):

    """
    A class to read input files to setup the
    :py:class:`~PQAnalysis.analysis.adf.adf.ADF` class.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    #: str: The input file key of the second target (ligand-2) selection.
    target_selection_2_key = "target_selection_2"
    #: str: The input file key of the number of angle bins.
    n_angle_bins_key = "n_angle_bins"
    #: str: The input file key of the angle bin width in degrees.
    delta_angle_key = "delta_angle"
    #: str: The input file key of the lower i-j gate radius.
    r_min_1_key = "r_min_1"
    #: str: The input file key of the upper i-j gate radius.
    r_max_1_key = "r_max_1"
    #: str: The input file key of the lower i-k gate radius.
    r_min_2_key = "r_min_2"
    #: str: The input file key of the upper i-k gate radius.
    r_max_2_key = "r_max_2"

    #: List[str]: The required keys of the input file
    required_keys = [
        Reader.traj_files_key,
        Reader.reference_selection_key,
        Reader.target_selection_key,
        Reader.out_file_key,
    ]

    #: List[str]: The optional keys of the input file
    optional_keys = required_keys + [
        target_selection_2_key,
        n_angle_bins_key,
        delta_angle_key,
        r_min_1_key,
        r_max_1_key,
        r_min_2_key,
        r_max_2_key,
        Reader.log_file_key,
        Reader.use_full_atom_info_key,
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
        """
        super().read()
        super().check_required_keys(self.required_keys)
        super().check_known_keys(self.required_keys + self.optional_keys)
        super().not_defined_optional_keys(self.optional_keys)

    @property
    def target_selection_2(self) -> str | None:
        """str | None: The second target (ligand-2) selection."""
        return _parse_string(self.dictionary, self.target_selection_2_key)

    @property
    def n_angle_bins(self) -> PositiveInt | None:
        """PositiveInt | None: The number of angle bins."""
        return _parse_positive_int(self.dictionary, self.n_angle_bins_key)

    @property
    def delta_angle(self) -> PositiveReal | None:
        """PositiveReal | None: The angle bin width in degrees."""
        return _parse_positive_real(self.dictionary, self.delta_angle_key)

    @property
    def r_min_1(self) -> PositiveReal | None:
        """PositiveReal | None: The lower i-j gate radius."""
        return _parse_positive_real(self.dictionary, self.r_min_1_key)

    @property
    def r_max_1(self) -> PositiveReal | None:
        """PositiveReal | None: The upper i-j gate radius."""
        return _parse_positive_real(self.dictionary, self.r_max_1_key)

    @property
    def r_min_2(self) -> PositiveReal | None:
        """PositiveReal | None: The lower i-k gate radius."""
        return _parse_positive_real(self.dictionary, self.r_min_2_key)

    @property
    def r_max_2(self) -> PositiveReal | None:
        """PositiveReal | None: The upper i-k gate radius."""
        return _parse_positive_real(self.dictionary, self.r_max_2_key)



# pylint: disable=invalid-name
input_keys_documentation = f"""

For the ADF analysis input file several keys are available of which some are required and some are optional. For more details on the grammar and syntax of the input file see :ref:`inputFile`.

.. list-table:: Required keys
    :header-rows: 1

    * - Key
        - Value
    * - {Reader.traj_files_key}
        - The trajectory files to read. This can be a single file or a list of files.
    * - {Reader.reference_selection_key}
        - The selection string to select the reference (center ``i``) atoms. For more details see :py:class:`~PQAnalysis.topology.selection.Selection`.
    * - {Reader.target_selection_key}
        - The selection string to select the first ligand (``j``) atoms. For more details see :py:class:`~PQAnalysis.topology.selection.Selection`.
    * - {Reader.out_file_key}
        - The output file to write the ADF data to.

.. list-table:: Optional keys
    :header-rows: 1

    * - Key
        - Value
    * - {ADFInputFileReader.target_selection_2_key}
        - The selection string to select the second ligand (``k``) atoms. If omitted, the first target selection is reused so that the ``j-i-k`` angle is taken within one ligand set.
    * - {ADFInputFileReader.n_angle_bins_key}
        - The number of angle bins spanning 0 to 180 degrees (default 180).
    * - {ADFInputFileReader.delta_angle_key}
        - The width (in degrees) of the angle bins. Mutually exclusive
          with {ADFInputFileReader.n_angle_bins_key}.
    * - {ADFInputFileReader.r_min_1_key}
        - The lower (inclusive) ``i-j`` gate radius. Setting either ``i-j`` gate bound activates the ``i-j`` radial gate.
    * - {ADFInputFileReader.r_max_1_key}
        - The upper (exclusive) ``i-j`` gate radius.
    * - {ADFInputFileReader.r_min_2_key}
        - The lower (inclusive) ``i-k`` gate radius. Setting either ``i-k`` gate bound activates the ``i-k`` radial gate.
    * - {ADFInputFileReader.r_max_2_key}
        - The upper (exclusive) ``i-k`` gate radius.
    * - {Reader.log_file_key}
        - The log file to write the log information to.
    * - {Reader.use_full_atom_info_key}
        - Whether to use full atom information for the selections.

Note
----
Optional keys does not mean that they are optional for the analysis.
They are optional in the input file, but they might be required for
the analysis. For example:

- {ADFInputFileReader.n_angle_bins_key} and {ADFInputFileReader.delta_angle_key}
  are mutually exclusive; if neither is given, {ADFInputFileReader.n_angle_bins_key}
  defaults to 180.
- The radial gates are optional; if a gate bound is omitted the gate
  extends to 0 (lower) or infinity (upper) on that side.

"""

ADFInputFileReader.__doc__ += input_keys_documentation
