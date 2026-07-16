"""
A module containing a class to read input files to setup the
:py:class:`~PQAnalysis.analysis.msd.msd.MSD` class.
"""
import logging

# local imports
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import PQAnalysisInputFileReader as Reader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis.io.input_file_reader.pq_analysis._parse import (
    _parse_int,
    _parse_positive_int,
    _parse_positive_real,
)
from PQAnalysis.types import PositiveInt, PositiveReal
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking



class MSDInputFileReader(Reader):

    """
    A class to read input files to setup the
    :py:class:`~PQAnalysis.analysis.msd.msd.MSD` class.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    window_key = "window"
    gap_key = "gap"
    first_frame_key = "first_frame"
    start_key = "start"
    time_step_key = "time_step"
    fit_window_key = "fit_window"

    #: List[str]: The required keys of the input file
    required_keys = [
        Reader.traj_files_key,
        Reader.target_selection_key,
        Reader.out_file_key,
    ]

    #: List[str]: The optional keys of the input file
    optional_keys = required_keys + [
        window_key,
        gap_key,
        first_frame_key,
        start_key,
        time_step_key,
        fit_window_key,
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
        InputFileError
            if both the first_frame and the start key are set
        InputFileError
            if the first_frame or start key is negative
        InputFileError
            if the time_step key is not positive
        """
        super().read()
        super().check_required_keys(self.required_keys)
        super().check_known_keys(self.required_keys + self.optional_keys)
        super().not_defined_optional_keys(self.optional_keys)

        first_frame = _parse_int(self.dictionary, self.first_frame_key)
        start = _parse_int(self.dictionary, self.start_key)

        if first_frame is not None and start is not None:
            self.logger.error(
                (
                    f"The keys '{self.first_frame_key}' and "
                    f"'{self.start_key}' are aliases and cannot be "
                    "used at the same time."
                ),
                exception=InputFileError,
            )

        if self.n_start is not None and self.n_start < 0:
            self.logger.error(
                (
                    f"The '{self.first_frame_key}'/'{self.start_key}' "
                    "value has to be a non-negative integer - "
                    f"It actually is {self.n_start}!"
                ),
                exception=InputFileError,
            )

        if self.time_step is not None and self.time_step <= 0.0:
            self.logger.error(
                (
                    f"The '{self.time_step_key}' value has to be a "
                    "positive real number - "
                    f"It actually is {self.time_step}!"
                ),
                exception=InputFileError,
            )

    @property
    def window(self) -> PositiveInt | None:
        """
        PositiveInt | None: The correlation window size in frames.
        """
        return _parse_positive_int(self.dictionary, self.window_key)

    @property
    def gap(self) -> PositiveInt | None:
        """
        PositiveInt | None: The gap between two time origins in frames.
        """
        return _parse_positive_int(self.dictionary, self.gap_key)

    @property
    def n_start(self) -> int | None:
        """
        int | None: The first frame at which processing starts.
        """
        first_frame = _parse_int(self.dictionary, self.first_frame_key)

        if first_frame is not None:
            return first_frame

        return _parse_int(self.dictionary, self.start_key)

    @property
    def time_step(self) -> PositiveReal | None:
        """
        PositiveReal | None: The time step between two frames in ps.
        """
        return _parse_positive_real(self.dictionary, self.time_step_key)

    @property
    def fit_window(self) -> PositiveInt | None:
        """
        PositiveInt | None: The number of trailing MSD points to fit.
        """
        return _parse_positive_int(self.dictionary, self.fit_window_key)



input_keys_documentation = f"""

For the MSD analysis input file several keys are available of which some are required and some are optional. For more details on the grammar and syntax of the input file see :ref:`inputFile`.

.. list-table:: Required keys
    :header-rows: 1

    * - Key
      - Value
    * - {Reader.traj_files_key}
      - The trajectory files to read. This can be a single file or a list of files.
    * - {Reader.target_selection_key}
      - The selection string to select the atoms for which the MSD is calculated. For more details see :py:class:`~PQAnalysis.topology.selection.Selection`.
    * - {Reader.out_file_key}
      - The output file to write the MSD data to. It must not exist yet.

.. list-table:: Optional keys
    :header-rows: 1

    * - Key
      - Value
    * - {MSDInputFileReader.window_key}
      - The correlation window size in frames. Default is 1000.
        It has to be an integer multiple of the gap.
    * - {MSDInputFileReader.gap_key}
      - The gap between two time origins in frames. Default is 10.
    * - {MSDInputFileReader.first_frame_key}
      - The first frame (1-based frame counter) at which processing
        starts. Default is 0. The key {MSDInputFileReader.start_key}
        is an alias for this key (legacy Diffcalc naming).
    * - {MSDInputFileReader.time_step_key}
      - The time step between two frames in ps. If given, diffusion
        coefficients are calculated from a linear fit of the MSD tail
        and written to the log output. It has to be positive.
    * - {MSDInputFileReader.fit_window_key}
      - The number of trailing MSD points used for the diffusion fit.
        Default is the last 20% of the window. It has to be at least
        2 and at most window + 1.
    * - {Reader.log_file_key}
      - The log file to write the log information to.
    * - {Reader.use_full_atom_info_key}
      - Whether to use full atom information for the selections.

Note
----
The MSD output file follows the legacy Diffcalc format: one row per
lag index with the columns lag index, MSD in x, MSD in y and MSD in z
(all in Angstrom^2).

"""

MSDInputFileReader.__doc__ += input_keys_documentation
