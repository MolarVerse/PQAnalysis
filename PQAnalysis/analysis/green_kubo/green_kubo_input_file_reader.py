"""
A module containing a class to read input files to setup the
:py:class:`~PQAnalysis.analysis.green_kubo.green_kubo.GreenKubo` class.
"""
import logging

# local imports
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import PQAnalysisInputFileReader as Reader
from PQAnalysis.io.input_file_reader.pq_analysis._parse import (
    _parse_positive_int,
    _parse_positive_real,
    _parse_string,
)
from PQAnalysis.types import PositiveInt, PositiveReal
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking



class GreenKuboInputFileReader(Reader):

    """
    A class to read input files to setup the
    :py:class:`~PQAnalysis.analysis.green_kubo.green_kubo.GreenKubo`
    class.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    window_key = "window"
    gap_key = "gap"
    time_step_key = "time_step"
    method_key = "method"
    coefficient_key = "coefficient"
    fit_start_key = "fit_start"
    fit_stop_key = "fit_stop"
    n_blocks_key = "n_blocks"

    #: List[str]: The required keys of the input file
    required_keys = [
        Reader.traj_files_key,
        Reader.target_selection_key,
        Reader.out_file_key,
        time_step_key,
    ]

    #: List[str]: The optional keys of the input file
    optional_keys = required_keys + [
        window_key,
        gap_key,
        method_key,
        coefficient_key,
        fit_start_key,
        fit_stop_key,
        n_blocks_key,
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
    def time_step(self) -> PositiveReal | None:
        """
        PositiveReal | None: The time step between two frames in ps.
        """
        return _parse_positive_real(self.dictionary, self.time_step_key)

    @property
    def method(self) -> str | None:
        """
        str | None: The auto-correlation estimator (fft or direct).
        """
        return _parse_string(self.dictionary, self.method_key)

    @property
    def coefficient(self) -> str | None:
        """
        str | None: The transport coefficient to calculate (diffusion).
        """
        return _parse_string(self.dictionary, self.coefficient_key)

    @property
    def fit_start(self) -> PositiveReal | None:
        """
        PositiveReal | None: The plateau fit window start fraction.
        """
        return _parse_positive_real(self.dictionary, self.fit_start_key)

    @property
    def fit_stop(self) -> PositiveReal | None:
        """
        PositiveReal | None: The plateau fit window stop fraction.
        """
        return _parse_positive_real(self.dictionary, self.fit_stop_key)

    @property
    def n_blocks(self) -> PositiveInt | None:
        """
        PositiveInt | None: The number of blocks for the block-averaged
        standard error of the diffusion coefficient.
        """
        return _parse_positive_int(self.dictionary, self.n_blocks_key)



# the lower_case module doc-string name is shared by all analysis
# input file readers (msd, vacf, ...) and imported by the cli tools
# pylint: disable=invalid-name
input_keys_documentation = f"""

For the Green-Kubo analysis input file several keys are available of which some are required and some are optional. For more details on the grammar and syntax of the input file see :ref:`inputFile`.

.. list-table:: Required keys
    :header-rows: 1

    * - Key
      - Value
    * - {Reader.traj_files_key}
      - The velocity trajectory files (.vel) to read. This can be a
        single file or a list of files. The velocities have to be given
        in Angstrom / s (the PQ velocity output unit).
    * - {Reader.target_selection_key}
      - The selection string to select the atoms for which the
        Green-Kubo diffusion coefficient is calculated. For more details
        see :py:class:`~PQAnalysis.topology.selection.Selection`.
    * - {Reader.out_file_key}
      - The output file to write the running-integral data to. It must
        not exist yet.
    * - {GreenKuboInputFileReader.time_step_key}
      - The time step between two frames in ps. It is used to build the
        lag-time axis and to integrate the velocity auto-correlation
        function.

.. list-table:: Optional keys
    :header-rows: 1

    * - Key
      - Value
    * - {GreenKuboInputFileReader.window_key}
      - The correlation window size in frames. Default is 1000 (capped
        to the number of frames minus one).
    * - {GreenKuboInputFileReader.gap_key}
      - The gap between two time origins in frames for the direct
        estimator. Default is 1. It is ignored by the fft estimator.
    * - {GreenKuboInputFileReader.method_key}
      - The auto-correlation estimator, either "fft" (default,
        Wiener-Khinchin) or "direct" (sliding time origins).
    * - {GreenKuboInputFileReader.coefficient_key}
      - The transport coefficient to calculate. Currently only
        "diffusion" is supported (default).
    * - {GreenKuboInputFileReader.fit_start_key}
      - The start of the plateau fit window as a fraction of the
        correlation window (0..1). Default is 0.5.
    * - {GreenKuboInputFileReader.fit_stop_key}
      - The end of the plateau fit window as a fraction of the
        correlation window (0..1). Default is 1.0.
    * - {GreenKuboInputFileReader.n_blocks_key}
      - The number of contiguous, non-overlapping blocks the trajectory
        is split into for the block-averaged standard error of the
        diffusion coefficient. Must be at least 2. Default is 5. It is
        clamped down (with a warning) if the blocks would be shorter
        than the correlation window.
    * - {Reader.log_file_key}
      - The log file to write the log information to.
    * - {Reader.use_full_atom_info_key}
      - Whether to use full atom information for the selections.

Note
----
The Green-Kubo output file contains one row per lag with the columns
lag time in ps, the un-normalized velocity auto-correlation function
Cvv in (Angstrom / s)^2 and the running diffusion coefficient
D_running in m^2 / s. The plateau diffusion coefficient with its
block-averaged standard error (over n_blocks blocks) is reported in
m^2 / s and cm^2 / s in the log output.

"""

GreenKuboInputFileReader.__doc__ += input_keys_documentation
