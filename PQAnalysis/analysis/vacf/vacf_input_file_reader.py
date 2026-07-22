"""
A module containing a class to read input files to setup the
:py:class:`~PQAnalysis.analysis.vacf.vacf.VACF` class.
"""
import logging

# 3rd party imports
from beartype.typing import List

# local imports
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import PQAnalysisInputFileReader as Reader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis.io.input_file_reader.pq_analysis._parse import (
    _parse_files,
    _parse_positive_int,
    _parse_positive_real,
    _parse_string,
)
from PQAnalysis.types import PositiveInt, PositiveReal
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

# local relative imports
from .spectrum import WINDOW_FUNCTIONS



class VACFInputFileReader(Reader):

    """
    A class to read input files to setup the
    :py:class:`~PQAnalysis.analysis.vacf.vacf.VACF` class.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    window_key = "window"
    gap_key = "gap"
    time_step_key = "time_step"
    method_key = "method"
    charge_file_key = "charge_file"
    charge_files_key = "charge_files"
    spectrum_file_key = "spectrum_file"
    ftsize_key = "ftsize"
    window_function_key = "window_function"
    window_param_key = "window_param"
    window_start_key = "window_start"
    window_stop_key = "window_stop"
    windowed_out_file_key = "windowed_out_file"

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
        charge_file_key,
        charge_files_key,
        spectrum_file_key,
        ftsize_key,
        window_function_key,
        window_param_key,
        window_start_key,
        window_stop_key,
        windowed_out_file_key,
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
            if both a static charge file and charge trajectory
            files are given
        InputFileError
            if a windowed output file is requested without a
            spectrum file
        InputFileError
            if the window function is unknown
        """
        super().read()
        super().check_required_keys(self.required_keys)
        super().check_known_keys(self.required_keys + self.optional_keys)
        super().not_defined_optional_keys(self.optional_keys)

        if self.charge_file is not None and self.charge_files is not None:
            self.logger.error(
                (
                    f"The keys '{self.charge_file_key}' and "
                    f"'{self.charge_files_key}' cannot be used at the "
                    "same time. Please provide only one charge source "
                    "for the charge-flux mode."
                ),
                exception=InputFileError,
            )

        if (
            self.windowed_out_file is not None and
            self.spectrum_file is None
        ):
            self.logger.error(
                (
                    f"The key '{self.windowed_out_file_key}' can only "
                    "be used together with the key "
                    f"'{self.spectrum_file_key}'."
                ),
                exception=InputFileError,
            )

        if (
            self.window_function is not None and
            self.window_function.lower() not in WINDOW_FUNCTIONS
        ):
            self.logger.error(
                (
                    f"Unknown window function '{self.window_function}'. "
                    "Possible window functions are: "
                    f"{', '.join(WINDOW_FUNCTIONS)}."
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
    def time_step(self) -> PositiveReal | None:
        """
        PositiveReal | None: The time step between two frames in ps.
        """
        return _parse_positive_real(self.dictionary, self.time_step_key)

    @property
    def method(self) -> str | None:
        """
        str | None: The VACF estimator method (direct or fft).
        """
        return _parse_string(self.dictionary, self.method_key)

    @property
    def charge_file(self) -> str | None:
        """
        str | None: The static charge file for the charge-flux mode.
        """
        return _parse_string(self.dictionary, self.charge_file_key)

    @property
    def charge_files(self) -> List[str] | None:
        """
        List[str] | None: The charge trajectory files for the charge-flux mode.
        """
        return _parse_files(self.dictionary, self.charge_files_key)

    @property
    def spectrum_file(self) -> str | None:
        """
        str | None: The output file for the VACF spectrum.
        """
        return _parse_string(self.dictionary, self.spectrum_file_key)

    @property
    def ftsize(self) -> PositiveInt | None:
        """
        PositiveInt | None: The Fourier transform point size.
        """
        return _parse_positive_int(self.dictionary, self.ftsize_key)

    @property
    def window_function(self) -> str | None:
        """
        str | None: The apodization window function of the spectrum.
        """
        return _parse_string(self.dictionary, self.window_function_key)

    @property
    def window_param(self) -> PositiveReal | None:
        """
        PositiveReal | None: The exponential window decay coefficient.
        """
        return _parse_positive_real(self.dictionary, self.window_param_key)

    @property
    def window_start(self) -> PositiveReal | None:
        """
        PositiveReal | None: The window start time in ps.
        """
        return _parse_positive_real(self.dictionary, self.window_start_key)

    @property
    def window_stop(self) -> PositiveReal | None:
        """
        PositiveReal | None: The window stop time in ps.
        """
        return _parse_positive_real(self.dictionary, self.window_stop_key)

    @property
    def windowed_out_file(self) -> str | None:
        """
        str | None: The output file for the windowed correlation function.
        """
        return _parse_string(self.dictionary, self.windowed_out_file_key)



input_keys_documentation = f"""

For the VACF analysis input file several keys are available of which some are required and some are optional. For more details on the grammar and syntax of the input file see :ref:`inputFile`.

.. list-table:: Required keys
    :header-rows: 1

    * - Key
      - Value
    * - {Reader.traj_files_key}
      - The velocity trajectory files to read. This can be a single
        file or a list of files.
    * - {Reader.target_selection_key}
      - The selection string to select the atoms for which the VACF is
        calculated. For more details see
        :py:class:`~PQAnalysis.topology.selection.Selection`.
    * - {Reader.out_file_key}
      - The output file to write the VACF data to. It must not exist yet.
    * - {VACFInputFileReader.time_step_key}
      - The time step between two frames in ps. It is used to build the
        time axis of the VACF and the frequency axis of the spectrum.

.. list-table:: Optional keys
    :header-rows: 1

    * - Key
      - Value
    * - {VACFInputFileReader.window_key}
      - The correlation window size in frames. Default is 1000.
        It has to be an integer multiple of the gap.
    * - {VACFInputFileReader.gap_key}
      - The gap between two time origins in frames. Default is 1.
    * - {VACFInputFileReader.method_key}
      - The VACF estimator, either "direct" (legacy-exact sliding
        origins, default) or "fft" (denser-origin Wiener-Khinchin
        estimator, ignores the gap).
    * - {VACFInputFileReader.charge_file_key}
      - A static charge file for the charge-flux mode. The file
        consists of a header line with the number of atoms, a comment
        line and one "name charge" line per atom.
    * - {VACFInputFileReader.charge_files_key}
      - Charge trajectory files (.chrg) for the charge-flux mode with
        time-dependent charges, read in lockstep with the velocity
        trajectory. Cannot be combined with
        {VACFInputFileReader.charge_file_key}.
    * - {VACFInputFileReader.spectrum_file_key}
      - The output file for the cosine-transform spectrum of the VACF.
        If not given, no spectrum is calculated.
    * - {VACFInputFileReader.ftsize_key}
      - The Fourier transform point size of the spectrum. The VACF is
        zero-padded (or truncated) to this size. Default is 2000.
    * - {VACFInputFileReader.window_function_key}
      - The apodization window function of the spectrum, one of
        "none", "exponential", "hann" and "blackman". Default is
        "none".
    * - {VACFInputFileReader.window_param_key}
      - The decay coefficient of the exponential window. Default is 4.0.
    * - {VACFInputFileReader.window_start_key}
      - The time in ps at which the apodization window starts to
        decay. Default is 0.0.
    * - {VACFInputFileReader.window_stop_key}
      - The time in ps at which the apodization window becomes zero.
        Default is 1000.0.
    * - {VACFInputFileReader.windowed_out_file_key}
      - The output file for the windowed VACF. Can only be used
        together with {VACFInputFileReader.spectrum_file_key}.
    * - {Reader.log_file_key}
      - The log file to write the log information to.
    * - {Reader.use_full_atom_info_key}
      - Whether to use full atom information for the selections.

Note
----
The VACF output file follows the legacy FreqCalc format: one row per
lag with the columns lag time in ps and normalized VACF. The spectrum
output file follows the legacy ft.f format: one row per frequency
index with the columns wavenumber in cm^-1 and amplitude.

"""

VACFInputFileReader.__doc__ += input_keys_documentation
