"""
Input-file reader for vibrational analysis.
"""

import logging

from PQAnalysis import __package_name__
from PQAnalysis.exceptions import PQKeyError
from PQAnalysis.io import PQAnalysisInputFileReader as Reader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis.io.input_file_reader.pq_analysis._parse import (
    _parse_positive_int,
    _parse_positive_real,
    _parse_string,
)
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.type_checking import runtime_type_checking



class VibrationalAnalysisInputFileReader(Reader):

    """
    A class to read input files for vibrational analysis.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    structure_file_key = "structure_file"
    hessian_file_key = "hessian_file"
    unit_key = "unit"
    hessian_sign_key = "hessian_sign"
    normal_modes_file_key = "normal_modes_file"
    modes_prefix_key = "modes_prefix"
    modes_file_key = "modes_file"
    modes_key = "modes"
    modes_frames_key = "modes_frames"
    modes_amplitude_key = "modes_amplitude"
    modes_temperature_key = "modes_temperature"
    modes_threshold_key = "modes_threshold"

    required_keys = [
        structure_file_key,
        hessian_file_key,
        Reader.out_file_key,
    ]

    optional_keys = required_keys + [
        Reader.moldescriptor_file_key,
        unit_key,
        hessian_sign_key,
        normal_modes_file_key,
        modes_prefix_key,
        modes_file_key,
        modes_key,
        modes_frames_key,
        modes_amplitude_key,
        modes_temperature_key,
        modes_threshold_key,
    ]

    @runtime_type_checking
    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            The input file.
        """
        self.filename = filename
        super().__init__(filename)

    def read(self) -> None:
        """
        Read and validate the input file.
        """
        super().read()
        super().check_required_keys(self.required_keys)
        super().check_known_keys(self.optional_keys)
        super().not_defined_optional_keys(self.optional_keys)

        if self.unit.lower() not in {"kcal", "hartree", "ev"}:
            self.logger.error(
                "The unit key must be one of: kcal, hartree, ev.",
                exception=InputFileError,
            )

        if self.hessian_sign.lower() not in {
            "auto", "positive", "negative", "1", "-1"
        }:
            self.logger.error(
                (
                    "The hessian_sign key must be one of: auto, positive, "
                    "negative, 1, -1."
                ),
                exception=InputFileError,
            )

        if isinstance(self.modes, str) and self.modes not in {
            "all", "nonzero", "positive"
        }:
            self.logger.error(
                (
                    "The modes key must be one of: all, nonzero, positive, "
                    "an integer, a list of integers or a range."
                ),
                exception=InputFileError,
            )

        if isinstance(self.modes,
                      list) and any(mode < 1 for mode in self.modes):
            self.logger.error(
                "Mode numbers must be positive and one-based.",
                exception=InputFileError,
            )

        if self.modes_temperature is not None and self.modes_temperature <= 0.0:
            self.logger.error(
                "The modes_temperature key must be positive.",
                exception=InputFileError,
            )

    @property
    def structure_file(self) -> str:
        """
        str: The structure file.
        """
        return _parse_string(self.dictionary, self.structure_file_key)

    @property
    def hessian_file(self) -> str:
        """
        str: The Hessian file.
        """
        return _parse_string(self.dictionary, self.hessian_file_key)

    @property
    def unit(self) -> str:
        """
        str: The Hessian unit.
        """
        unit = _parse_string(self.dictionary, self.unit_key)
        return "kcal" if unit is None else unit

    @property
    def hessian_sign(self) -> str:
        """
        str: The Hessian sign convention.
        """
        hessian_sign = _parse_string(self.dictionary, self.hessian_sign_key)
        return "auto" if hessian_sign is None else hessian_sign

    @property
    def normal_modes_file(self) -> str | None:
        """
        str | None: The normal-mode matrix output file.
        """
        return _parse_string(self.dictionary, self.normal_modes_file_key)

    @property
    def modes_prefix(self) -> str | None:
        """
        str | None: The XYZ mode file prefix.
        """
        return _parse_string(self.dictionary, self.modes_prefix_key)

    @property
    def modes_file(self) -> str | None:
        """
        str | None: The extended XYZ mode file.
        """
        return _parse_string(self.dictionary, self.modes_file_key)

    @property
    def modes(self) -> str | list[int]:
        """
        str | list[int]: The normal modes selected for visualization.
        """
        modes = _parse_modes(self.dictionary, self.modes_key)
        return "all" if modes is None else modes

    @property
    def modes_frames(self) -> int:
        """
        int: The number of frames per animated mode trajectory.
        """
        modes_frames = _parse_positive_int(
            self.dictionary,
            self.modes_frames_key,
        )
        return 30 if modes_frames is None else modes_frames

    @property
    def modes_amplitude(self) -> float:
        """
        float: The maximum displacement in animated mode trajectories.
        """
        modes_amplitude = _parse_positive_real(
            self.dictionary,
            self.modes_amplitude_key,
        )
        return 0.25 if modes_amplitude is None else float(modes_amplitude)

    @property
    def modes_temperature(self) -> float | None:
        """
        float | None: The ASE-style temperature used to scale modes.
        """
        modes_temperature = _parse_positive_real(
            self.dictionary,
            self.modes_temperature_key,
        )
        return None if modes_temperature is None else float(modes_temperature)

    @property
    def modes_threshold(self) -> float:
        """
        float: Wavenumber threshold for named mode selections.
        """
        modes_threshold = _parse_positive_real(
            self.dictionary,
            self.modes_threshold_key,
        )
        return 1.0e-8 if modes_threshold is None else float(modes_threshold)



def _parse_modes(input_dict, key: str) -> str | list[int] | None:
    """
    Parse mode selection input.
    """
    try:
        data = input_dict[key]
    except PQKeyError:
        return None

    value, data_type, _ = data

    if data_type == "None":
        return None

    if data_type == "str":
        return str(value).lower()

    if data_type == "int":
        return [int(value)]

    if data_type == "list(int)":
        return [int(mode) for mode in value]

    if data_type == "range":
        return [int(mode) for mode in value]

    VibrationalAnalysisInputFileReader.logger.error(
        (
            "The modes key must be one of: all, nonzero, positive, "
            "an integer, a list of integers or a range."
        ),
        exception=InputFileError,
    )
    return None



input_keys_documentation = f"""

For the vibrational analysis input file several keys are available.

.. list-table:: Required keys
    :header-rows: 1

    * - Key
      - Value
    * - {VibrationalAnalysisInputFileReader.structure_file_key}
      - The restart or single-frame XYZ structure file.
    * - {VibrationalAnalysisInputFileReader.hessian_file_key}
      - The plain square Hessian matrix file.
    * - {Reader.out_file_key}
      - The tabular output file.

.. list-table:: Optional keys
    :header-rows: 1

    * - Key
      - Value
    * - {Reader.moldescriptor_file_key}
      - The moldescriptor file used for IR intensities.
    * - {VibrationalAnalysisInputFileReader.unit_key}
      - Hessian unit. Options are kcal, hartree and ev. Default is kcal.
    * - {VibrationalAnalysisInputFileReader.hessian_sign_key}
      - Hessian sign convention. Options are auto, positive, negative, 1 and -1.
    * - {VibrationalAnalysisInputFileReader.normal_modes_file_key}
      - Optional matrix-format normal-mode output file.
    * - {VibrationalAnalysisInputFileReader.modes_prefix_key}
      - Optional prefix for one sinusoidal XYZ animation file per selected mode.
    * - {VibrationalAnalysisInputFileReader.modes_file_key}
      - Optional extended XYZ file with all selected mode vectors and metadata.
    * - {VibrationalAnalysisInputFileReader.modes_key}
      - Modes to write. Options are all, nonzero, positive, one integer,
        a list of integers or a range. Explicit mode numbers are one-based.
    * - {VibrationalAnalysisInputFileReader.modes_frames_key}
      - Number of frames per animated mode trajectory. Default is 30.
    * - {VibrationalAnalysisInputFileReader.modes_amplitude_key}
      - Maximum displacement in Angstrom for fixed-amplitude animations.
        Default is 0.25.
    * - {VibrationalAnalysisInputFileReader.modes_temperature_key}
      - Optional temperature in Kelvin for ASE-style energy-scaled animations.
    * - {VibrationalAnalysisInputFileReader.modes_threshold_key}
      - Wavenumber threshold in cm-1 for nonzero and positive selections.
        Default is 1.0e-8.

"""

VibrationalAnalysisInputFileReader.__doc__ += input_keys_documentation
