"""
A module to read input files to setup the :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion` class.
"""
import logging

# local imports
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import PQAnalysisInputFileReader as Reader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking


class ThermalExpansionInputFileReader(Reader):

    """
    A class to read input files to setup the :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion` class.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    #: List[str]: The required keys of the input file
    required_keys = [
        Reader.temperature_points_key,
        Reader.box_files_key,
        Reader.out_file_key,
    ]

    #: List[str]: The optional keys of the input file
    optional_keys = required_keys + [
        Reader.log_file_key,
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


input_keys_documentation = f"""

For the linear or volumetric thermal expansion coefficient analysis the following keys are required:
.. list-table::
    :header-rows: 1

    * - Key
        - Value
    * - {Reader.temperature_points_key}
        - The temperatures which corrosponds to the box_files.
    * - {Reader.box_files_key}
        - The file containing the box dimensions: :math:`a`, :math:`b`, :math:`c`, :math:`\\alpha`, :math:`\\beta`, :math:`\\gamma`.
    * - {Reader.out_file_key}
        - The output file to write linear thermal expansion coefficients to.

.. list-table:: Optional keys
    :header-rows: 1

    * - Key
        - Value
    * - {Reader.log_file_key}
        - The log file to write the log information to.
Note
----
Optional keys does not mean that they are optional for the analysis.
They are optional in the input file, but they might be required for
the analysis. This means that if an optional keyword is specified
other keywords might be required.
- :code:`{Reader.log_file_key}` is optional for the analysis.
(for more information see :py:class:`~PQAnalysis.io.input_file_reader.pq_analysis.thermal_expansion.thermal_expansion.ThermalExpansion`).

"""

ThermalExpansionInputFileReader.__doc__ += input_keys_documentation
