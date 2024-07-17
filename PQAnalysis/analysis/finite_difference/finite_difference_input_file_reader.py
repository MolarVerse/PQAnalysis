"""
A module to read input files to setup the :py:class:`~PQAnalysis.analysis.finite_difference.finite_difference.FiniteDifference` class.
"""
import logging

# local imports
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import PQAnalysisInputFileReader as Reader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking


class FiniteDifferenceInputFileReader(Reader):

    """
    A class to read input files to setup the :py:class:`~PQAnalysis.analysis.finite_difference.finite_difference.FiniteDifference` class.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    #: List[str]: The required keys of the input file
    required_keys = [
        Reader.temperature_points_key,
        Reader.finite_difference_points_key,
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
            If a required key is missing or an unknown key is present.
        InputFileWarning
            if unknown keys are set in the input file
        """

        super().read()
        super().check_required_keys(self.required_keys)
        super().check_known_keys(self.required_keys + self.optional_keys)


input_keys_documentation = f"""

For the finite difference analysis the following keys are required:
.. list-table::
    :header-rows: 1

    * - Key
        - Value
    * - {Reader.temperature_points_key}
        - The temperature points to use for the finite difference analysis.
    * - {Reader.finite_difference_points_key}
        - The number of points to use for the finite difference analysis.

    * - {Reader.out_file_key}
        - The output file to write finite differences to.

.. list-table:: Optional keys
    :header-rows: 1

    * - Key
        - Value
    * - {Reader.std_points_key}
        - The standard deviation points to use for the finite difference analysis.
    * - {Reader.log_file_key}
        - The log file to write the log information to.

"""

FiniteDifferenceInputFileReader.__doc__ += input_keys_documentation
