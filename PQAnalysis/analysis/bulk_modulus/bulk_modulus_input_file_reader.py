"""
A module containing a class to read input files to setup the 
:py:class:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus.BulkModulus` class.
"""
import logging

# local imports
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import PQAnalysisInputFileReader as Reader
from PQAnalysis.io.input_file_reader.exceptions import InputFileError
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking


class BulkModulusInputFileReader(Reader):

    """
    A class to read input files to setup the :py:class:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus.BulkModulus` class.
    """

    #: List[str]: The required keys of the input file
    required_keys = [
        Reader.energy_perturbation_files_key,
        Reader.volumes_perturbation_key,
        Reader.volume_equilibrium_key,
        Reader.out_file_key,
    ]

    #: List[str]: The optional keys of the input file
    optional_keys = required_keys + [
        Reader.mode_key,
        Reader.info_file_key,
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

For the BulkModulus analysis input file several keys are available of which some are required and some are optional.
For more details on the grammar and syntax of the input file see :ref:`inputFile`.

.. list-table:: Required keys
    :header-rows: 1

    * - Key
        - Value
    * - {Reader.energy_perturbation_files_key}
        - The energy files with small perturbations from the equilibrium state.
        - Each file is a different perturbation.
    * - {Reader.volumes_perturbation_key}
        - The volumes which show a small pertubation.
    * - {Reader.volume_equilibrium_key}
        - The equilibrium volume.
    * - {Reader.out_file_key}
        - The output file to write the BulkModulus data to.

.. list-table:: Optional keys
    :header-rows: 1

    * - Key
        - Value
    * - {Reader.mode_key}
        - The mode of the BulkModulus analysis.
        - Default is "simple": finite difference method of two points.
        - Possible values are:
            - "BMEOS": Birch-Murnaghan equation of state
            - "MEOS": Murnaghan equation of state
            - "simple": finite difference method of two points
    * - {Reader.info_file_key}
        - The info file of the energy files.
    * - {Reader.log_file_key}
        - The log file to write the log information to.


Note
----
Optional keys does not mean that they are optional for the analysis.
They are optional in the input file, but they might be required for
the analysis. This means that if an optional keyword is specified
other keywords might be required. For example:


(for more information see :py:class:`~PQAnalysis.analysis.bulk_modulus.bulk_modulus.BulkModulus`).

"""

BulkModulusInputFileReader.__doc__ += input_keys_documentation
