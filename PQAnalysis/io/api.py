"""
This module provides API functions for input/output 
handling of molecular dynamics simulations.
"""

import logging

from PQAnalysis.types import PositiveReal
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.exceptions import PQNotImplementedError

from .input_file_reader import PQInputFileReader as Reader
from .input_file_reader.formats import InputFileFormat

logger = logging.getLogger(__package_name__).getChild("io.api")
logger = setup_logger(logger)



@runtime_type_checking
def continue_input_file(
    input_file: str,
    n: PositiveReal = 1,
    input_format: InputFileFormat | str = InputFileFormat.PQ
) -> None:
    """
    API function for continuing an input file.

    This function reads the input file and continues it n times.
    This means that the number in the filename is increased by one
    and all other numbers in the start- and output-file keys within 
    the input file are increased by one as well.

    Parameters
    ----------
    input_file : str
        the path to the input file, which should be continued
    n : PositiveReal, optional
        the number of times the input file should be continued, by default 1
    input_format : InputFileFormat | str, optional
        the format of the input file, by default InputFileFormat.PQ

    Raises
    ------
    NotImplementedError
        if the input format is not PQ
    """
    input_format = InputFileFormat(input_format)

    if input_format != InputFileFormat.PQ:
        logger.error(
            (
            f"Format {input_format} not implemented "
            "yet for continuing input file."
            ),
            exception=PQNotImplementedError
        )

    reader = Reader(input_file)
    reader.read()
    reader.continue_input_file(n)
