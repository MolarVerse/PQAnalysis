"""
This module provides API functions for input/output handling of molecular dynamics simulations.
"""

from .inputFileReader import PQ_InputFileReader as Reader
from .inputFileReader.formats import InputFileFormat
from .write_api import write, write_box
from .conversion_api import rst2xyz, traj2box, traj2qmcfc

from PQAnalysis.types import PositiveReal


def continue_input_file(input_file: str,
                        n: PositiveReal = 1,
                        input_format: InputFileFormat | str = InputFileFormat.PQ
                        ) -> None:
    """
    API function for continuing an input file.

    This function reads the input file and continues it n times. This means that the number in the filename is increased by one and all other numbers in the start- and output-file keys within the input file are increased by one as well.

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
        raise NotImplementedError(
            f"Format {input_format} not implemented yet for continuing input file.")

    reader = Reader(input_file)
    reader.read()
    reader.continue_input_file(n)
