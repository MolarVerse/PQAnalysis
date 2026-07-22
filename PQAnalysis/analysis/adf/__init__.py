"""
A package containing classes and functions to handle angular distribution functions.

Classes
-------
:py:class:`~PQAnalysis.analysis.adf.adf.ADF`
    A class to handle angular distribution functions.
:py:class:`~PQAnalysis.analysis.adf.adf_input_file_reader.ADFInputFileReader`
    A class to read ADFs from input files.
:py:class:`~PQAnalysis.analysis.adf.adf_output_file_writer.ADFDataWriter`
    A class to write ADFs to output files.
:py:class:`~PQAnalysis.analysis.adf.adf_output_file_writer.ADFLogWriter`
    A class to write log files.

Functions
---------
:py:func:`~PQAnalysis.analysis.adf.api.adf`
    A function to create ADFs from an input file.
"""

from .api import adf
from .adf import ADF
from .adf_input_file_reader import ADFInputFileReader
from .adf_output_file_writer import ADFDataWriter, ADFLogWriter
from .exceptions import ADFError, ADFWarning

__all__ = [
    "adf",
    "ADF",
    "ADFInputFileReader",
    "ADFDataWriter",
    "ADFLogWriter",
    "ADFError",
    "ADFWarning",
]
