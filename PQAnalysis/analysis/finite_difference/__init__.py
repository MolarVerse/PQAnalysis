"""
A package containing classes and functions to handle finite difference.

Classes
-------
:py:class:`~PQAnalysis.analysis.finite_difference.finite_difference.FiniteDifference`
    A class to handle finite difference.
:py:class:`~PQAnalysis.analysis.finite_difference.finite_difference_input_file_reader.FiniteDifferenceInputFileReader`
    A class to read finite differences from input files.
:py:class:`~PQAnalysis.analysis.finite_difference.finite_difference_output_file_writer.FiniteDifferenceDataWriter`
    A class to write finite differences to output files.
:py:class:`~PQAnalysis.analysis.finite_difference.finite_difference_output_file_writer.FiniteDifferenceLogWriter`
    A class to write log files.

Functions
---------
:py:func:`~PQAnalysis.analysis.finite_difference.finite_difference.finite_difference`
    A function to create finite differences from an input file. 

"""

from .api import finite_difference
from .finite_difference import FiniteDifference
from .finite_difference_input_file_reader import FiniteDifferenceInputFileReader
from .finite_difference_output_file_writer import FiniteDifferenceDataWriter, FiniteDifferenceLogWriter
from .exceptions import *
