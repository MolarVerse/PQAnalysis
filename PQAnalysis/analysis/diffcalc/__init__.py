"""
A package containing classes and functions to handle radial distribution functions.

Classes
-------
:py:class:`~PQAnalysis.analysis.diffcalc.diffcalc.diffcalc`
    A class to handle radial distribution functions.
:py:class:`~PQAnalysis.analysis.diffcalc.diffcalcInputFileReader.diffcalcInputFileReader`
    A class to read diffcalcs from input files.
:py:class:`~PQAnalysis.analysis.diffcalc.diffcalcOutputFileWriter.diffcalcDataWriter`
    A class to write diffcalcs to output files.
:py:class:`~PQAnalysis.analysis.diffcalc.diffcalcOutputFileWriter.diffcalcLogWriter`
    A class to write log files.
   
Functions
---------
:py:func:`~PQAnalysis.analysis.diffcalc.diffcalc.diffcalc`
    A function to create diffcalcs from an input file. 
"""

from .api import diffcalc
from .diffcalc import diffcalc
from .diffcalcInputFileReader import DiffcalcInputFileReader
from .diffcalcOutputFileWriter import DiffcalcDataWriter, DiffcalcLogWriter
from .exceptions import *
