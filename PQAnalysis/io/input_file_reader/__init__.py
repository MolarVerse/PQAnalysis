"""
A package containing modules, classes and functions to 
parse and read input files of PQAnalysis itself and the md engines it supports.
"""

from .input_file_parser import InputFileParser, InputDictionary
from .pq import PQInputFileReader
from .pq_analysis import PQAnalysisInputFileReader
from .formats import InputFileFormat
