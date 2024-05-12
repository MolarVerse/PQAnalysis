"""
A package containing classes and functions to handle radial distribution functions.

Classes
-------
:py:class:`~PQAnalysis.analysis.rdf.rdf.RDF`
    A class to handle radial distribution functions.
:py:class:`~PQAnalysis.analysis.rdf.rdfInputFileReader.RDFInputFileReader`
    A class to read RDFs from input files.
:py:class:`~PQAnalysis.analysis.rdf.rdfOutputFileWriter.RDFDataWriter`
    A class to write RDFs to output files.
:py:class:`~PQAnalysis.analysis.rdf.rdfOutputFileWriter.RDFLogWriter`
    A class to write log files.
   
Functions
---------
:py:func:`~PQAnalysis.analysis.rdf.rdf.rdf`
    A function to create RDFs from an input file. 
"""

from .api import rdf
from .rdf import RDF
from .rdf_input_file_reader import RDFInputFileReader
from .rdf_output_file_writer import RDFDataWriter, RDFLogWriter
from .exceptions import *
