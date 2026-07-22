"""
A package containing classes and functions to handle mean square
displacement (MSD) analyses.

Classes
-------
:py:class:`~PQAnalysis.analysis.msd.msd.MSD`
    A class to calculate mean square displacements.
:py:class:`~PQAnalysis.analysis.msd.msd.MSDDiffusionFit`
    A container for the result of a linear diffusion fit.
:py:class:`~PQAnalysis.analysis.msd.msd_input_file_reader.MSDInputFileReader`
    A class to read MSD setups from input files.
:py:class:`~PQAnalysis.analysis.msd.msd_output_file_writer.MSDDataWriter`
    A class to write MSD data to output files.
:py:class:`~PQAnalysis.analysis.msd.msd_output_file_writer.MSDLogWriter`
    A class to write log files.

Functions
---------
:py:func:`~PQAnalysis.analysis.msd.api.msd`
    A function to calculate MSDs from an input file.
"""

from .api import msd
from .msd import MSD, MSDDiffusionFit
from .msd_input_file_reader import MSDInputFileReader
from .msd_output_file_writer import MSDDataWriter, MSDLogWriter
from .exceptions import MSDError, MSDWarning

__all__ = [
    "MSD",
    "MSDDataWriter",
    "MSDDiffusionFit",
    "MSDError",
    "MSDInputFileReader",
    "MSDLogWriter",
    "MSDWarning",
    "msd",
]
