"""
This module provides API functions for the finite difference analysis.
"""
# 3rd party imports
from beartype.typing import Tuple
import numpy as np

from PQAnalysis.io import TrajectoryReader, RestartFileReader, MoldescriptorReader
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.topology import Topology
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.types import Np1DNumberArray

from .finite_difference import FiniteDifference
from .finite_difference_input_file_reader import FiniteDifferenceInputFileReader
from .finite_difference_output_file_writer import FiniteDifferenceDataWriter, FiniteDifferenceLogWriter


@runtime_type_checking
def finite_difference(input_file: str, md_format: MDEngineFormat | str = MDEngineFormat.PQ):
    """
    Calculates the finite difference using a given input file.

    This is just a wrapper function combining the underlying classes and functions.

    For more information on the input file keys please
    visit :py:mod:`~PQAnalysis.analysis.finite_difference.finite_difference_input_file_reader`.
    For more information on the exact calculation of
    the finite difference please visit :py:class:`~PQAnalysis.analysis.finite_difference.finite_difference.FiniteDifference`.

    Parameters
    ----------
    input_file : str
        The input file. For more information on the input file
        keys please visit :py:mod:`~PQAnalysis.analysis.finite_difference.finite_difference_input_file_reader`.
    md_format : MDEngineFormat | str, optional
        the format of the input trajectory. Default is "PQ".
        For more information on the supported formats please visit
        :py:class:`~PQAnalysis.traj.formats.MDEngineFormat`.    
    """
    md_format = MDEngineFormat(md_format)

    input_reader = FiniteDifferenceInputFileReader(input_file)
    input_reader.read()
    finite_difference_points = input_reader.finite_difference_points

    temperature_points = input_reader.temperature_points

    if (np.all(np.diff(temperature_points) == temperature_points[1] - temperature_points[0])):
        raise ValueError(
            "The temperature points must be equi-partition steps.")
    if input_reader.std_points_key is not None:
        std_points = input_reader.std_points
    else:
        std_points = None
    _finite_difference = FiniteDifference(temperature_points=temperature_points,
                                          finite_difference_points=finite_difference_points,
                                          std_points=std_points)

    data_writer = FiniteDifferenceDataWriter(input_reader.out_file)
    log_writer = FiniteDifferenceLogWriter(input_reader.log_file)
    log_writer.write_before_run(_finite_difference)

    finite_difference_data = _finite_difference.run()

    data_writer.write(finite_difference_data)
    log_writer.write_after_run(finite_difference_data)
