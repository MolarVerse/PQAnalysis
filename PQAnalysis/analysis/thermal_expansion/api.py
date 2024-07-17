"""
This module provides API functions for the linear thermal expansion coefficient analysis.
"""
import numpy as np

from PQAnalysis.io import BoxFileReader
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.type_checking import runtime_type_checking


from .thermal_expansion import ThermalExpansion
from .thermal_expansion_input_file_reader import ThermalExpansionInputFileReader
from .thermal_expansion_output_file_writer import ThermalExpansionDataWriter
from .thermal_expansion_output_file_writer import ThermalExpansionLogWriter


@runtime_type_checking
def thermal_expansion(input_file: str, md_format: MDEngineFormat | str = MDEngineFormat.PQ):
    """
    Calculates the linear thermal expansion coefficient using a given input file.

    This is just a wrapper function combining the underlying classes and functions.

    For more information on the input file keys please
    visit :py:mod:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion_input_file_reader`.
    For more information on the exact calculation of
    the linear thermal expansion coefficient please visit :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion`.

    Parameters
    ----------
    input_file : str
        The input file. For more information on the input file
        keys please visit :py:mod:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion_input_file_reader`.
    md_format : MDEngineFormat | str, optional
        the format of the input trajectory. Default is "PQ".
        For more information on the supported formats please visit
        :py:class:`~PQAnalysis.traj.formats.MDEngineFormat`.    
    """
    md_format = MDEngineFormat(md_format)

    input_reader = ThermalExpansionInputFileReader(input_file)
    input_reader.read()

    if len(input_reader.box_files) != len(input_reader.temperature_points):
        raise ValueError(
            "The number of box files and temperature points must be the same.")

    # check if input_reader.temperature_points is equi partitioned steps between each temperature point should be the same
    temperature_points = input_reader.temperature_points

    if (np.all(np.diff(temperature_points) == temperature_points[1] - temperature_points[0])):
        raise ValueError(
            "The temperature points must be equi-partition steps.")

    box_data = []
    for box_file in input_reader.box_files:
        box_reader = BoxFileReader(box_file)
        box = box_reader.read()
        box_data.append(box)

    _thermal_expansion = ThermalExpansion(
        box_data, temperature_points)

    data_writer = ThermalExpansionDataWriter(input_reader.out_file_key)
    log_writer = ThermalExpansionLogWriter(input_reader.log_file)
    log_writer.write_before_run(_thermal_expansion)

    box_av_data, box_std_data, thermal_expansion_data = _thermal_expansion.run()

    data_writer.write(box_av_data, box_std_data, thermal_expansion_data)
    log_writer.write_after_run(_thermal_expansion)
