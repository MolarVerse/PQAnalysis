"""
This module provides API functions for
linear or volumetric thermal expansion coefficient analysis.
"""

import numpy as np

from PQAnalysis.io.box_file_reader import BoxFileReader
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.io.formats import FileWritingMode

from .thermal_expansion import ThermalExpansion
from .thermal_expansion_input_file_reader import ThermalExpansionInputFileReader
from .thermal_expansion_output_file_writer import ThermalExpansionDataWriter
from .thermal_expansion_output_file_writer import ThermalExpansionLogWriter



@runtime_type_checking
def thermal_expansion(
    input_file: str,
    md_format: MDEngineFormat | str = MDEngineFormat.PQ,
    mode: FileWritingMode | str = "w"
):
    """
    Calculates the thermal expansion coefficient using a given input file.

    This is just a wrapper function combining the underlying classes and functions.

    For more information on the input file keys please
    visit 
    :py:mod:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion_input_file_reader`.
    For more information on the exact calculation of the thermal expansion coefficient please visit 
    :py:class:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion.ThermalExpansion`.

    Parameters
    ----------
    input_file : str
        The input file. For more information on the input file
        keys please visit 
        :py:mod:`~PQAnalysis.analysis.thermal_expansion.thermal_expansion_input_file_reader`.
    md_format : MDEngineFormat | str, optional
        the format of the input trajectory. Default is "PQ".
        For more information on the supported formats please visit
        :py:class:`~PQAnalysis.traj.formats.MDEngineFormat`.    
    mode : FileWritingMode | str, optional
        The writing mode. Default is "w".
        The writing mode can be either a string or a FileWritingMode enum value.
        Possible values are:
        - "w" or FileWritingMode.WRITE: write mode (default, no overwrite)
        - "a" or FileWritingMode.APPEND: append mode
        - "o" or FileWritingMode.OVERWRITE: overwrite mode
    """
    md_format = MDEngineFormat(md_format)

    input_reader = ThermalExpansionInputFileReader(input_file)
    input_reader.read()

    temperature_points = np.array(input_reader.temperature_points)

    box_data = []
    for i, box_file in enumerate(input_reader.box_files):
        print(
            f"Reading box file: {box_file} at temperature: {temperature_points[i]} K"
        )
        box_reader = BoxFileReader(filename=box_file, engine_format=md_format)
        box = box_reader.read()
        box_data.append(box)

    _thermal_expansion = ThermalExpansion(
        temperature_points=temperature_points, boxes=box_data
    )

    data_writer = ThermalExpansionDataWriter(
        filename=input_reader.out_file_key, mode=mode
    )

    log_writer = ThermalExpansionLogWriter(
        filename=input_reader.log_file, mode=mode
    )

    log_writer.write_before_run(_thermal_expansion)

    boxes_avg_data, boxes_std_data, thermal_expansion_data = _thermal_expansion.run()

    data_writer.write(
        temperature_points=temperature_points,
        boxes_avg_data=boxes_avg_data,
        boxes_std_data=boxes_std_data,
        thermal_expansion_data=thermal_expansion_data
    )
    log_writer.write_after_run(_thermal_expansion)
