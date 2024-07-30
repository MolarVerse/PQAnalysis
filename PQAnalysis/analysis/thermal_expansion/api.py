"""
This module provides API functions for
linear or volumetric thermal expansion coefficient analysis.
"""


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
    """
    md_format = MDEngineFormat(md_format)

    input_reader = ThermalExpansionInputFileReader(input_file)
    input_reader.read()

    temperature_points = input_reader.temperature_points

    box_data = []
    for box_file in input_reader.box_files:
        box_reader = BoxFileReader(
            filename=box_file,
            engine_format=md_format,
            unit=input_reader.unit
        )
        box = box_reader.read()
        box_data.append(box)

    _thermal_expansion = ThermalExpansion(
        temperature_points=temperature_points,
        boxes=box_data
    )

    data_writer = ThermalExpansionDataWriter(
        filename=input_reader.out_file_key
    )

    log_writer = ThermalExpansionLogWriter(filename=input_reader.log_file)

    log_writer.write_before_run(_thermal_expansion)

    box_avg_data, box_std_data, thermal_expansion_data = _thermal_expansion.run()

    data_writer.write(
        temperature_points=temperature_points,
        box_avg_data=box_avg_data,
        box_std_data=box_std_data,
        thermal_expansion_data=thermal_expansion_data
    )
    log_writer.write_after_run(_thermal_expansion)
