"""
This module provides API functions for the angular distribution function (ADF) analysis.
"""

from PQAnalysis.io import TrajectoryReader
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.type_checking import runtime_type_checking

from .adf import ADF
from .adf_input_file_reader import ADFInputFileReader
from .adf_output_file_writer import ADFDataWriter, ADFLogWriter



@runtime_type_checking
def adf(input_file: str, md_format: MDEngineFormat | str = MDEngineFormat.PQ):
    """
    Calculates the angular distribution function (ADF) using a given input file.

    This is just a wrapper function combining the underlying classes and functions.

    For more information on the input file keys please
    visit :py:mod:`~PQAnalysis.analysis.adf.adf_input_file_reader`.
    For more information on the exact calculation of the ADF please
    visit :py:class:`~PQAnalysis.analysis.adf.adf.ADF`.

    Parameters
    ----------
    input_file : str
        The input file. For more information on the input file
        keys please visit :py:mod:`~PQAnalysis.analysis.adf.adf_input_file_reader`.
    md_format : MDEngineFormat | str, optional
        the format of the input trajectory. Default is "PQ".
        For more information on the supported formats please visit
        :py:class:`~PQAnalysis.traj.formats.MDEngineFormat`.
    """

    md_format = MDEngineFormat(md_format)

    input_reader = ADFInputFileReader(input_file)
    input_reader.read()

    traj_reader = TrajectoryReader(
        input_reader.traj_files,
        md_format=md_format,
    )

    _adf = ADF(
        traj=traj_reader,
        reference_species=input_reader.reference_selection,
        target_species=input_reader.target_selection,
        target_species_2=input_reader.target_selection_2,
        use_full_atom_info=input_reader.use_full_atom_info,
        n_angle_bins=input_reader.n_angle_bins,
        delta_angle=input_reader.delta_angle,
        r_min_1=input_reader.r_min_1,
        r_max_1=input_reader.r_max_1,
        r_min_2=input_reader.r_min_2,
        r_max_2=input_reader.r_max_2,
    )

    data_writer = ADFDataWriter(input_reader.out_file)
    log_writer = ADFLogWriter(input_reader.log_file)
    log_writer.write_before_run(_adf)

    adf_data = _adf.run()

    data_writer.write(adf_data)
    log_writer.write_after_run(_adf)
