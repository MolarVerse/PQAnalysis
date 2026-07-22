"""
This module provides API functions for the mean square displacement (MSD) analysis.
"""

from PQAnalysis.io import TrajectoryReader
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.type_checking import runtime_type_checking

from .msd import MSD
from .msd_input_file_reader import MSDInputFileReader
from .msd_output_file_writer import MSDDataWriter, MSDLogWriter



@runtime_type_checking
def msd(input_file: str, md_format: MDEngineFormat | str = MDEngineFormat.PQ):
    """
    Calculates the mean square displacement (MSD) using a given input file.

    This is just a wrapper function combining the underlying classes and functions.

    For more information on the input file keys please
    visit :py:mod:`~PQAnalysis.analysis.msd.msd_input_file_reader`.
    For more information on the exact calculation of
    the MSD please visit :py:class:`~PQAnalysis.analysis.msd.msd.MSD`.

    Parameters
    ----------
    input_file : str
        The input file. For more information on the input file
        keys please visit :py:mod:`~PQAnalysis.analysis.msd.msd_input_file_reader`.
    md_format : MDEngineFormat | str, optional
        the format of the input trajectory. Default is "PQ".
        For more information on the supported formats please visit
        :py:class:`~PQAnalysis.traj.formats.MDEngineFormat`.
    """

    md_format = MDEngineFormat(md_format)

    input_reader = MSDInputFileReader(input_file)
    input_reader.read()

    traj_reader = TrajectoryReader(
        input_reader.traj_files,
        md_format=md_format
    )

    _msd = MSD(
        traj=traj_reader,
        target_species=input_reader.target_selection,
        use_full_atom_info=input_reader.use_full_atom_info,
        window=input_reader.window,
        gap=input_reader.gap,
        n_start=input_reader.n_start,
        time_step=input_reader.time_step,
        fit_window=input_reader.fit_window,
    )

    data_writer = MSDDataWriter(input_reader.out_file)
    log_writer = MSDLogWriter(input_reader.log_file)
    log_writer.write_before_run(_msd)

    msd_data = _msd.run()

    data_writer.write(msd_data)
    log_writer.write_after_run(_msd)
