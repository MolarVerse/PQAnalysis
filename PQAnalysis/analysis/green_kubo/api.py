"""
This module provides API functions for the Green-Kubo transport
coefficient analysis.
"""

from PQAnalysis.io import TrajectoryReader
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.type_checking import runtime_type_checking

from .green_kubo import GreenKubo
from .green_kubo_input_file_reader import GreenKuboInputFileReader
from .green_kubo_output_file_writer import (
    GreenKuboDataWriter,
    GreenKuboLogWriter,
)



@runtime_type_checking
def green_kubo(
    input_file: str,
    md_format: MDEngineFormat | str = MDEngineFormat.PQ,
):
    """
    Calculates a Green-Kubo transport coefficient (currently the
    self-diffusion coefficient) using a given input file.

    This is just a wrapper function combining the underlying classes
    and functions.

    For more information on the input file keys please visit
    :py:mod:`~PQAnalysis.analysis.green_kubo.green_kubo_input_file_reader`.
    For more information on the exact calculation of the Green-Kubo
    diffusion coefficient please visit
    :py:class:`~PQAnalysis.analysis.green_kubo.green_kubo.GreenKubo`.

    Parameters
    ----------
    input_file : str
        The input file. For more information on the input file keys
        please visit
        :py:mod:`~PQAnalysis.analysis.green_kubo.green_kubo_input_file_reader`.
    md_format : MDEngineFormat | str, optional
        the format of the input trajectory. Default is "PQ".
        For more information on the supported formats please visit
        :py:class:`~PQAnalysis.traj.formats.MDEngineFormat`.
    """

    md_format = MDEngineFormat(md_format)

    input_reader = GreenKuboInputFileReader(input_file)
    input_reader.read()

    traj_reader = TrajectoryReader(
        input_reader.traj_files,
        md_format=md_format,
    )

    _green_kubo = GreenKubo(
        traj=traj_reader,
        time_step=input_reader.time_step,
        window_size=input_reader.window,
        target_species=input_reader.target_selection,
        gap=input_reader.gap,
        fit_start=input_reader.fit_start,
        fit_stop=input_reader.fit_stop,
        method=input_reader.method,
        coefficient=input_reader.coefficient,
        n_blocks=input_reader.n_blocks,
        use_full_atom_info=input_reader.use_full_atom_info,
    )

    # the output writers are constructed before the run so that a
    # pre-existing output file aborts before the expensive analysis
    data_writer = GreenKuboDataWriter(input_reader.out_file)
    log_writer = GreenKuboLogWriter(input_reader.log_file)

    log_writer.write_before_run(_green_kubo)

    data = _green_kubo.run()

    data_writer.write(data)
    log_writer.write_after_run(_green_kubo)
