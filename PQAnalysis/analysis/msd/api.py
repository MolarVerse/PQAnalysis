"""
This module provides API functions for the radial distribution function (MSD) analysis.
"""

from .msd import MSD
from .msdInputFileReader import MSDInputFileReader
from .msdOutputFileWriter import MSDDataWriter, MSDLogWriter
from PQAnalysis.io import TrajectoryReader, RestartFileReader, MoldescriptorReader
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.topology import Topology


def msd(input_file: str, md_format: MDEngineFormat | str = MDEngineFormat.PQ):
    """
    Calculates the self-diffusion coefficient function (MSD) using a given input file.

    This is just a wrapper function combining the underlying classes and functions.

    For more information on the input file keys please visit :py:mod:`~PQAnalysis.analysis.msd.msdInputFileReader`.
    For more information on the exact calculation of the MSD please visit :py:class:`~PQAnalysis.analysis.msd.msd.MSD`.

    Parameters
    ----------
    input_file : str
        The input file. For more information on the input file keys please visit :py:mod:`~PQAnalysis.analysis.msd.msdInputFileReader`.
    md_format : MDEngineFormat | str, optional
        the format of the input trajectory. Default is "PQ". For more information on the supported formats please visit :py:class:`~PQAnalysis.traj.formats.MDEngineFormat`.
    """
    md_format = MDEngineFormat(md_format)

    input_reader = MSDInputFileReader(input_file)
    input_reader.read()

    if input_reader.restart_file is not None:
        restart_reader = RestartFileReader(input_reader.restart_file)
        restart_frame = restart_reader.read()
        topology = restart_frame.topology
    else:
        topology = None

    if input_reader.moldescriptor_file is not None:
        moldescriptor_reader = MoldescriptorReader(
            input_reader.moldescriptor_file)
        reference_residues = moldescriptor_reader.read()
    else:
        reference_residues = None

    if topology is not None:
        topology = Topology(atoms=topology.atoms, residue_ids=topology.residue_ids,
                            reference_residues=reference_residues)

    traj_reader = TrajectoryReader(
        input_reader.traj_files, md_format=md_format, topology=topology)

    msd = MSD(
        traj=traj_reader,
        reference_species=input_reader.reference_selection,
        target_species=input_reader.target_selection,
        use_full_atom_info=input_reader.use_full_atom_info,
        no_intra_molecular=input_reader.no_intra_molecular,
        n_bins=input_reader.n_bins,
        delta_r=input_reader.delta_r,
        r_max=input_reader.r_max,
        r_min=input_reader.r_min,
    )

    data_writer = MSDDataWriter(input_reader.out_file)
    log_writer = MSDLogWriter(input_reader.log_file)
    log_writer.write_before_run(msd)

    msd_data = msd.run()

    data_writer.write(msd_data)
    log_writer.write_after_run(msd)
