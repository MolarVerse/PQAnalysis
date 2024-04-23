"""
This module provides API functions for the radial distribution function (Diffcalc) analysis.
"""

from .diffcalc import Diffcalc
from .diffcalcInputFileReader import DiffcalcInputFileReader
from .diffcalcOutputFileWriter import DiffcalcDataWriter, DiffcalcLogWriter
from PQAnalysis.io import TrajectoryReader, RestartFileReader, MoldescriptorReader
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.topology import Topology


def diffcalc(input_file: str, md_format: MDEngineFormat | str = MDEngineFormat.PQ):
    """
    Calculates the self-diffusion coefficient function (Diffcalc) using a given input file.

    This is just a wrapper function combining the underlying classes and functions.

    For more information on the input file keys please visit :py:mod:`~PQAnalysis.analysis.diffcalc.diffcalcInputFileReader`.
    For more information on the exact calculation of the Diffcalc please visit :py:class:`~PQAnalysis.analysis.diffcalc.diffcalc.Diffcalc`.

    Parameters
    ----------
    input_file : str
        The input file. For more information on the input file keys please visit :py:mod:`~PQAnalysis.analysis.diffcalc.diffcalcInputFileReader`.
    md_format : MDEngineFormat | str, optional
        the format of the input trajectory. Default is "PQ". For more information on the supported formats please visit :py:class:`~PQAnalysis.traj.formats.MDEngineFormat`.
    """
    md_format = MDEngineFormat(md_format)

    input_reader = DiffcalcInputFileReader(input_file)
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

    diffcalc = Diffcalc(
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

    data_writer = DiffcalcDataWriter(input_reader.out_file)
    log_writer = DiffcalcLogWriter(input_reader.log_file)
    log_writer.write_before_run(diffcalc)

    diffcalc_data = diffcalc.run()

    data_writer.write(diffcalc_data)
    log_writer.write_after_run(diffcalc)
