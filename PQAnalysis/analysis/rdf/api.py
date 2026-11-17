"""
This module provides API functions for the radial distribution function (RDF) analysis.
"""

from beartype.typing import Sequence

from PQAnalysis.io import TrajectoryReader, RestartFileReader, MoldescriptorReader
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.topology import Topology
from PQAnalysis.type_checking import runtime_type_checking

from .rdf import RDF
from .rdf_input_file_reader import RDFInputFileReader
from .rdf_output_file_writer import RDFDataWriter, RDFLogWriter



@runtime_type_checking
def rdf(
    input_file: str,
    md_format: MDEngineFormat | str = MDEngineFormat.PQ,
    export_files: Sequence[str] | None = None,
    mode: str | FileWritingMode = FileWritingMode.WRITE,
):
    """
    Calculates the radial distribution function (RDF) using a given input file.

    This is just a wrapper function combining the underlying classes and functions.

    For more information on the input file keys please
    visit :py:mod:`~PQAnalysis.analysis.rdf.rdfInputFileReader`.
    For more information on the exact calculation of
    the RDF please visit :py:class:`~PQAnalysis.analysis.rdf.rdf.RDF`.

    Parameters
    ----------
    input_file : str
        The input file. For more information on the input file
        keys please visit :py:mod:`~PQAnalysis.analysis.rdf.rdfInputFileReader`.
    md_format : MDEngineFormat | str, optional
        the format of the input trajectory. Default is "PQ".
        For more information on the supported formats please visit
        :py:class:`~PQAnalysis.traj.formats.MDEngineFormat`.
    export_files : Sequence[str] | None, optional
        Additional table outputs. ``.csv``, ``.tsv`` and ``.xvg`` select
        those formats; other extensions use native PQAnalysis text.
    mode : str | FileWritingMode, optional
        The writing mode for all output files: "w" stops if a file
        exists, "a" appends, "o" overwrites, by default "w".
    """

    md_format = MDEngineFormat(md_format)

    input_reader = RDFInputFileReader(input_file)
    input_reader.read()

    # writers check for existing output files, so fail before the slow setup
    data_writer = RDFDataWriter(
        input_reader.out_file,
        export_files=export_files,
        mode=mode,
    )
    log_writer = RDFLogWriter(input_reader.log_file, mode=mode)

    if input_reader.restart_file is not None:
        restart_reader = RestartFileReader(input_reader.restart_file)
        restart_frame = restart_reader.read()
        topology = restart_frame.topology
    else:
        topology = None

    if input_reader.moldescriptor_file is not None:
        moldescriptor_reader = MoldescriptorReader(
            input_reader.moldescriptor_file
        )
        reference_residues = moldescriptor_reader.read()
    else:
        reference_residues = None

    if topology is not None:
        topology = Topology(
            atoms=topology.atoms,
            residue_ids=topology.residue_ids,
            reference_residues=reference_residues
        )

    traj_reader = TrajectoryReader(
        input_reader.traj_files, md_format=md_format, topology=topology
    )

    _rdf = RDF(
        traj=traj_reader,
        reference_species=input_reader.reference_selection,
        target_species=input_reader.target_selection,
        use_full_atom_info=input_reader.use_full_atom_info,
        no_intra_molecular=input_reader.no_intra_molecular,
        reference_com=bool(input_reader.reference_com),
        target_com=bool(input_reader.target_com),
        n_bins=input_reader.n_bins,
        delta_r=input_reader.delta_r,
        r_max=input_reader.r_max,
        r_min=input_reader.r_min,
    )

    log_writer.write_before_run(_rdf)

    rdf_data = _rdf.run()

    data_writer.write(rdf_data)
    log_writer.write_after_run(_rdf)
