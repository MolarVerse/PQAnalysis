"""
A module containing different API functions for reading and writing restart files.
"""

from PQAnalysis.core import Residues
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.type_checking import runtime_type_checking

from .restart_reader import RestartFileReader
from .restart_writer import RestartFileWriter


@runtime_type_checking
def read_restart_file(
    filename: str,
    moldescriptor_filename: str | None = None,
    reference_residues: Residues | None = None,
    md_engine_format: MDEngineFormat | str = MDEngineFormat.PQ
) -> AtomicSystem:
    """
    API function for reading a restart file.

    Parameters
    ----------
    filename : str
        The filename of the restart file.
    moldescriptor_filename : str | None, optional
        The filename of the moldescriptor file that is read by the MoldescriptorReader 
        to obtain the reference residues of the system,
        by default None
    reference_residues : Residues | None, optional
        The reference residues of the system, in general these are obtained by the 
        MoldescriptorReader - only used if moldescriptor_filename is None,
        by default None
    md_engine_format : MDEngineFormat | str, optional
        The format of the restart file,
        by default MDEngineFormat.PQ

    Returns
    -------
    AtomicSystem
        The AtomicSystem object including the Topology with the molecular types.
    """

    reader = RestartFileReader(
        filename=filename,
        moldescriptor_filename=moldescriptor_filename,
        reference_residues=reference_residues,
        md_engine_format=md_engine_format
    )

    return reader.read()


@runtime_type_checking
def write_restart_file(
    atomic_system: AtomicSystem,
    filename: str | None = None,
    md_engine_format: MDEngineFormat | str = MDEngineFormat.PQ,
    mode: FileWritingMode | str = 'w'
) -> None:
    """
    API function for reading a restart file.

    Parameters
    ----------
    atomic_system : AtomicSystem
        The AtomicSystem object to write.
    filename : str
        The filename of the restart file.
    md_engine_format : MDEngineFormat | str, optional
        The format of the restart file,
        by default MDEngineFormat.PQ
    mode : FileWritingMode | str, optional
        The mode of the file. Either 'w' for write, 
        'a' for append or 'o' for overwrite. The default is 'w'.
    """

    writer = RestartFileWriter(
        filename=filename,
        md_engine_format=md_engine_format,
        mode=mode
    )

    writer.write(atomic_system)
