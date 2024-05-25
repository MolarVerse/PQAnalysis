"""
A module containing different functions to read and write .gen files.
"""

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.type_checking import runtime_type_checking

from .gen_file_reader import GenFileReader
from .gen_file_writer import GenFileWriter



@runtime_type_checking
def read_gen_file(filename: str) -> AtomicSystem:
    """
    Function to read a gen file.

    Parameters
    ----------
    filename : str
        The filename of the gen file.

    Returns
    -------
    AtomicSystem
        The AtomicSystem including the Cell object.
    """

    return GenFileReader(filename).read()



@runtime_type_checking
def write_gen_file(
    system: AtomicSystem,
    filename: str | None = None,
    periodic: bool | None = None,
    mode: FileWritingMode | str = "w"
) -> None:
    """
    Function to write a gen file.

    Parameters
    ----------
    system : AtomicSystem
        The system to write.
    filename : str | None, optional
        The filename of the gen file. If None, the file is not opened, 
        by default None.
    periodic : bool, optional
        The periodicity of the system. If True, the system is considered periodic. 
        If False, the system is considered non-periodic. If None, the periodicity
        is inferred from the system, by default None.
    mode : FileWritingMode | str, optional
        The writing mode, by default "w". The following modes are available:
        - "w": write
        - "a": append
        - "o": overwrite
    """

    GenFileWriter(filename, mode=mode).write(system, periodic)
