from .genFileReader import GenFileReader
from .genFileWriter import GenFileWriter
from PQAnalysis.atomicSystem import AtomicSystem


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


def write_gen_file(filename: str,
                   system: AtomicSystem,
                   periodic: bool | None = None
                   ) -> None:
    """
    Function to write a gen file.

    Parameters
    ----------
    filename : str
        The filename of the gen file.
    system : AtomicSystem
        The system to write.
    periodic : bool, optional
        The periodicity of the system. If True, the system is considered periodic. If False, the system is considered non-periodic. If None, the periodicity is inferred from the system, by default None.
    """

    GenFileWriter(filename).write(system, periodic)
