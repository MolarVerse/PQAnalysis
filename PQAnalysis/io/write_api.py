from beartype.typing import Any

from . import BoxWriter, FileWritingMode
from .traj_file.api import write_trajectory

from PQAnalysis.traj import Trajectory
from PQAnalysis.atomicSystem import AtomicSystem


def write(object_to_write: Any,
          filename: str | None = None,
          mode: FileWritingMode | str = FileWritingMode.WRITE,
          **kwargs,
          ) -> None:
    """
    API write wrapper function for writing different objects to a file.

    It can call the following specialized write functions (depending on the object_to_write):

    write_trajectory: Writes a trajectory to a file.
        - Trajectory
        - AtomicSystem


    Parameters
    ----------
    object_to_write : Any
        _description_
    filename : str | None, optional
        _description_, by default None
    mode : FileWritingMode | str, optional
        _description_, by default FileWritingMode.WRITE
    kwargs : dict
        kwargs dictionary which is passed to the specialized write function for the object.
    """

    if isinstance(object_to_write, Trajectory):
        write_trajectory(object_to_write, filename, format, type, mode)
    elif isinstance(object_to_write, AtomicSystem):
        write_trajectory(Trajectory(object_to_write),
                         filename, format, type, mode)
    else:
        raise NotImplementedError(
            f"Writing object of type {type(object_to_write)} is not implemented yet.")


def write_box(traj: Trajectory,
              filename: str | None = None,
              output_format: str | None = None
              ) -> None:
    '''
    Writes the given trajectory to the file in a selected box file format.



    Parameters
    ----------
    traj : Trajectory
        The trajectory to write.
    filename : str, optional
        The name of the file to write to. If None, the output is printed to stdout.
    output_format : str, optional
        The format of the file. If None, the format is inferred as a data file format.
        (see BoxWriter.formats for available formats)
    '''

    writer = BoxWriter(filename, output_format)
    writer.write(traj)
