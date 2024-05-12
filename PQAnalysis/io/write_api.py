"""
A module containing API functions for writing different objects to a file.
"""

import logging

from beartype.typing import Any

from PQAnalysis.traj import Trajectory
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.exceptions import PQNotImplementedError
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__

from .box_writer import BoxWriter
from .formats import FileWritingMode
from .traj_file.api import write_trajectory

logger = logging.getLogger(__package_name__).getChild("io.write_api")
logger = setup_logger(logger)



@runtime_type_checking
def write(
    object_to_write: Any,
    filename: str | None = None,
    mode: FileWritingMode | str = FileWritingMode.WRITE,
) -> None:
    """
    API write wrapper function for writing different objects to a file.

    It can call the following specialized write functions
    (depending on the object_to_write):

    write_trajectory: Writes a trajectory to a file.
        - Trajectory
        - AtomicSystem


    Parameters
    ----------
    object_to_write : Any
        The object to write to a file.
    filename : str | None, optional
        The name of the file to write to. If None, the output is printed to stdout.
    mode : FileWritingMode | str, optional
        The writing mode. If a string is given, it is converted to a FileWritingMode enum.
    """

    if isinstance(object_to_write, (Trajectory, AtomicSystem)):

        write_trajectory(  # TODO: still missing important kwargs
            object_to_write,
            filename,
            mode=mode
        )

    else:

        logger.error(
            (
            f"Writing object of type {type(object_to_write)} "
            "is not implemented yet."
            ),
            exception=PQNotImplementedError
        )



@runtime_type_checking
def write_box(
    traj: Trajectory,
    filename: str | None = None,
    output_format: str | None = None
) -> None:
    """
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
    """

    writer = BoxWriter(filename, output_format)
    writer.write(traj)
