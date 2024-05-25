"""
A module containing API functions to convert between different file formats.
"""

from beartype.typing import List

from PQAnalysis.core import Cell
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.type_checking import runtime_type_checking

from .traj_file import (
    TrajectoryWriter,
    TrajectoryReader,
)
from .box_writer import BoxWriter
from .formats import BoxFileFormat
from .gen_file import (
    write_gen_file,
    read_gen_file,
)
from .restart_file.api import read_restart_file
from .traj_file.api import write_trajectory



@runtime_type_checking
def gen2xyz(
    gen_file: str,
    output: str | None = None,
    md_format: MDEngineFormat | str = MDEngineFormat.PQ,
    print_box: bool = True,
    mode: FileWritingMode | str = "w"
) -> None:
    """
    Converts a gen file to a xyz file and prints it to stdout or writes it to a file.

    Parameters
    ----------
    gen_file : str
        The gen file to be converted.
    output : str | None
        The output file. If not specified, the output is printed to stdout.
    md_format : MDEngineFormat | str, optional
        The format of the md engine for the output file. The default is MDEngineFormat.PQ.
    print_box : bool, optional
        If True, the box is not printed. If False, the box is printed. Default is False.
    mode : FileWritingMode | str, optional
        The writing mode, by default "w". The following modes are available:
        - "w": write
        - "a": append
        - "o": overwrite
    """
    system = read_gen_file(gen_file)

    if not print_box:
        system.cell = Cell()

    write_trajectory(
        system,
        output,
        engine_format=md_format,
        traj_type="xyz",
        mode=mode,
    )



@runtime_type_checking
def xyz2gen(
    xyz_file: str,
    output: str | None = None,
    periodic: bool | None = None,
    mode: FileWritingMode | str = "w",
    md_format: MDEngineFormat | str = MDEngineFormat.PQ,
) -> None:
    """
    Converts a xyz file to a gen file and prints it to stdout or writes it to a file.

    Parameters
    ----------
    xyz_file : str
        The xyz file to be converted.
    output : str | None
        The output file. If not specified, the output is printed to stdout.
    periodic : bool | None, optional
        The periodicity of the system. If True, the system is considered periodic.
        If False, the system is considered non-periodic. If None, the periodicity 
        is inferred from the system, by default None.
    mode : FileWritingMode | str, optional
        The writing mode, by default "w". The following modes are available:
        - "w": write
        - "a": append
        - "o": overwrite
    md_format : MDEngineFormat | str, optional
        The format of the md engine for the output file. The default is MDEngineFormat.PQ.
    """

    system = TrajectoryReader(
        xyz_file,
        md_format=md_format,
        traj_format="xyz",
    ).read()

    write_gen_file(
        filename=output,
        system=system[0],
        periodic=periodic,
        mode=mode,
    )



@runtime_type_checking
def rst2xyz(
    restart_file: str,
    output: str | None = None,
    print_box: bool = True,
    md_format: MDEngineFormat | str = MDEngineFormat.PQ,
    mode: FileWritingMode | str = "w"
) -> None:
    """
    Converts a restart file to a xyz file and prints it to stdout or writes it to a file.

    When the print_box flag is set to True, the box is printed as well.
    This means that after the number of atoms the box is printed in the
    same line in the format a b c alpha beta gamma.

    Parameters
    ----------
    restart_file : str
        The restart file to be converted.
    output : str | None
        The output file. If not specified, the output is printed to stdout.
    print_box : bool
        If True, the box is printed. If False, the box is not printed. Default is True.
    md_format : MDEngineFormat | str, optional
        The format of the md engine for the output file. The default is MDEngineFormat.PQ.
    mode : FileWritingMode | str, optional
        The writing mode, by default "w". The following modes are available:
        - "w": write
        - "a": append
        - "o": overwrite
    """
    system = read_restart_file(restart_file)

    if not print_box:
        system.cell = Cell()

    write_trajectory(
        system, output, engine_format=md_format, traj_type="xyz", mode=mode
    )



@runtime_type_checking
def traj2box(
    trajectory_files: List[str],
    vmd: bool,
    output: str | None = None,
    mode: FileWritingMode | str = "w"
) -> None:
    """
    Converts multiple trajectory files to a box file and 
    prints it to stdout or writes it to a file.

    Without the vmd option the output is printed in a data file format.
    The first column represents the step starting from 1, the second to fourth column
    represent the box vectors a, b, c, the fifth to seventh column represent the box angles.

    With the vmd option the output is printed in a VMD file format. Meaning the output is
    in xyz format with 8 particle entries representing the vertices of the box. The comment
    line contains the information about the box dimensions a, b and c and the box angles.

    Parameters
    ----------
    trajectory_files : list of str
        The trajectory file(s) to be converted.
    vmd : bool
        Output in VMD format.
    output : str | None, optional
        The output file. If not specified, the output is printed to stdout. Default is None.
    mode : FileWritingMode | str, optional
        The writing mode, by default "w". The following modes are available:
        - "w": write
        - "a": append
        - "o": overwrite
    """

    if vmd:
        output_format = BoxFileFormat.VMD
    else:
        output_format = BoxFileFormat.DATA

    writer = BoxWriter(filename=output, output_format=output_format, mode=mode)

    for filename in trajectory_files:
        reader = TrajectoryReader(filename)
        trajectory = reader.read()

        writer.write(trajectory, reset_counter=False)



@runtime_type_checking
def traj2qmcfc(
    trajectory_files: List[str],
    output: str | None = None,
    mode: FileWritingMode | str = "w"
) -> None:
    """
    Converts multiple trajectory files from a PQ format to a 
    QMCFC format and prints it to stdout or writes it to a file.

    Parameters
    ----------
    trajectory_files : list of str
        The trajectory file(s) to be converted.
    output : str, optional
        The output file. If not specified, the output is printed to stdout.
    mode : FileWritingMode | str, optional
        The writing mode, by default "w". The following modes are available:
        - "w": write
        - "a": append
        - "o": overwrite
    """

    writer = TrajectoryWriter(
        filename=output, engine_format="qmcfc", mode=mode
    )

    for filename in trajectory_files:
        reader = TrajectoryReader(filename)
        trajectory = reader.read()

        writer.write(trajectory)
