from beartype.typing import List

from . import RestartFileReader, TrajectoryWriter, BoxWriter, TrajectoryReader, BoxFileFormat

from PQAnalysis.core import Cell
from PQAnalysis.traj import MDEngineFormat


def rst2xyz(restart_file: str,
            output: str | None = None,
            print_box: bool = True,
            md_format: MDEngineFormat | str = MDEngineFormat.PQ
            ):
    """
    Converts a restart file to a xyz file and prints it to stdout or writes it to a file.

    When the print_box flag is set to True, the box is printed as well. This means that after the number of atoms the box is printed in the same line in the format a b c alpha beta gamma.

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
    """
    reader = RestartFileReader(restart_file)
    frame = reader.read()

    if not print_box:
        frame.cell = Cell()

    writer = TrajectoryWriter(filename=output, format=md_format)
    writer.write(frame, type="xyz")


def traj2box(trajectory_files: List[str], vmd: bool, output: str | None = None) -> None:
    """
    Converts multiple trajectory files to a box file and prints it to stdout or writes it to a file.

    Without the vmd option the output is printed in a data file format.
    The first column represents the step starting from 1, the second to fourth column
    represent the box vectors a, b, c, the fifth to seventh column represent the box angles.

    With the vmd option the output is printed in a VMD file format. Meaning the output is
    in xyz format with 8 particle entries representing the vertices of the box. The comment
    line contains the information about the box dimensions a, b and c and the box angles.

    Parameters
    ----------
    trajectory_file : list of str
        The trajectory file(s) to be converted.
    vmd : bool
        Output in VMD format.
    output : str | None, optional
        The output file. If not specified, the output is printed to stdout. Default is None.
    """

    if vmd:
        output_format = BoxFileFormat.VMD
    else:
        output_format = BoxFileFormat.DATA

    writer = BoxWriter(filename=output, output_format=output_format)

    for filename in trajectory_files:
        reader = TrajectoryReader(filename)
        trajectory = reader.read()

        writer.write(trajectory, reset_counter=False)


def traj2qmcfc(trajectory_files: List[str], output: str | None = None):
    """
    Converts multiple trajectory files from a PQ format to a QMCFC format and prints it to stdout or writes it to a file.

    Parameters
    ----------
    trajectory_file : list of str
        The trajectory file(s) to be converted.
    output : str, optional
        The output file. If not specified, the output is printed to stdout.
    """

    writer = TrajectoryWriter(filename=output, format="qmcfc")

    for filename in trajectory_files:
        reader = TrajectoryReader(filename)
        trajectory = reader.read()

        writer.write(trajectory)
