"""
This module provides API functions for input/output handling of molecular dynamics simulations.
"""

from beartype.typing import List

from . import RestartFileReader, TrajectoryWriter, BoxWriter, TrajectoryReader, BoxFileFormat
from .inputFileReader import PIMD_QMCF_InputFileReader as Reader
from .inputFileReader.formats import InputFileFormat

from PQAnalysis.types import PositiveReal
from PQAnalysis.core import Cell
from PQAnalysis.traj import Trajectory, TrajectoryFormat, MDEngineFormat


def continue_input_file(input_file: str, n: PositiveReal = 1, input_format: InputFileFormat | str = InputFileFormat.PIMD_QMCF):
    """
    API function for continuing an input file.

    This function reads the input file and continues it n times. This means that the number in the filename is increased by one and all other numbers in the start- and output-file keys within the input file are increased by one as well.

    Parameters
    ----------
    input_file : str
        the path to the input file, which should be continued
    n : PositiveReal, optional
        the number of times the input file should be continued, by default 1
    input_format : InputFileFormat | str, optional
        the format of the input file, by default InputFileFormat.PIMD_QMCF

    Raises
    ------
    NotImplementedError
        if the input format is not pimd-qmcf
    """
    input_format = InputFileFormat(input_format)

    if input_format != InputFileFormat.PIMD_QMCF:
        raise NotImplementedError(
            f"Format {input_format} not implemented yet for continuing input file.")

    reader = Reader(input_file)
    reader.read()
    reader.continue_input_file(n)


def rst2xyz(restart_file: str, output: str | None = None, print_box: bool = True):
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
    """
    reader = RestartFileReader(restart_file)
    frame = reader.read()

    if not print_box:
        frame.cell = Cell()

    writer = TrajectoryWriter(filename=output)
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


def write_trajectory(traj,
                     filename: str | None = None,
                     format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF,
                     type: TrajectoryFormat | str = TrajectoryFormat.XYZ
                     ) -> None:
    """Wrapper for TrajectoryWriter to write a trajectory to a file.

    if format is None, the default PIMD-QMCF format is used. (see TrajectoryWriter.formats for available formats)
    if format is 'qmcfc', the QMCFC format is used (see TrajectoryWriter.formats for more information).

    Parameters
    ----------
    traj : Trajectory
        The trajectory to write.
    filename : str, optional
        The name of the file to write to. If None, the output is printed to stdout.
    format : MDEngineFormat | str, optional
        The format of the md engine for the output file. The default is MDEngineFormat.PIMD_QMCF.
    type : TrajectoryFormat | str, optional
        The type of the data to write to the file. Default is TrajectoryFormat.XYZ.

    """

    writer = TrajectoryWriter(filename, format=format)
    writer.write(traj, type=type)


def traj2qmcfc(trajectory_files: List[str], output: str | None = None):
    """
    Converts multiple trajectory files from a PIMD-QMCF format to a QMCFC format and prints it to stdout or writes it to a file.

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
