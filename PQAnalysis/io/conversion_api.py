"""
A module containing API functions to convert between different file formats.
"""

import numpy as np

from beartype.typing import List

from . import (
    TrajectoryWriter,
    BoxWriter,
    TrajectoryReader,
    BoxFileFormat,
    write_gen_file,
    read_gen_file,
    read_restart_file,
    write_trajectory,
)

from PQAnalysis.core import Cell
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.io.formats import FileWritingMode


def gen2xyz(gen_file: str,
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
    no_box : bool, optional
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

    write_trajectory(system, output, format=md_format, type="xyz", mode=mode)


def xyz2gen(xyz_file: str,
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
        The periodicity of the system. If True, the system is considered periodic. If False, the system is considered non-periodic. If None, the periodicity is inferred from the system, by default None.
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
        traj_format="xyz"
    ).read()

    write_gen_file(output, system, periodic, mode)


def rst2xyz(restart_file: str,
            output: str | None = None,
            print_box: bool = True,
            md_format: MDEngineFormat | str = MDEngineFormat.PQ,
            mode: FileWritingMode | str = "w"
            ) -> None:
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
    mode : FileWritingMode | str, optional
        The writing mode, by default "w". The following modes are available:
        - "w": write
        - "a": append
        - "o": overwrite
    """
    system = read_restart_file(restart_file)

    if not print_box:
        system.cell = Cell()

    write_trajectory(system, output, format=md_format, type="xyz", mode=mode)


def traj2box(trajectory_files: List[str],
             vmd: bool,
             output: str | None = None,
             mode: FileWritingMode | str = "w"
             ) -> None:
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


def traj2qmcfc(trajectory_files: List[str],
               output: str | None = None,
               mode: FileWritingMode | str = "w"
               ) -> None:
    """
    Converts multiple trajectory files from a PQ format to a QMCFC format and prints it to stdout or writes it to a file.

    Parameters
    ----------
    trajectory_file : list of str
        The trajectory file(s) to be converted.
    output : str, optional
        The output file. If not specified, the output is printed to stdout.
    mode : FileWritingMode | str, optional
        The writing mode, by default "w". The following modes are available:
        - "w": write
        - "a": append
        - "o": overwrite
    """

    writer = TrajectoryWriter(filename=output, format="qmcfc", mode=mode)

    for filename in trajectory_files:
        reader = TrajectoryReader(filename)
        trajectory = reader.read()

        writer.write(trajectory)


def create_nep_traj_files(file_prefixes: List[str] = None,
                          trajectory_files: List[str] | str = None,
                          force_files: List[str] | str = None,
                          stress_files: List[str] | str = None,
                          virial_files: List[str] | str = None,
                          energy_files: List[str] | str = None,
                          info_files: List[str] | str = None,
                          ):
    # convert all files to list

    trajectory_files = list(np.atleast_1d(trajectory_files))
    force_files = list(np.atleast_1d(force_files))

    if stress_files is not None:
        stress_files = list(np.atleast_1d(stress_files))

    if virial_files is not None:
        virial_files = list(np.atleast_1d(virial_files))

    if energy_files is not None:
        energy_files = list(np.atleast_1d(energy_files))

    if info_files is not None:
        info_files = list(np.atleast_1d(info_files))

    # check if len of all files is the same
