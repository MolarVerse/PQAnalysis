"""
A module containing the TrajectoryWriter class and its associated methods.

...

Classes
-------
TrajectoryWriter
    A class for writing a trajectory to a file.

Methods
-------
write_trajectory(traj, filename=None, format=None)
    Wrapper for TrajectoryWriter to write a trajectory to a file.
"""

import numpy as np

from typing import Union, List

from PQAnalysis.io.base import BaseWriter
from PQAnalysis.traj.trajectory import Trajectory
from PQAnalysis.pbc.cell import Cell


def write_trajectory(traj, filename: str = None, format: str = None):
    '''
    Wrapper for TrajectoryWriter to write a trajectory to a file.

    if format is None, the default PIMD-QMCF format is used. (see TrajectoryWriter.formats for available formats)
    if format is 'qmcfc', the QMCFC format is used (see TrajectoryWriter.formats for more information).

    Parameters
    ----------
    traj : Trajectory
        The trajectory to write.
    filename : str, optional
        The name of the file to write to. If None, the output is printed to stdout.
    '''

    writer = TrajectoryWriter(filename, format)
    writer.write(traj)


class TrajectoryWriter(BaseWriter):
    """
    A class for writing a trajectory to a file.
    Inherits from BaseWriter. See BaseWriter for more information.

    It can write a trajectory to a file in either a PIMD-QMCF format or a QMCFC format.

    ...

    Attributes
    ----------
    formats : list of str
        The available formats for the trajectory file.

            PIMD-QMCF format for one frame:
                header line containing the number of atoms and the cell information (if available)
                arbitrary comment line
                coordinates of the atoms in the format 'element x y z'

            QMCFC format for one frame:
                header line containing the number of atoms and the cell information (if available)
                arbitrary comment line
                X 0.0 0.0 0.0
                coordinates of the atoms in the format 'element x y z'

    format : str
        The format of the file.

    Methods
    -------
    write(traj)
        Writes the trajectory to the file.
    """

    formats = [None, 'qmcfc']

    def __init__(self, filename: str = None, format: str = None, mode: str = 'w'):
        """
        It sets the file to write to - either a file or stdout (if filename is None) - and the mode of the file.

        Parameters
        ----------
        filename : str, optional
            The name of the file to write to. If None, the output is printed to stdout.
        format : str, optional
            The format of the file. If None, the default PIMD-QMCF format is used.
            (see TrajectoryWriter.formats for available formats)
        mode : str, optional
            The mode of the file. Either 'w' for write or 'a' for append.
        """

        super().__init__(filename, mode)
        if format not in self.formats:
            raise ValueError(
                'Invalid format. Has to be either \'qmcfc\' or \'None\'.')

        self.format = format

    def write(self, trajectory: Trajectory):
        """
        Writes the trajectory to the file.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to write.
        """
        self.open()
        for frame in trajectory:
            self.__write_header__(frame.n_atoms, frame.cell)
            self.__write_coordinates__(frame.xyz, frame.atoms)
        self.close()

    def __write_header__(self, n_atoms: int, cell: Cell):
        """
        Writes the header line of the frame to the file.

        Parameters
        ----------
        n_atoms : int
            The number of atoms in the frame.
        cell : Cell
            The cell of the frame. If None, only the number of atoms is written.
        """

        if cell is not None:
            print(
                f"{n_atoms} {cell.x} {cell.y} {cell.z} {cell.alpha} {cell.beta} {cell.gamma}\n", file=self.file)
        else:
            print(f"{n_atoms}\n", file=self.file)

    def __write_coordinates__(self, xyz: np.array, atoms: List[str]):
        """
        Writes the coordinates of the frame to the file.

        If format is 'qmcfc', an additional X 0.0 0.0 0.0 line is written.

        Parameters
        ----------
        xyz : np.array
            The xyz coordinates of the atoms.
        atoms : list of str
            The names of the atoms.
        """

        if self.format == "qmcfc":
            print("X   0.0 0.0 0.0", file=self.file)

        for i in range(len(atoms)):
            print(
                f"{atoms[i]} {xyz[i][0]} {xyz[i][1]} {xyz[i][2]}", file=self.file)
