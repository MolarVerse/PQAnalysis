"""
A module containing the TrajectoryWriter class and its associated methods.

...

Classes
-------
TrajectoryWriter
    A class for writing a trajectory to a file.
"""

import numpy as np

from beartype.typing import List

from .base import BaseWriter
from ..traj.trajectory import Trajectory
from ..traj.formats import TrajectoryFormat, MDEngineFormat
from ..traj.frame import Frame
from ..core.cell import Cell
from ..core.atom import Atom
from ..types import Numpy2DFloatArray, Numpy1DFloatArray


def write_trajectory(traj,
                     filename: str | None = None,
                     format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF,
                     type: TrajectoryFormat | str = TrajectoryFormat.XYZ
                     ) -> None:
    """
    Wrapper for TrajectoryWriter to write a trajectory to a file.

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


class TrajectoryWriter(BaseWriter):
    """
    A class for writing a trajectory to a file.
    Inherits from BaseWriter. See BaseWriter for more information.

    It can write a trajectory to a file in either a PIMD-QMCF format or a QMCFC format.

    ...

    Class Attributes
    ----------------
    formats : list of str
        The available formats for the trajectory file.

            #TODO: put this description into formats!!!
            PIMD-QMCF format for one frame:
                header line containing the number of atoms and the cell information (if available)
                arbitrary comment line
                coordinates of the atoms in the format 'element x y z'

            QMCFC format for one frame:
                header line containing the number of atoms and the cell information (if available)
                arbitrary comment line
                X 0.0 0.0 0.0
                coordinates of the atoms in the format 'element x y z'

    _type : TrajectoryFormat
        The type of the data to write to the file. Default is TrajectoryFormat.XYZ.

    Attributes
    ----------
    format : MDEngineFormat
        The format of the md engine for the output file. The default is MDEngineFormat.PIMD_QMCF.
    """

    _format: MDEngineFormat
    _type: TrajectoryFormat = TrajectoryFormat.XYZ

    def __init__(self,
                 filename: str | None = None,
                 format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF,
                 mode: str = 'w'
                 ) -> None:
        """
        It sets the file to write to - either a file or stdout (if filename is None) - and the mode of the file.

        Parameters
        ----------
        filename : str, optional
            The name of the file to write to. If None, the output is printed to stdout.
        format : MDEngineFormat | str, optional
            The format of the md engine for the output file. The default is MDEngineFormat.PIMD_QMCF.
        mode : str, optional
            The mode of the file. Either 'w' for write or 'a' for append.
        """

        super().__init__(filename, mode)

        self.format = MDEngineFormat(format)

    def write(self, trajectory: Trajectory, type: TrajectoryFormat | str = TrajectoryFormat.XYZ) -> None:
        """
        Writes the trajectory to the file.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to write.
        """
        self._type = TrajectoryFormat(type)
        if self._type == TrajectoryFormat.XYZ:
            self.write_positions(trajectory)
        elif self._type == TrajectoryFormat.VEL:
            self.write_velocities(trajectory)
        elif self._type == TrajectoryFormat.FORCE:
            self.write_forces(trajectory)
        elif self._type == TrajectoryFormat.CHARGE:
            self.write_charges(trajectory)

        self.close()

    def write_positions(self, trajectory: Trajectory) -> None:
        """
        Writes the positions of the trajectory to the file.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to write.
        """
        self._type = TrajectoryFormat.XYZ
        self.open()
        for frame in trajectory:
            self._write_header(frame.n_atoms, frame.cell)
            self._write_comment(frame)
            self._write_xyz(frame.pos, frame.atoms)

        self.close()

    def write_velocities(self, trajectory: Trajectory) -> None:
        """
        Writes the velocities of the trajectory to the file.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to write.
        """
        self._type = TrajectoryFormat.VEL
        self.open()
        for frame in trajectory:
            self._write_header(frame.n_atoms, frame.cell)
            self._write_comment(frame)
            self._write_xyz(frame.vel, frame.atoms)

        self.close()

    def write_forces(self, trajectory: Trajectory) -> None:
        """
        Writes the forces of the trajectory to the file.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to write.
        """
        self._type = TrajectoryFormat.FORCE
        self.open()
        for frame in trajectory:
            self._write_header(frame.n_atoms, frame.cell)
            self._write_comment(frame)
            self._write_xyz(frame.forces, frame.atoms)

        self.close()

    def write_charges(self, trajectory: Trajectory) -> None:
        """
        Writes the charges of the trajectory to the file.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to write.
        """
        self._type = TrajectoryFormat.CHARGE
        self.open()
        for frame in trajectory:
            self._write_header(frame.n_atoms, frame.cell)
            self._write_comment(frame)
            self._write_scalar(frame.charges, frame.atoms)

        self.close()

    def _write_header(self, n_atoms: int, cell: Cell | None = None) -> None:
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
                f"{n_atoms} {cell.x} {cell.y} {cell.z} {cell.alpha} {cell.beta} {cell.gamma}", file=self.file)
        else:
            print(f"{n_atoms}", file=self.file)

    def _write_comment(self, frame: Frame) -> None:
        """
        Writes the comment line of the frame to the file.

        Parameters
        ----------
        frame : Frame
            The frame to write the comment line of.
        """

        if self._type == TrajectoryFormat.FORCE:
            sum_forces = sum(frame.forces)
            print(
                f"sum of forces: {sum_forces[0]} {sum_forces[1]} {sum_forces[2]}", file=self.file)
        else:
            print("", file=self.file)

    def _write_xyz(self, xyz: Numpy2DFloatArray, atoms: List[Atom]) -> None:
        """
        Writes the xyz of the frame to the file.

        If format is 'qmcfc', an additional X 0.0 0.0 0.0 line is written.

        Parameters
        ----------
        xyz : np.array
            The xyz data of the atoms (either positions, velocities or forces).
        atoms : Elements
            The elements of the frame.
        """

        if self.format == MDEngineFormat.QMCFC and self._type == TrajectoryFormat.XYZ:
            print("X   0.0 0.0 0.0", file=self.file)

        for i in range(len(atoms)):
            print(
                f"{atoms[i].name} {xyz[i][0]} {xyz[i][1]} {xyz[i][2]}", file=self.file)

    def _write_scalar(self, scalar: Numpy1DFloatArray, atoms: List[Atom]) -> None:
        """
        Writes the charges of the frame to the file.

        Parameters
        ----------
        scalar : np.array
            scalar data of the atoms (atm only charges).
        atoms : Elements
            The elements of the frame.
        """

        for i in range(len(atoms)):
            print(
                f"{atoms[i].name} {scalar[i]}", file=self.file)

    @property
    def format(self) -> MDEngineFormat:
        return self._format

    @format.setter
    def format(self, format: MDEngineFormat | str) -> None:
        self._format = MDEngineFormat(format)
