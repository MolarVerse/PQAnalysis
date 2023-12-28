"""
A module for writing restart files.

Classes
-------
RestartFileWriter
    A class for writing restart files.
"""

import numpy as np

from . import BaseWriter, RestartFileWriterError
from ..traj import MDEngineFormat, Frame
from ..core import Cell


class RestartFileWriter(BaseWriter):
    """
    A class for writing restart files.

    Inherits from the BaseWriter class.
    """

    def __init__(self,
                 filename: str | None = None,
                 format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF,
                 mode: str = 'w'
                 ) -> None:
        """
        Initializes the RestartFileWriter with the given parameters.

        Parameters
        ----------
        filename : str | None, optional
            The filename of the restart file, by default None (stdout)
        format : MDEngineFormat | str, optional
            The format of the restart file, by default MDEngineFormat.PIMD_QMCF
        mode : str, optional
            The mode of the file, by default 'w'
        """
        super().__init__(filename, mode)

        self.format = MDEngineFormat(format)

    def write(self, frame: Frame) -> None:
        """
        Writes the frame to the file.

        Parameters
        ----------
        frame : Frame
            The frame to write.
        """
        self.open()
        self._write_box(frame.cell)
        self._write_atoms(frame)
        self.close()

    def _write_box(self, cell: Cell) -> None:
        """
        Writes the box to the file.

        Parameters
        ----------
        frame : Frame
            The frame to write.
        """
        print(
            f"Box  {cell.x} {cell.y} {cell.z}  {cell.alpha} {cell.beta} {cell.gamma}", file=self.file)

    def _write_atoms(self, frame: Frame) -> None:
        """
        Writes the atoms to the file.

        Parameters
        ----------
        frame : Frame
            The frame to write.
        """
        mol_types = frame.topology.moltype_ids

        for i in range(frame.n_atoms):
            atom = frame.system.atoms[i]
            pos = frame.system.pos[i]

            try:
                vel = frame.system.vel[i]
            except:
                vel = np.zeros(3)

            try:
                force = frame.system.forces[i]
            except:
                force = np.zeros(3)

            mol_type = mol_types[i]
            print(f"{atom.name}    {i}    {mol_type}",
                  file=self.file, end="    ")
            print(
                f"{pos[0]} {pos[1]} {pos[2]}", file=self.file, end=" ")
            print(
                f"{vel[0]} {vel[1]} {vel[2]}", file=self.file, end=" ")
            print(
                f"{force[0]} {force[1]} {force[2]}", file=self.file, end=" ")

            if self.format == MDEngineFormat.PIMD_QMCF:
                print(file=self.file)
            else:
                print(f"{pos[0]} {pos[1]} {pos[2]}", file=self.file, end=" ")
                print(f"{vel[0]} {vel[1]} {vel[2]}", file=self.file, end=" ")
                print(f"{force[0]} {force[1]} {force[2]}", file=self.file)
