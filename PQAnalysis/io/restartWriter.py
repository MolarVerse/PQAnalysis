import numpy as np

from .base import BaseWriter
from ..traj.formats import MDEngineFormat
from ..traj.frame import Frame
from ..core.cell import Cell


class RestartFileWriter(BaseWriter):
    def __init__(self,
                 filename: str | None = None,
                 format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF,
                 mode: str = 'w'
                 ) -> None:
        super().__init__(filename, mode)

        self.format = MDEngineFormat(format)

        if self.format == MDEngineFormat.AUTO:
            raise ValueError(
                "The format AUTO is not supported for writing restart files.")

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
        if frame.topology == None or frame.topology.mol_types is None:
            mol_types = np.zeros(frame.system.n_atoms, dtype=int)
        else:
            if frame.n_atoms != frame.topology.mol_types:
                raise ValueError(
                    "The number of mol_types does not match the number of atoms.")
            mol_types = frame.topology.mol_types

        for i in range(frame.n_atoms):
            atom = frame.system.atoms[i]
            pos = frame.system.pos[i]
            vel = frame.system.vel[i]
            force = frame.system.forces[i]
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
