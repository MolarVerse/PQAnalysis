"""
A module for writing restart files.
"""

import numpy as np

from . import BaseWriter
from PQAnalysis.traj import MDEngineFormat, Frame
from PQAnalysis.core import Cell


class RestartFileWriter(BaseWriter):
    """
    A class for writing restart files.

    Inherits from the BaseWriter class.
    It writes restart files from the following molecular dynamics engines:
    - PIMD_QMCF
    - QMCFC

    The restart file of these two md-engines are very similar. Both contain one line including the step information and one line including the box information. The following lines contain the atom information. The atom information is different for the two md-engines. The atom information of the PIMD_QMCF restart file contains the atom name, atom id, residue id, x position, y position, z position, x velocity, y velocity, z velocity, x force, y force and z force. The atom information of the QMCFC restart file contains the atom name, atom id, residue id, x position, y position, z position, x velocity, y velocity, z velocity, x force, y force, z force, x pos of previous step, y pos of previous step, z pos of previous step, x vel of previous step, y vel of previous step, z vel of previous step, x force of previous step, y force of previous step and z force of previous step. The old values are ignored.

    For more information on how the restart file of PIMD_QMCF simulations is structured, see the corresponding documentation of the `PIMD-QMCF <https://molarverse.github.io/pimd_qmcf>`_ code.
    """

    def __init__(self,
                 filename: str | None = None,
                 md_engine_format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF,
                 mode: str = 'w'
                 ) -> None:
        """
        Parameters
        ----------
        filename : str | None, optional
            The filename of the restart file, by default None (stdout)
        md_engine_format : MDEngineFormat | str, optional
            The format of the restart file, by default MDEngineFormat.PIMD_QMCF
        mode : str, optional
            The mode of the file, by default 'w'
        """
        super().__init__(filename, mode)

        self.md_engine_format = MDEngineFormat(md_engine_format)

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
        residues = frame.topology.residue_ids

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

            residue = residues[i]
            print(f"{atom.name}    {i}    {residue}",
                  file=self.file, end="    ")
            print(
                f"{pos[0]} {pos[1]} {pos[2]}", file=self.file, end=" ")
            print(
                f"{vel[0]} {vel[1]} {vel[2]}", file=self.file, end=" ")
            print(
                f"{force[0]} {force[1]} {force[2]}", file=self.file, end=" ")

            if self.md_engine_format == MDEngineFormat.PIMD_QMCF:
                print(file=self.file)
            else:
                print(f"{pos[0]} {pos[1]} {pos[2]}", file=self.file, end=" ")
                print(f"{vel[0]} {vel[1]} {vel[2]}", file=self.file, end=" ")
                print(f"{force[0]} {force[1]} {force[2]}", file=self.file)
