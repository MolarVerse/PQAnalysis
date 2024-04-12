"""
A module for writing restart files.
"""

import numpy as np

from beartype.typing import List

from .. import BaseWriter, FileWritingMode
from PQAnalysis.traj import MDEngineFormat, Frame
from PQAnalysis.core import Cell
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.atomicSystem import AtomicSystem


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
                 mode: FileWritingMode | str = 'w'
                 ) -> None:
        """
        Parameters
        ----------
        filename : str | None, optional
            The filename of the restart file, by default None (stdout)
        md_engine_format : MDEngineFormat | str, optional
            The format of the restart file, by default MDEngineFormat.PIMD_QMCF
        mode : FileWritingMode | str, optional
            The writing mode, by default 'w' - Possible values are 'w' (write), 'a' (append) and 'o' (overwrite).
        """
        super().__init__(filename, mode)

        self.md_engine_format = MDEngineFormat(md_engine_format)

    def write(self,
              frame: Frame | AtomicSystem,
              atom_counter: int | Np1DNumberArray | None = None,
              ) -> None:
        """
        Writes the frame to the file.

        Parameters
        ----------
        frame : Frame | AtomicSystem
            The frame to write.
        atom_counter : int | Np1DNumberArray | None, optional
            The atom counter, by default None. If only a single integer is given, this number will be used as the atom counter for all atoms. If an array is given, this array has to have the same length as the number of atoms in the frame.
        """

        lines = self.get_lines(frame, atom_counter)

        self.write_lines_to_file(lines)

    def write_lines_to_file(self, lines: List[str]) -> None:
        """
        Writes the lines to the file.

        Parameters
        ----------
        lines : List[str]
            The lines to write.
        """

        self.open()

        for line in lines:
            print(line, file=self.file)

        self.close()

    def get_lines(self,
                  frame: Frame | AtomicSystem,
                  atom_counter: int | Np1DNumberArray | None = None,
                  ) -> List[str]:
        """
        Collects the lines to write to the file.

        Parameters
        ----------
        frame : Frame | AtomicSystem
            The frame to write.
        atom_counter : int | Np1DNumberArray | None, optional
            The atom counter, by default None. If only a single integer is given, this number will be used as the atom counter for all atoms. If an array is given, this array has to have the same length as the number of atoms in the frame.
        """

        lines = []
        lines.append(self._write_box(frame.cell))
        lines += self._write_atoms(frame, atom_counter)

        return lines

    def _write_box(self, cell: Cell) -> str:
        """
        Writes the box to the file.

        Parameters
        ----------
        frame : Frame
            The frame to write.
        """
        return f"Box  {cell.x} {cell.y} {cell.z}  {cell.alpha} {cell.beta} {cell.gamma}"

    def _write_atoms(self,
                     frame: Frame | AtomicSystem,
                     atom_counter: int | Np1DNumberArray | None = None
                     ) -> List[str]:
        """
        Writes the atoms to the file.

        Parameters
        ----------
        frame : Frame | AtomicSystem
            The frame to write.
        atom_counter : int | Np1DNumberArray | None, optional
            The atom counter, by default None. If only a single integer is given, this number will be used as the atom counter for all atoms. If an array is given, this array has to have the same length as the number of atoms in the frame.
        """

        if atom_counter is not None and not isinstance(atom_counter, int):
            if len(atom_counter) != frame.n_atoms:
                raise ValueError(
                    "The atom counter has to have the same length as the number of atoms in the frame."
                )
        elif isinstance(atom_counter, int):
            atom_counter = [atom_counter] * frame.n_atoms
        elif atom_counter is None:
            atom_counter = range(frame.n_atoms)

        residues = frame.topology.residue_ids

        lines = []

        for i in range(frame.n_atoms):
            atom = frame.atoms[i]
            pos = frame.pos[i]

            try:
                vel = frame.vel[i]
            except Exception:
                vel = np.zeros(3)

            try:
                force = frame.forces[i]
            except Exception:
                force = np.zeros(3)

            residue = residues[i]

            line = ""

            line += f"{atom.name}    {atom_counter[i]}    {residue}    "
            line += f"{pos[0]} {pos[1]} {pos[2]} "
            line += f"{vel[0]} {vel[1]} {vel[2]} "
            line += f"{force[0]} {force[1]} {force[2]}"

            if self.md_engine_format != MDEngineFormat.PIMD_QMCF:
                line += f" {pos[0]} {pos[1]} {pos[2]} "
                line += f"{vel[0]} {vel[1]} {vel[2]} "
                line += f"{force[0]} {force[1]} {force[2]}"

            lines.append(line)

        return lines
