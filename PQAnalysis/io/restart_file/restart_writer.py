"""
A module for writing restart files.
"""
import logging

import numpy as np

from beartype.typing import List

from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.core import Cell
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.io.base import BaseWriter
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

from .exceptions import RestartFileWriterError



class RestartFileWriter(BaseWriter):

    """
    A class for writing restart files.

    Inherits from the BaseWriter class.
    It writes restart files from the following molecular dynamics engines:
    - PQ
    - QMCFC

    The restart file of these two md-engines are very similar.
    Both contain one line including the step information and one
    line including the box information. The following lines 
    contain the atom information. The atom information is 
    different for the two md-engines. The atom information of
    the PQ restart file contains the atom name, atom id, 
    residue id, x position, y position, z position, x velocity,
    y velocity, z velocity, x force, y force and z force. 
    The atom information of the QMCFC restart file contains
    the atom name, atom id, residue id, x position, y position,
    z position, x velocity, y velocity, z velocity, x force, 
    y force, z force, x pos of previous step, y pos of previous
    step, z pos of previous step, x vel of previous step, 
    y vel of previous step, z vel of previous step, 
    x force of previous step, y force of previous step and
    z force of previous step. The old values are ignored.

    For more information on how the restart file of PQ simulations
    is structured, see the corresponding documentation of the
    `PQ <https://molarverse.github.io/PQ>`_ code.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        filename: str | None = None,
        md_engine_format: MDEngineFormat | str = MDEngineFormat.PQ,
        mode: FileWritingMode | str = 'w'
    ) -> None:
        """
        Parameters
        ----------
        filename : str | None, optional
            The filename of the restart file, by default None (stdout)
        md_engine_format : MDEngineFormat | str, optional
            The format of the restart file, by default MDEngineFormat.PQ
        mode : FileWritingMode | str, optional
            The writing mode, by default 'w' - 
            Possible values are 'w' (write), 'a' (append) and 'o' (overwrite).
        """
        super().__init__(filename, mode)

        self.md_engine_format = MDEngineFormat(md_engine_format)

    @runtime_type_checking
    def write(
        self,
        frame: AtomicSystem,
        atom_counter: int | Np1DNumberArray | None = None,
    ) -> None:
        """
        Writes the frame to the file.

        Parameters
        ----------
        frame : AtomicSystem
            The frame to write.
        atom_counter : int | Np1DNumberArray | None, optional
            The atom counter, by default None. If only a single
            integer is given, this number will be used as the atom 
            counter for all atoms. If an array is given, this array
            has to have the same length as the number of atoms 
            in the frame.
        """

        lines = self._get_lines(frame, atom_counter)

        self._write_lines_to_file(lines)

    def _write_lines_to_file(self, lines: List[str]) -> None:
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

    def _get_lines(
        self,
        frame: AtomicSystem,
        atom_counter: int | Np1DNumberArray | None = None,
    ) -> List[str]:
        """
        Collects the lines to write to the file.

        Parameters
        ----------
        frame : AtomicSystem
            The frame to write.
        atom_counter : int | Np1DNumberArray | None, optional
            The atom counter, by default None. If only a single integer is
            given, this number will be used as the atom counter for all 
            atoms. If an array is given, this array has to have the same 
            length as the number of atoms in the frame.
        """

        lines = []
        lines.append(self._get_box_line(frame.cell))

        lines += self._get_atom_lines(
            frame,
            atom_counter,
            self.md_engine_format
        )

        return lines

    def _get_box_line(self, cell: Cell) -> str:
        """
        Writes the box to the file.

        Parameters
        ----------
        frame : Frame
            The frame to write.
        """
        return (
            f"Box  {cell.x} {cell.y} {cell.z}  "
            f"{cell.alpha} {cell.beta} {cell.gamma}"
        )

    @classmethod
    def _get_atom_lines(
        cls,
        frame: AtomicSystem,
        atom_counter: int | Np1DNumberArray | None = None,
        md_engine_format: MDEngineFormat | str = MDEngineFormat.PQ
    ) -> List[str]:
        """
        Writes the atoms to the file.

        Parameters
        ----------
        frame : AtomicSystem
            The frame to write.
        atom_counter : int | Np1DNumberArray | None, optional
            The atom counter, by default None. If only a single integer is 
            given, this number will be used as the atom counter for all 
            atoms. If an array is given, this array has to have the same 
            length as the number of atoms in the frame.
        md_engine_format : MDEngineFormat | str, optional
            The format of the restart file, by default MDEngineFormat.PQ
        """

        if atom_counter is not None and not isinstance(atom_counter, int):
            if len(atom_counter) != frame.n_atoms:
                cls.logger.error(
                    (
                    "The atom counter has to have the same length as "
                    "the number of atoms in the frame if it is given "
                    "as an array."
                    ),
                    exception=RestartFileWriterError
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

            vel = frame.vel[i] if frame.has_vel else np.zeros(3)
            force = frame.forces[i] if frame.has_forces else np.zeros(3)

            residue = residues[i]

            line = ""

            line += f"{atom.name}    {atom_counter[i]}    {residue}    "
            line += f"{pos[0]} {pos[1]} {pos[2]} "
            line += f"{vel[0]} {vel[1]} {vel[2]} "
            line += f"{force[0]} {force[1]} {force[2]}"

            if md_engine_format != MDEngineFormat.PQ:
                line += f" {pos[0]} {pos[1]} {pos[2]} "
                line += f"{vel[0]} {vel[1]} {vel[2]} "
                line += f"{force[0]} {force[1]} {force[2]}"

            lines.append(line)

        return lines
