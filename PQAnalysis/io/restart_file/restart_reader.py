"""
A module containing the RestartFileReader class
"""

import logging

import numpy as np

from beartype.typing import List

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell, Residues
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.topology import Topology
from PQAnalysis.io.base import BaseReader
from PQAnalysis.io.moldescriptor_reader import MoldescriptorReader
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

from .exceptions import RestartFileReaderError



class RestartFileReader(BaseReader):

    """
    A class for reading restart files.

    Inherits from the BaseReader class.
    It reads restart files from the following molecular dynamics engines:
    - PQ
    - QMCFC

    The restart file of these two md-engines are very similar. Both contain 
    one line including the step information and one line including the box
    information. The following lines contain the atom information. The atom
    information is different for the two md-engines. The atom information 
    of the PQ restart file contains the atom name, atom id, residue id, 
    x position, y position, z position, x velocity, y velocity, z velocity,
    x force, y force and z force. The atom information of the QMCFC restart
    file contains the atom name, atom id, residue id, x position, y position,
    z position, x velocity, y velocity, z velocity, x force, y force,
    z force, x pos of previous step, y pos of previous step, z pos of 
    previous step, x vel of previous step, y vel of previous step, 
    z vel of previous step, x force of previous step, y force of previous 
    step and z force of previous step. The old values are ignored.

    For more information on how the restart file of PQ simulations is 
    structured, see the corresponding documentation of the 
    `PQ <https://molarverse.github.io/PQ>`_ code.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        filename: str,
        moldescriptor_filename: str | None = None,
        reference_residues: Residues | None = None,
        md_engine_format: MDEngineFormat | str = MDEngineFormat.PQ
    ) -> None:
        """
        Parameters
        ----------
        filename : str
            The filename of the restart file.
        moldescriptor_filename : str, optional
            The filename of the moldescriptor file that is read by the
            MoldescriptorReader to obtain the reference residues of 
            the system, by default None
        reference_residues : Residues, optional
            The reference residues of the system, in general these are
            obtained by the MoldescriptorReader - only used if 
            moldescriptor_filename is None, by default None
        md_engine_format : MDEngineFormat | str, optional
            The format of the restart file, by default MDEngineFormat.PQ
        """
        super().__init__(filename)

        if moldescriptor_filename is not None and reference_residues is not None:
            self.logger.error(
                (
                    "Both moldescriptor_filename and reference_residues "
                    "are given. They are mutually exclusive."
                ),
                exception=RestartFileReaderError
            )

        self.moldescriptor_filename = moldescriptor_filename
        self.reference_residues = reference_residues

        self.md_engine_format = MDEngineFormat(md_engine_format)

    def read(self) -> AtomicSystem:
        """
        Reads the restart file and returns an AtomicSystem and
        an Np1DIntArray containing the molecular types.

        It reads the restart file and extracts the box information
        and the atom information. The atom information is then used
        to create an AtomicSystem object. The box information is used
        to create a Cell object. The AtomicSystem and the Cell are
        then used to create a Frame object.

        If a moldescriptor file is given, the molecular types are 
        read from the moldescriptor file and added to the Topology 
        of the Frame.

        Returns
        -------
        AtomicSystem:
            The AtomicSystem object including the Topology with the molecular types.
        """

        if self.moldescriptor_filename is not None:
            self.reference_residues = MoldescriptorReader(
                self.moldescriptor_filename
            ).read()

        cell = Cell()
        atom_lines = []

        with open(self.filename, 'r', encoding="utf-8") as file:
            for line in file.readlines():
                if line.strip().startswith('#'):
                    continue

                if len(line.strip().split()) == 0:
                    continue

                line = line.strip().split()

                if line[0].lower() == "box":
                    cell = self._parse_box(line)
                elif line[0].lower() == "step":
                    continue
                elif line[0].lower() == "chi":
                    continue
                else:
                    atom_lines.append(" ".join(line))

        return self._parse_atoms(atom_lines, cell, self.reference_residues)

    @classmethod
    def _parse_box(cls, line: List[str]) -> Cell:
        """
        Parses the box line of the restart file.

        The box line can have 1, 4 or 7 arguments.
        If it has 1 argument, the box is assumed to be a vacuum box.
        If it has 4 arguments, the box is assumed to be a orthorhombic box.
        If it has 7 arguments, the box is assumed to be a triclinic box.

        Parameters
        ----------
        line : List[str]
            The box line of the restart file.

        Returns
        -------
        Cell
            The Cell object.

        Raises
        ------
        RestartFileReaderError
            If the number of arguments is not 1, 4 or 7.
        """
        if len(line) == 1:
            return Cell()

        if len(line) == 4:
            box_lengths = [float(l) for l in line[1:]]
            return Cell(*box_lengths)

        if len(line) == 7:
            box_lengths = [float(l) for l in line[1:4]]
            box_angles = [float(a) for a in line[4:]]
            return Cell(*box_lengths, *box_angles)

        cls.logger.error(
            f"Invalid number of arguments for box: {len(line)}",
            exception=RestartFileReaderError
        )

    @classmethod
    def _parse_atoms(
        cls,
        lines: List[str],
        cell: Cell = Cell(),
        reference_residues: Residues | None = None
    ) -> AtomicSystem:
        """
        Parses the atom lines of the restart file.

        An atom line can have 12 or 21 arguments.
        if it has 12 arguments, the atom line is assumed to be from a PQ restart file.
        if it has 21 arguments, the atom line is assumed to be from a QMCFC restart file.

        For the PQ restart file, the arguments are:
            - atom name
            - atom id
            - residue id
            - x position
            - y position
            - z position
            - x velocity
            - y velocity
            - z velocity
            - x force
            - y force
            - z force

        For the QMCFC restart file, the arguments are:
            - atom name
            - atom id
            - residue id
            - x position
            - y position
            - z position
            - x velocity
            - y velocity
            - z velocity
            - x force
            - y force
            - z force
            - x pos of previous step
            - y pos of previous step
            - z pos of previous step
            - x vel of previous step
            - y vel of previous step
            - z vel of previous step
            - x force of previous step
            - y force of previous step
            - z force of previous step

            where all old values are ignored.

        Parameters
        ----------
        lines : List[str]
            The atom lines of the restart file.
        cell : Cell, optional
            The cell of the AtomicSystem, by default Cell()

        Returns
        -------
        AtomicSystem:
            The Frame object including the Topology with the molecular types.

        Raises
        ------
        RestartFileReaderError
            If the number of arguments is not 12 or 21.
        RestartFileReaderError
            If no atoms are found in the restart file.
        """
        atoms = []
        positions = []
        velocities = []
        forces = []
        residues = []

        for line in lines:
            line = line.strip().split()

            if len(line) != 12 and len(line) != 21:
                cls.logger.error(
                    f"Invalid number of arguments for atom: {len(line)}",
                    exception=RestartFileReaderError
                )

            atoms.append(Atom(line[0]))
            residues.append(int(line[2]))
            positions.append(np.array([float(l) for l in line[3:6]]))
            velocities.append(np.array([float(l) for l in line[6:9]]))
            forces.append(np.array([float(l) for l in line[9:12]]))

        if not atoms:
            cls.logger.error(
                "No atoms found in restart file.",
                exception=RestartFileReaderError
            )

        topology = Topology(
            atoms=atoms,
            residue_ids=np.array(residues),
            reference_residues=reference_residues
        )

        return AtomicSystem(
            pos=np.array(positions),
            vel=np.array(velocities),
            forces=np.array(forces),
            cell=cell,
            topology=topology
        )
