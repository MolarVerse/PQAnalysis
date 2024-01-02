"""
A module containing the RestartFileReader class

Classes
-------
RestartFileReader
    A class for reading restart files.
"""

import numpy as np

from beartype.typing import List

from . import BaseReader, RestartFileReaderError
from ..core import AtomicSystem, Atom, Cell
from ..traj import MDEngineFormat, Frame
from ..topology import Topology


class RestartFileReader(BaseReader):
    """
    A class for reading restart files.

    Inherits from the BaseReader class.
    It reads restart files from the following molecular dynamics engines:
        - PIMD_QMCF
        - QMCFC


    Attributes
    ----------
    filename : str
        The filename of the restart file.
    format : MDEngineFormat
        The format of the restart file.
    """

    def __init__(self, filename: str, format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF) -> None:
        """
        Initializes the RestartFileReader with the given parameters.

        It automatically checks if the file exists.

        Parameters
        ----------
        filename : str
            The filename of the restart file.
        format : MDEngineFormat | str, optional
            The format of the restart file, by default MDEngineFormat.PIMD_QMCF
        """
        super().__init__(filename)

        self.format = MDEngineFormat(format)

    def read(self) -> Frame:
        """
        Reads the restart file and returns an AtomicSystem and a Np1DIntArray containing the molecular types.

        Returns
        -------
        Frame:
            The Frame object including the AtomicSystem and the Topology with the molecular types.
        """

        cell = Cell()
        atom_lines = []

        with open(self.filename, 'r') as file:
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

        return self._parse_atoms(atom_lines, cell)

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
        elif len(line) == 4:
            box_lengths = [float(l) for l in line[1:]]
            return Cell(*box_lengths)
        elif len(line) == 7:
            box_lengths = [float(l) for l in line[1:4]]
            box_angles = [float(a) for a in line[4:]]
            return Cell(*box_lengths, *box_angles)
        else:
            raise RestartFileReaderError(
                f"Invalid number of arguments for box: {len(line)}")

    @classmethod
    def _parse_atoms(cls, lines: List[str], cell: Cell = Cell()) -> Frame:
        """
        Parses the atom lines of the restart file.

        An atom line can have 12 or 21 arguments.
        if it has 12 arguments, the atom line is assumed to be from a PIMD_QMCF restart file.
        if it has 21 arguments, the atom line is assumed to be from a QMCFC restart file.

        For the PIMD_QMCF restart file, the arguments are:
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
        Frame:
            The Frame object including the AtomicSystem and the Topology with the molecular types.

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
                raise RestartFileReaderError(
                    f"Invalid number of arguments for atom: {len(line)}")

            atoms.append(Atom(line[0], use_guess_element=False))
            residues.append(int(line[2]))
            positions.append(np.array([float(l) for l in line[3:6]]))
            velocities.append(np.array([float(l) for l in line[6:9]]))
            forces.append(np.array([float(l) for l in line[9:12]]))

        if atoms == []:
            raise RestartFileReaderError("No atoms found in restart file.")

        topology = Topology(atoms=atoms, residue_ids=np.array(residues))

        system = AtomicSystem(pos=np.array(positions), vel=np.array(
            velocities), forces=np.array(forces), cell=cell, topology=topology)

        return Frame(system=system)
