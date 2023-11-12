"""
A module containing classes for reading a frame from a string.

...

Classes
-------
FrameReader
    A class for reading a frame from a string.
"""

from __future__ import annotations

import numpy as np

from beartype.typing import List, Tuple

from ..core.atomicSystem import AtomicSystem
from ..core.atom import Atom
from ..core.cell import Cell
from ..types import Numpy2DFloatArray, Numpy1DFloatArray
from ..traj.frame import Frame
from ..traj.formats import TrajectoryFormat
from ..exceptions import ElementNotFoundError


class FrameReader:
    """
    FrameReader reads a frame from a string.
    """

    def read(self, frame_string: str, format: TrajectoryFormat | str = TrajectoryFormat.XYZ) -> Frame:
        """
        Reads a frame from a string.

        Parameters
        ----------
        frame_string : str
            The string to read the frame from.
        format : TrajectoryFormat | str, optional
            The format of the trajectory. Default is TrajectoryFormat.XYZ.

        Returns
        -------
        Frame
            The frame read from the string.

        Raises
        ------
        ValueError
            If the given format is not valid.
        """

        # Note: TrajectoryFormat(format) automatically gives an error if format is not a valid TrajectoryFormat

        if TrajectoryFormat(format) is TrajectoryFormat.XYZ:
            return self.read_positions(frame_string)
        elif TrajectoryFormat(format) is TrajectoryFormat.VEL:
            return self.read_velocities(frame_string)
        elif TrajectoryFormat(format) is TrajectoryFormat.FORCE:
            return self.read_forces(frame_string)
        elif TrajectoryFormat(format) is TrajectoryFormat.CHARGE:
            return self.read_charges(frame_string)

    def read_positions(self, frame_string: str) -> Frame:
        """
        Reads the positions of the atoms in a frame from a string.

        Parameters
        ----------
        frame_string : str
            The string to read the frame from.

        Returns
        -------
        Frame
            The frame read from the string.

        Raises
        ------
        TypeError
            If the given frame_string is not a string.
        """

        splitted_frame_string = frame_string.split('\n')
        header_line = splitted_frame_string[0]

        n_atoms, cell = self._read_header_line(header_line)

        xyz, atoms = self._read_xyz(splitted_frame_string, n_atoms)

        try:
            atoms = [Atom(atom) for atom in atoms]
        except ElementNotFoundError:
            atoms = [Atom(atom, use_guess_element=False) for atom in atoms]

        return Frame(AtomicSystem(atoms=atoms, pos=xyz, cell=cell))

    def read_velocities(self, frame_string: str) -> Frame:
        """
        Reads the velocities of the atoms in a frame from a string.

        Parameters
        ----------
        frame_string : str
            The string to read the frame from.

        Returns
        -------
        Frame
            The frame read from the string.

        Raises
        ------
        TypeError
            If the given frame_string is not a string.
        """

        splitted_frame_string = frame_string.split('\n')
        header_line = splitted_frame_string[0]

        n_atoms, cell = self._read_header_line(header_line)

        vel, atoms = self._read_xyz(splitted_frame_string, n_atoms)

        try:
            atoms = [Atom(atom) for atom in atoms]
        except ElementNotFoundError:
            atoms = [Atom(atom, use_guess_element=False) for atom in atoms]

        return Frame(AtomicSystem(atoms=atoms, vel=vel, cell=cell))

    def read_forces(self, frame_string: str) -> Frame:
        """
        Reads the forces of the atoms in a frame from a string.

        Parameters
        ----------
        frame_string : str
            The string to read the frame from.

        Returns
        -------
        Frame
            The frame read from the string.

        Raises
        ------
        TypeError
            If the given frame_string is not a string.
        """

        splitted_frame_string = frame_string.split('\n')
        header_line = splitted_frame_string[0]

        n_atoms, cell = self._read_header_line(header_line)

        forces, atoms = self._read_xyz(splitted_frame_string, n_atoms)

        try:
            atoms = [Atom(atom) for atom in atoms]
        except ElementNotFoundError:
            atoms = [Atom(atom, use_guess_element=False) for atom in atoms]

        return Frame(AtomicSystem(atoms=atoms, forces=forces, cell=cell))

    def read_charges(self, frame_string: str) -> Frame:
        """
        Reads the charge values of the atoms in a frame from a string.

        Parameters
        ----------
        frame_string : str
            The string to read the frame from.

        Returns
        -------
        Frame
            The frame read from the string.

        Raises
        ------
        TypeError
            If the given frame_string is not a string.
        """

        splitted_frame_string = frame_string.split('\n')
        header_line = splitted_frame_string[0]

        n_atoms, cell = self._read_header_line(header_line)

        charges, atoms = self._read_scalar(splitted_frame_string, n_atoms)

        try:
            atoms = [Atom(atom) for atom in atoms]
        except ElementNotFoundError:
            atoms = [Atom(atom, use_guess_element=False) for atom in atoms]

        return Frame(AtomicSystem(atoms=atoms, charges=charges, cell=cell))

    def _read_header_line(self, header_line: str) -> Tuple[int, Cell | None]:
        """
        Reads the header line of a frame.

        It reads the number of atoms and the cell information from the header line.
        If the header line contains only the number of atoms, the cell is set to None.
        If the header line contains only the number of atoms and the box dimensions,
         the cell is set to a Cell object with the given box dimensions and box angles set to 90°.

        Parameters
        ----------
        header_line : str
            The header line to read.

        Returns
        -------
        n_atoms : int
            The number of atoms in the frame.
        cell : Cell
            The cell of the frame.

        Raises
        ------
        ValueError
            If the header line is not valid. Either it contains too many or too few values.
        """

        header_line = header_line.split()

        if len(header_line) == 4:
            n_atoms = int(header_line[0])
            a, b, c = map(float, header_line[1:4])
            cell = Cell(a, b, c)
        elif len(header_line) == 7:
            n_atoms = int(header_line[0])
            a, b, c, alpha, beta, gamma = map(float, header_line[1:7])
            cell = Cell(a, b, c, alpha, beta, gamma)
        elif len(header_line) == 1:
            n_atoms = int(header_line[0])
            cell = None
        else:
            raise ValueError('Invalid file format in header line of Frame.')

        return n_atoms, cell

    def _read_xyz(self, splitted_frame_string: List[str], n_atoms: int) -> Tuple[Numpy2DFloatArray, List[str]]:
        """
        Reads the xyz coordinates and the atom names from the given string.

        Parameters
        ----------
        splitted_frame_string : str
            The string to read the xyz coordinates and the atom names from.
        n_atoms : int
            The number of atoms in the frame.

        Returns
        -------
        xyz : np.array
            The xyz coordinates of the atoms.
        atoms : list of str
            The names of the atoms.

        Raises
        ------
        ValueError
            If the given string does not contain the correct number of lines.
        """

        xyz = np.zeros((n_atoms, 3))
        atoms = []
        for i in range(n_atoms):
            line = splitted_frame_string[2+i]

            if len(line.split()) != 4:
                raise ValueError(
                    'Invalid file format in xyz coordinates of Frame.')

            xyz[i] = np.array([float(x) for x in line.split()[1:4]])
            atoms.append(line.split()[0])

        return xyz, atoms

    def _read_scalar(self, splitted_frame_string: List[str], n_atoms: int) -> Tuple[Numpy1DFloatArray, List[str]]:
        """
        Reads the scalar values and the atom names from the given string.

        Parameters
        ----------
        splitted_frame_string : str
            The string to read the scalar values and the atom names from.
        n_atoms : int
            The number of atoms in the frame.

        Returns
        -------
        scalar : np.array
            The scalar values of the atoms.
        atoms : list of str
            The names of the atoms.

        Raises
        ------
        ValueError
            If the given string does not contain the correct number of lines.
        """

        scalar = np.zeros((n_atoms))
        atoms = []
        for i in range(n_atoms):
            line = splitted_frame_string[2+i]

            if len(line.split()) != 2:
                raise ValueError(
                    'Invalid file format in scalar values of Frame.')

            scalar[i] = float(line.split()[1])
            atoms.append(line.split()[0])

        return scalar, atoms
