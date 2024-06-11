"""
A module containing classes for reading a frame from a string.
"""

import logging

import numpy as np

from beartype.typing import List, Tuple

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell, ElementNotFoundError
from PQAnalysis.types import Np2DNumberArray, Np1DNumberArray
from PQAnalysis.traj import TrajectoryFormat, MDEngineFormat
from PQAnalysis.topology import Topology
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__

from .exceptions import FrameReaderError
from .process_lines import process_lines_with_atoms  # pylint: disable=import-error



class _FrameReader:

    """
    This class provides methods for reading a frame from a string.
    The string can be a single frame or a whole trajectory.
    The format of the string can be specified with the format parameter.
    Generally, this class should not be used directly.
    Instead, the :py:class:`~PQAnalysis.io.trajectoryReader.TrajectoryReader`
    class can be used to read a whole trajectory as well as
    a single frame from a file.

    For more information about the format of the string,
    see :py:class:`~PQAnalysis.traj.formats.TrajectoryFormat`.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    def __init__(
        self, md_format: MDEngineFormat | str = MDEngineFormat.PQ
    ) -> None:
        """
        Parameters
        ----------
        md_format : MDEngineFormat | str, optional
            The format of the MD engine. Default is "PQ".
        """
        self.md_format = MDEngineFormat(md_format)
        self.topology = None

    def read(
        self,
        frame_string: str,
        topology: Topology | None = None,
        traj_format: TrajectoryFormat | str = TrajectoryFormat.XYZ
    ) -> AtomicSystem:
        """
        Reads a frame from a string.

        Parameters
        ----------
        frame_string : str
            The string to read the frame from.
        topology : Topology | None, optional
            The topology of the frame. Default is None.
        traj_format : TrajectoryFormat | str, optional
            The format of the trajectory. Default is TrajectoryFormat.XYZ.

        Returns
        -------
        AtomicSystem
            The frame read from the string.

        Raises
        ------
        FrameReaderError
            If the given format is not valid.
        """

        # Note: TrajectoryFormat(format) automatically gives
        # an error if format is not a valid TrajectoryFormat

        self.topology = topology

        if TrajectoryFormat(traj_format) is TrajectoryFormat.XYZ:
            return self.read_positions(frame_string)

        if TrajectoryFormat(traj_format) is TrajectoryFormat.VEL:
            return self.read_velocities(frame_string)

        if TrajectoryFormat(traj_format) is TrajectoryFormat.FORCE:
            return self.read_forces(frame_string)

        if TrajectoryFormat(traj_format) is TrajectoryFormat.CHARGE:
            return self.read_charges(frame_string)

        # This should never happen - only for safety
        self.logger.error(
            f'Invalid TrajectoryFormat given. {traj_format=}',
            exception=FrameReaderError
        )

    def read_positions(self, frame_string: str) -> AtomicSystem:
        """
        Reads the positions of the atoms in a frame from a string.

        Parameters
        ----------
        frame_string : str
            The string to read the frame from.

        Returns
        -------
        AtomicSystem
            The frame read from the string.
        """

        splitted_frame_string = frame_string.split('\n')
        header_line = splitted_frame_string[0]

        n_atoms, cell = self._read_header_line(header_line)

        xyz, atoms = self._read_xyz(splitted_frame_string, n_atoms)
        xyz, atoms = self._check_qmcfc(atoms, xyz)

        topology = self._get_topology(atoms, self.topology)

        return AtomicSystem(topology=topology, pos=xyz, cell=cell)

    def read_velocities(self, frame_string: str) -> AtomicSystem:
        """
        Reads the velocities of the atoms in a frame from a string.

        Parameters
        ----------
        frame_string : str
            The string to read the frame from.

        Returns
        -------
        AtomicSystem
            The frame read from the string.
        """

        splitted_frame_string = frame_string.split('\n')
        header_line = splitted_frame_string[0]

        n_atoms, cell = self._read_header_line(header_line)

        vel, atoms = self._read_xyz(splitted_frame_string, n_atoms)
        vel, atoms = self._check_qmcfc(atoms, vel)

        topology = self._get_topology(atoms, self.topology)

        return AtomicSystem(topology=topology, vel=vel, cell=cell)

    def read_forces(self, frame_string: str) -> AtomicSystem:
        """
        Reads the forces of the atoms in a frame from a string.

        Parameters
        ----------
        frame_string : str
            The string to read the frame from.

        Returns
        -------
        AtomicSystem
            The frame read from the string.
        """

        splitted_frame_string = frame_string.split('\n')
        header_line = splitted_frame_string[0]

        n_atoms, cell = self._read_header_line(header_line)

        forces, atoms = self._read_xyz(splitted_frame_string, n_atoms)
        forces, atoms = self._check_qmcfc(atoms, forces)

        topology = self._get_topology(atoms, self.topology)

        return AtomicSystem(topology=topology, forces=forces, cell=cell)

    def read_charges(self, frame_string: str) -> AtomicSystem:
        """
        Reads the charge values of the atoms in a frame from a string.

        Parameters
        ----------
        frame_string : str
            The string to read the frame from.

        Returns
        -------
        AtomicSystem
            The frame read from the string.
        """

        splitted_frame_string = frame_string.split('\n')
        header_line = splitted_frame_string[0]

        n_atoms, cell = self._read_header_line(header_line)

        charges, atoms = self._read_scalar(splitted_frame_string, n_atoms)
        charges, atoms = self._check_qmcfc(atoms, charges)

        topology = self._get_topology(atoms, self.topology)

        return AtomicSystem(topology=topology, charges=charges, cell=cell)

    def _check_qmcfc(
        self,
        atoms: List[str],
        value: Np1DNumberArray | Np2DNumberArray,
    ) -> Tuple[Np1DNumberArray | Np2DNumberArray, List[str]]:
        """
        Check if the first atom is X for QMCFC. If it is, remove it from the list and array.

        Parameters
        ----------
        atoms : List[str]
            The list of atoms.
        value : Np1DNumberArray | Np2DNumberArray
            The array of values.

        Returns
        -------
        Tuple[List[str], Np1DNumberArray | Np2DNumberArray]
            The list of atoms and the array of values.

        Raises
        ------
        FrameReaderError
            If the first atom is not X for QMCFC.
        """

        if self.md_format == MDEngineFormat.QMCFC:
            if atoms[0].upper() != 'X':
                self.logger.error(
                    (
                        'The first atom in one of the frames is not X. '
                        'Please use PQ (default) md engine instead'
                    ),
                    exception=FrameReaderError
                )
            value = value[1:]
            atoms = atoms[1:]

        return value, atoms

    def _get_topology(
        self,
        atoms: List[str],
        topology: Topology | None,
    ) -> Topology:
        """
        Returns the topology of the frame.

        Returns
        -------
        Topology
            The topology of the frame.
        """

        if topology is None:
            try:

                topology = Topology(
                    atoms=[
                        Atom(atom, disable_type_checking=True)
                        for atom in atoms
                    ]
                )

            except ElementNotFoundError:

                topology = Topology(
                    atoms=[
                        Atom(
                            atom,
                            use_guess_element=False,
                            disable_type_checking=True
                        ) for atom in atoms
                    ]
                )

        return topology

    def _read_header_line(self, header_line: str) -> Tuple[int, Cell]:
        """
        Reads the header line of a frame.

        It reads the number of atoms and the cell information from the header line.
        If the header line contains only the number of atoms, the cell is set Cell().
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
        FrameReaderError
            If the header line is not valid. Either it contains too many or too few values.
        """

        n_atoms = 0  # default value
        cell = None  # default value

        header_line = header_line.split()

        if len(header_line) == 4:
            n_atoms = int(header_line[0])
            a, b, c = (float(x) for x in header_line[1:4])
            cell = Cell(a, b, c)
        elif len(header_line) == 7:
            n_atoms = int(header_line[0])
            a, b, c, alpha, beta, gamma = (float(x) for x in header_line[1:7])
            cell = Cell(a, b, c, alpha, beta, gamma)
        elif len(header_line) == 1:
            n_atoms = int(header_line[0])
            cell = Cell()
        else:
            self.logger.error(
                'Invalid file format in header line of Frame.',
                exception=FrameReaderError
            )

        return n_atoms, cell

    def _read_xyz(
        self,
        splitted_frame_string: List[str],
        n_atoms: int,
    ) -> Tuple[Np2DNumberArray, List[str]]:
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
        FrameReaderError
            If the given string does not contain the correct number of lines.
        """
        try:
            atoms, xyz = process_lines_with_atoms(splitted_frame_string[2:], n_atoms)

            return xyz, atoms
        except ValueError:
            self.logger.error(
                'Invalid file format in xyz coordinates of Frame.',
                exception=FrameReaderError
            )

    def _read_scalar(
        self,
        splitted_frame_string: List[str],
        n_atoms: int,
    ) -> Tuple[Np1DNumberArray, List[str]]:
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
        FrameReaderError
            If the given string does not contain the correct number of lines.
        """

        scalar = np.zeros((n_atoms))
        atoms = []
        for i in range(n_atoms):
            line = splitted_frame_string[2 + i]

            if len(line.split()) != 2:
                self.logger.error(
                    'Invalid file format in scalar values of Frame.',
                    exception=FrameReaderError
                )

            scalar[i] = float(line.split()[1])
            atoms.append(line.split()[0])

        return scalar, atoms
