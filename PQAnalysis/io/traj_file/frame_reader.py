"""
A module containing classes for reading a frame from a string.
"""

import logging
import shlex
from abc import ABC, abstractmethod

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

try:
    from .process_lines import process_lines_with_atoms  # pylint: disable=import-error
except ModuleNotFoundError:
    from ._process_lines_py import process_lines_with_atoms



class BaseFrameReader(ABC):

    """
    Base class for frame readers.

    Subclasses implement the parsing for a concrete frame format while this
    class provides shared MD-engine handling and topology construction.
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

    @abstractmethod
    def read(
        self,
        frame_string: str,
        topology: Topology | None = None,
        traj_format: TrajectoryFormat | str = TrajectoryFormat.XYZ
    ) -> AtomicSystem:
        """
        Reads a frame from a string.
        """

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


class XYZFrameReader(BaseFrameReader):

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
        message = f'Invalid TrajectoryFormat given. {traj_format=}'
        self.logger.error(message, exception=FrameReaderError)
        raise FrameReaderError(message)

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


class _FrameReader(XYZFrameReader):

    """
    Backward-compatible alias for the default xyz-family frame reader.
    """


class ExtXYZFrameReader(BaseFrameReader):

    """
    A frame reader for extended xyz position frames.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    def read(
        self,
        frame_string: str,
        topology: Topology | None = None,
        traj_format: TrajectoryFormat | str = TrajectoryFormat.EXTXYZ
    ) -> AtomicSystem:
        """
        Reads an extended xyz frame from a string.
        """
        traj_format = TrajectoryFormat(traj_format)
        if traj_format is not TrajectoryFormat.EXTXYZ:
            message = f'Invalid TrajectoryFormat given. {traj_format=}'
            self.logger.error(message, exception=FrameReaderError)

        self.topology = topology

        lines = frame_string.split('\n')
        n_atoms = self._read_atom_count(lines[0])
        metadata = self._read_metadata(lines[1])
        properties = self._read_properties(metadata)

        values, atoms = self._read_property_values(lines, n_atoms, properties)
        values, atoms = self._check_qmcfc_properties(values, atoms)

        topology = self._get_topology(atoms, self.topology)

        return AtomicSystem(
            topology=topology,
            pos=values["pos"],
            vel=values.get("vel", None),
            forces=values.get("forces", None),
            charges=values.get("charges", None),
            energy=self._read_float_metadata(metadata, "energy"),
            virial=self._read_matrix_metadata(metadata, "virial"),
            stress=self._read_matrix_metadata(metadata, "stress"),
            cell=self._read_cell(metadata),
        )

    def _read_atom_count(  # pylint: disable=inconsistent-return-statements
        self,
        line: str
    ) -> int:
        """
        Reads the number of atoms from the first line.
        """
        try:
            return int(line.split()[0])
        except (ValueError, IndexError):
            self.logger.error(
                'Invalid number of atoms in extended xyz frame.',
                exception=FrameReaderError
            )

    def _read_metadata(self, line: str) -> dict[str, str]:
        """
        Parses extended xyz key-value metadata from the comment line.
        """
        metadata = {}
        tokens = shlex.split(line)
        i = 0
        while i < len(tokens):
            token = tokens[i]
            if i + 2 < len(tokens) and tokens[i + 1] == "=":
                metadata[token.lower()] = tokens[i + 2]
                i += 3
            elif "=" in token and token != "=":
                key, value = token.split("=", 1)
                if value == "" and i + 1 < len(tokens):
                    i += 1
                    value = tokens[i]
                metadata[key.lower()] = value
                i += 1
            else:
                i += 1

        return metadata

    def _read_properties(
        self,
        metadata: dict[str, str],
    ) -> dict[str, tuple[int, int]]:
        """
        Reads the extxyz Properties metadata.
        """
        if "properties" not in metadata:
            self.logger.error(
                "Extended xyz frame does not define Properties metadata.",
                exception=FrameReaderError
            )

        property_values = metadata["properties"].split(":")
        if len(property_values) % 3 != 0:
            self._raise_invalid_properties()

        properties = {}
        column_index = 0
        for i in range(0, len(property_values), 3):
            name = property_values[i].lower()

            try:
                column_count = int(property_values[i + 2])
            except ValueError:
                self._raise_invalid_properties()

            properties[name] = (column_index, column_count)
            column_index += column_count

        if "species" not in properties or "pos" not in properties:
            self.logger.error(
                "Extended xyz frame requires species and pos Properties.",
                exception=FrameReaderError
            )

        if properties["species"][1] != 1 or properties["pos"][1] != 3:
            self._raise_invalid_properties()

        return properties

    def _raise_invalid_properties(self) -> None:
        """
        Raises a consistent error for invalid Properties metadata.
        """
        self.logger.error(
            "Invalid Properties metadata in extended xyz frame.",
            exception=FrameReaderError
        )

    def _read_property_values(
        self,
        lines: list[str],
        n_atoms: int,
        properties: dict[str, tuple[int, int]],
    ) -> tuple[dict[str, np.ndarray], list[str]]:
        """
        Reads atom and numeric property values from frame lines.
        """
        atoms = []
        values = {
            "pos": np.zeros((n_atoms, 3)),
        }

        optional_properties = {
            "vel": ("vel", 3),
            "velocity": ("vel", 3),
            "forces": ("forces", 3),
            "force": ("forces", 3),
            "charge": ("charges", 1),
            "charges": ("charges", 1),
        }

        for extxyz_name, (property_name, column_count) in optional_properties.items():
            if extxyz_name in properties:
                if properties[extxyz_name][1] != column_count:
                    self._raise_invalid_properties()

                values[property_name] = self._empty_property_array(
                    n_atoms,
                    column_count,
                )

        for i in range(n_atoms):
            try:
                line_values = lines[2 + i].split()
            except IndexError:
                self.logger.error(
                    "Invalid atom line in extended xyz frame.",
                    exception=FrameReaderError
                )

            if len(line_values) < self._properties_column_count(properties):
                self.logger.error(
                    "Invalid atom line in extended xyz frame.",
                    exception=FrameReaderError
                )

            atoms.append(line_values[properties["species"][0]])
            values["pos"][i] = self._read_numeric_property(
                line_values,
                properties["pos"],
            )

            for extxyz_name, (property_name, _) in optional_properties.items():
                if extxyz_name in properties:
                    values[property_name][i] = self._read_numeric_property(
                        line_values,
                        properties[extxyz_name],
                    )

        return values, atoms

    def _check_qmcfc_properties(
        self,
        values: dict[str, np.ndarray],
        atoms: list[str],
    ) -> tuple[dict[str, np.ndarray], list[str]]:
        """
        Removes the QMCFC dummy atom and matching per-atom values.
        """
        if self.md_format != MDEngineFormat.QMCFC:
            return values, atoms

        if not atoms or atoms[0].upper() != "X":
            self.logger.error(
                (
                    'The first atom in one of the frames is not X. '
                    'Please use PQ (default) md engine instead'
                ),
                exception=FrameReaderError
            )

        return {key: value[1:] for key, value in values.items()}, atoms[1:]

    def _properties_column_count(
        self,
        properties: dict[str, tuple[int, int]]
    ) -> int:
        """
        Gets the total number of columns described by Properties metadata.
        """
        return max(
            column_index + column_count
            for column_index, column_count in properties.values()
        )

    def _empty_property_array(
        self,
        n_atoms: int,
        column_count: int,
    ) -> np.ndarray:
        """
        Creates an empty property array.
        """
        if column_count == 1:
            return np.zeros(n_atoms)

        return np.zeros((n_atoms, column_count))

    def _read_numeric_property(
        self,
        line_values: list[str],
        property_columns: tuple[int, int],
    ) -> float | list[float]:
        """
        Reads a numeric property from atom line values.
        """
        column_index, column_count = property_columns
        try:
            values = [
                float(value)
                for value in line_values[column_index:column_index + column_count]
            ]
        except ValueError:
            self.logger.error(
                "Invalid numeric value in extended xyz frame.",
                exception=FrameReaderError
            )

        if column_count == 1:
            return values[0]

        return values

    def _read_cell(self, metadata: dict[str, str]) -> Cell:
        """
        Reads the cell from Lattice metadata.
        """
        if "lattice" not in metadata:
            return Cell()

        try:
            lattice = np.array(
                [float(value) for value in metadata["lattice"].split()]
            ).reshape((3, 3))
        except ValueError:
            self.logger.error(
                "Invalid Lattice metadata in extended xyz frame.",
                exception=FrameReaderError
            )

        return self._cell_from_lattice(lattice)

    def _cell_from_lattice(self, lattice: np.ndarray) -> Cell:
        """
        Builds a cell from the three extended xyz lattice vectors.
        """
        vector_a, vector_b, vector_c = lattice

        a = np.linalg.norm(vector_a)
        b = np.linalg.norm(vector_b)
        c = np.linalg.norm(vector_c)

        alpha = self._angle_between(vector_b, vector_c)
        beta = self._angle_between(vector_a, vector_c)
        gamma = self._angle_between(vector_a, vector_b)

        return Cell(a, b, c, alpha, beta, gamma)

    def _angle_between(self, vector_a: np.ndarray, vector_b: np.ndarray) -> float:
        """
        Calculates the angle between two vectors in degrees.
        """
        denominator = np.linalg.norm(vector_a) * np.linalg.norm(vector_b)
        if denominator == 0.0:
            self.logger.error(
                "Invalid Lattice metadata in extended xyz frame.",
                exception=FrameReaderError
            )

        cosine = np.dot(vector_a, vector_b) / denominator
        return np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0)))

    def _read_float_metadata(  # pylint: disable=inconsistent-return-statements
        self,
        metadata: dict[str, str],
        key: str,
    ) -> float | None:
        """
        Reads a float metadata value if present.
        """
        if key not in metadata:
            return None

        try:
            return float(metadata[key])
        except ValueError:
            self.logger.error(
                f"Invalid {key} metadata in extended xyz frame.",
                exception=FrameReaderError
            )

    def _read_matrix_metadata(  # pylint: disable=inconsistent-return-statements
        self,
        metadata: dict[str, str],
        key: str,
    ) -> np.ndarray | None:
        """
        Reads a 3x3 matrix metadata value if present.
        """
        if key not in metadata:
            return None

        try:
            return np.array(
                [float(value) for value in metadata[key].split()]
            ).reshape((3, 3), order="F")
        except ValueError:
            self.logger.error(
                f"Invalid {key} metadata in extended xyz frame.",
                exception=FrameReaderError
            )


XYZ_FRAME_READER_TRAJ_FORMATS = (
    TrajectoryFormat.XYZ,
    TrajectoryFormat.VEL,
    TrajectoryFormat.FORCE,
    TrajectoryFormat.CHARGE,
)


def get_frame_reader(
    traj_format: TrajectoryFormat | str,
    md_format: MDEngineFormat | str = MDEngineFormat.PQ,
):
    """
    Return the frame reader implementation for a trajectory format.
    """

    traj_format = TrajectoryFormat(traj_format)
    if traj_format in XYZ_FRAME_READER_TRAJ_FORMATS:
        return XYZFrameReader(md_format=md_format)

    if traj_format == TrajectoryFormat.EXTXYZ:
        return ExtXYZFrameReader(md_format=md_format)

    # This should never happen - only for safety
    message = f'Invalid TrajectoryFormat given. {traj_format=}'
    BaseFrameReader.logger.error(message, exception=FrameReaderError)
    raise FrameReaderError(message)
