"""
A module containing classes for reading a trajectory from a file.

...

Classes
-------
TrajectoryReader
    A class for reading a trajectory from a file.
FrameReader
    A class for reading a frame from a string.
"""

import numpy as np

from beartype.typing import Tuple, List

from .base import BaseReader
from ..traj.frame import Frame
from ..traj.trajectory import Trajectory
from ..core.cell import Cell
from ..core.atom import Atom
from ..core.atomicSystem import AtomicSystem
from ..utils.exceptions import ElementNotFoundError
from ..utils.mytypes import Numpy2DFloatArray


class TrajectoryReader(BaseReader):
    """
    A class for reading a trajectory from a file.

    Inherited from BaseReader.

    ...

    Attributes
    ----------
    filename : str
        The name of the file to read from.
    frames : list of Frame
        The list of frames read from the file.
    """

    def __init__(self, filename: str) -> None:
        """
        Initializes the TrajectoryReader with the given filename.

        Parameters
        ----------
        filename : str
            The name of the file to read from.
        """
        super().__init__(filename)
        self.frames = []

    def read(self) -> Trajectory:
        """
        Reads the trajectory from the file.

        It reads the trajectory from the file and concatenates the lines of the same frame.
        The frame information is then read from the concatenated string with the FrameReader class and
        a Frame object is created.

        In order to read the cell information given in the file, the cell information of the last frame is used for
        all following frames that do not have cell information.

        Returns
        -------
        Trajectory
            The trajectory read from the file.
        """
        frame_reader = FrameReader()
        with open(self.filename, 'r') as f:

            # Concatenate lines of the same frame
            frame_string = ''
            for line in f:
                if line.strip() == '':
                    frame_string += line
                elif line.split()[0].isdigit():
                    if frame_string != '':
                        self.frames.append(frame_reader.read(frame_string))
                    frame_string = line
                else:
                    frame_string += line

            # Read the last frame and append it to the list of frames
            self.frames.append(frame_reader.read(frame_string))

            # If the read frame does not have cell information, use the cell information of the previous frame
            if self.frames[-1].cell is None:
                self.frames[-1].cell = self.frames[-2].cell

        return Trajectory(self.frames)


class FrameReader:
    """
    FrameReader reads a frame from a string.
    """

    def read(self, frame_string: str) -> Frame:
        """
        Reads a frame from a string.

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

        n_atoms, cell = self.__read_header_line__(header_line)

        xyz, atoms = self.__read_xyz__(splitted_frame_string, n_atoms)

        try:
            atoms = [Atom(atom) for atom in atoms]
        except ElementNotFoundError:
            atoms = [Atom(atom, use_guess_element=False) for atom in atoms]

        return Frame(AtomicSystem(atoms=atoms, pos=xyz, cell=cell))

    def __read_header_line__(self, header_line: str) -> Tuple[int, Cell | None]:
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

    def __read_xyz__(self, splitted_frame_string: List[str], n_atoms: int) -> Tuple[Numpy2DFloatArray, List[str]]:
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
