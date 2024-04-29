"""
A module containing classes for reading a trajectory from a file.
"""

from __future__ import annotations

# 3rd party modules
from beartype.typing import List, Generator
from tqdm.auto import tqdm

# Local absolute imports
from PQAnalysis.config import with_progress_bar
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.traj import Trajectory, TrajectoryFormat, MDEngineFormat
from PQAnalysis.core import Cell
from PQAnalysis.topology import Topology

# Local relative modules
from .. import BaseReader
from .exceptions import TrajectoryReaderError
from .frame_reader import FrameReader


class TrajectoryReader(BaseReader):
    """
    A class for reading a trajectory from a file.

    Inherited from BaseReader.
    """

    def __init__(self,
                 filename: str | List[str],
                 traj_format: TrajectoryFormat | str = TrajectoryFormat.AUTO,
                 md_format: MDEngineFormat | str = MDEngineFormat.PQ,
                 topology: Topology | None = None,
                 constant_topology: bool = True
                 ) -> None:
        """
        Parameters
        ----------
        filename : str or list of str
            The name of the file to read from or a list of filenames to read from.
        traj_format : TrajectoryFormat | str, optional
            The format of the trajectory. Default is TrajectoryFormat.AUTO.
            The format is inferred from the file extension.
        md_format : MDEngineFormat | str, optional
            The format of the trajectory. Default is MDEngineFormat.PQ.
        topology : Topology, optional
            The topology of the trajectory. Default is None.
        constant_topology : bool, optional
            Whether the topology is constant over the trajectory or does change. Default is True.
        """
        super().__init__(filename)

        if not self.multiple_files:
            self.filenames = [self.filename]

        self.file = None

        self.frames = []
        self.topology = topology
        self.constant_topology = constant_topology

        self.traj_format = TrajectoryFormat((traj_format, self.filenames[0]))

        self.md_format = MDEngineFormat(md_format)
        self.frame_reader = FrameReader(md_format=self.md_format)

        # NOTE: Progress bar is disabled by default
        #       This way the frame_generator can be used in other functions
        #       without the need to disable the progress bar
        #       The global config.with_progress_bar is only set in the read function
        self.with_progress_bar = False

    def read(self, topology: Topology | None = None) -> Trajectory:
        """
        Reads the trajectory from the file.

        It reads the trajectory from the file and concatenates the lines of the
        same frame. The frame information is then read from the concatenated
        string with the FrameReader class and a Frame object is created.

        In order to read the cell information given in the file, the cell 
        information of the last frame is used for all following frames that
        do not have cell information.

        If the trajectory is split into multiple files, the files are read one 
        after another and the frames are concatenated into a single trajectory.

        Parameters
        ----------
        topology : Topology, optional
            The topology of the trajectory. Default is None.

        Returns
        -------
        Trajectory
            The trajectory read from the file.
        """

        self.with_progress_bar = with_progress_bar
        self.topology = topology

        traj = Trajectory()
        for frame in self.frame_generator():
            traj.append(frame)

        return traj

    def frame_generator(self) -> Generator[AtomicSystem]:
        """
        A generator that yields the frames of the trajectory.

        The difference to the read method is that the read method returns the whole
        trajectory at once, while the frame_generator yields the frames one after 
        another. This is useful if the trajectory is very large and cannot be stored in memory.

        This method is used to read the trajectory from the file. It reads the 
        trajectory from the file and concatenates the lines of the same frame. 
        The frame information is then read from the concatenated string with the
        FrameReader class and a Frame object is created.

        Yields
        ------
        Generator[AtomicSystem]
            The frames of the trajectory.
        """
        last_cell = None
        for filename in self.filenames:
            with open(filename, 'r', encoding='utf-8') as self.file:
                sum_lines = sum(1 for _ in self.file)

            with open(filename, 'r', encoding='utf-8') as self.file:
                frame_lines = []
                for line in tqdm(self.file, total=sum_lines, disable=not self.with_progress_bar):
                    stripped_line = line.strip()
                    if stripped_line == '' or not stripped_line[0].isdigit():
                        frame_lines.append(line)
                    else:
                        if frame_lines:
                            frame = self._read_single_frame(
                                ''.join(frame_lines), self.topology)
                            if frame.cell.is_vacuum and last_cell is not None:
                                frame.cell = last_cell
                            last_cell = frame.cell
                            yield frame
                            if self.constant_topology and self.topology is not None:
                                self.topology = frame.topology
                        frame_lines = [line]

                if frame_lines:
                    frame = self._read_single_frame(
                        ''.join(frame_lines), self.topology)
                    if frame.cell.is_vacuum and last_cell is not None:
                        frame.cell = last_cell
                    last_cell = frame.cell
                    yield frame
                if self.constant_topology and self.topology is not None:
                    self.topology = frame.topology

    def calculate_number_of_frames(self) -> int:
        """
        Calculates the number of frames in the trajectory file.

        Returns
        -------
        int
            The number of frames in the trajectory file.
        """

        n_frames = 0

        for filename in self.filenames:
            with open(filename, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                n_frames += int(len(lines) / int(lines[0].split()[0]))

        return n_frames

    @property
    def cells(self) -> list[Cell]:
        """
        Returns the cells of the trajectory.

        Returns
        -------
        list of Cell
            The list of cells of the trajectory.
        """
        return list(self._cell_generator())

    def _cell_generator(self) -> Generator[List[Cell]]:
        """
        A generator that yields the cells of the trajectory.

        This method is used to read the cells from the file. 
        It reads the cells from the file and yields them one after
        another. If the cell information is not given in the 
        file, the cell information of the last frame is used 
        for all following frames that do not have cell information.

        Yields
        ------
        list of Cell
            The list of cells read from the trajectory.
        """
        last_cell = None
        with open(self.filenames[0], 'r', encoding='utf-8') as f:
            line = f.readline()
            n_atoms = int(line.split()[0])

        for filename in self.filenames:
            with open(filename, 'r', encoding='utf-8') as f:
                for line in f:
                    stripped_line = line.strip()
                    splitted_line = stripped_line.split()

                    if len(splitted_line) == 1 and cell is None:
                        cell = Cell()

                        if last_cell is not None:
                            cell = last_cell

                        yield cell

                    elif len(splitted_line) == 4:

                        cell = Cell(
                            float(splitted_line[1]),
                            float(splitted_line[2]),
                            float(splitted_line[3])
                        )

                        yield cell

                    elif len(splitted_line) == 7:

                        cell = Cell(
                            float(splitted_line[1]),
                            float(splitted_line[2]),
                            float(splitted_line[3]),
                            float(splitted_line[4]),
                            float(splitted_line[5]),
                            float(splitted_line[6])
                        )

                        yield cell

                    else:

                        raise TrajectoryReaderError(
                            "Invalid number of arguments for box: "
                            f"{len(splitted_line)} encountered in file "
                            f"{filename} {stripped_line}."
                        )

                    last_cell = cell

                    for _ in range(n_atoms+1):
                        next(f, None)  # Skip the next n_atoms+1 lines

    def _read_single_frame(self,
                           frame_string: str,
                           topology: Topology | None = None
                           ) -> AtomicSystem:
        """
        Reads a single frame from the given string.

        Parameters
        ----------
        frame_string : str
            The string containing the frame information.
        topology : Topology, optional
            The topology of the frame. Default is None.

        Returns
        -------
        AtomicSystem
            The AtomicSystem object of the frame.

        Raises
        ------
        TrajectoryReaderError
            If the first atom in the frame is not X for QMCFC.
        """
        return self.frame_reader.read(
            frame_string, traj_format=self.traj_format, topology=topology)
