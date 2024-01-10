"""
A module containing classes for reading a trajectory from a file.

...

Classes
-------
TrajectoryReader
    A class for reading a trajectory from a file.
"""

from __future__ import annotations

from beartype.typing import List, Generator
from tqdm.auto import tqdm

import PQAnalysis.config as config

from . import BaseReader, FrameReader, TrajectoryReaderError
from ..traj import Trajectory, TrajectoryFormat, MDEngineFormat, Frame
from ..core import Cell
from ..topology import Topology


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

    def __init__(self,
                 filename: str | List[str],
                 format: TrajectoryFormat | str = TrajectoryFormat.XYZ,
                 md_format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF,
                 topology: Topology | None = None,
                 constant_topology: bool = True
                 ) -> None:
        """
        Initializes the TrajectoryReader with the given filename.

        Parameters
        ----------
        filename : str or list of str
            The name of the file to read from or a list of filenames to read from.
        """
        super().__init__(filename)
        if not self.multiple_files:
            self.filenames = [self.filename]

        self.frames = []
        self.format = TrajectoryFormat(format)
        self.topology = topology
        self.constant_topology = constant_topology
        self.md_format = MDEngineFormat(md_format)
        self.frame_reader = FrameReader(md_format=self.md_format)

    def read(self, topology: Topology | None = None) -> Trajectory:
        """
        Reads the trajectory from the file.

        It reads the trajectory from the file and concatenates the lines of the same frame.
        The frame information is then read from the concatenated string with the FrameReader class and
        a Frame object is created.

        In order to read the cell information given in the file, the cell information of the last frame is used for
        all following frames that do not have cell information.

        If the trajectory is split into multiple files, the files are read one after another and the frames are
        concatenated into a single trajectory.

        Parameters
        ----------
        md_format : MDEngineFormat | str, optional
            The format of the trajectory. Default is MDEngineFormat.PIMD_QMCF.
        constant_topology : bool, optional
            Whether the topology is constant over the trajectory or does change. Default is True.

        Returns
        -------
        Trajectory
            The trajectory read from the file.
        """

        self.topology = topology

        traj = Trajectory()
        for frame in self.frame_generator():
            traj.append(frame)

        return traj

    def frame_generator(self) -> Generator[Frame]:
        last_cell = None
        for filename in self.filenames:
            with open(filename, 'r') as self.file:
                sum_lines = sum(1 for _ in self.file)

            with open(filename, 'r') as self.file:
                frame_lines = []
                for line in tqdm(self.file, total=sum_lines, disable=not config.with_progress_bar):
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
        Reads the cells from the trajectory.

        Returns
        -------
        list of Cell
            The list of cells read from the trajectory.
        """
        last_cell = None
        with open(self.filenames[0], 'r') as f:
            line = f.readline()
            n_atoms = int(line.strip()[0])

        for filename in self.filenames:
            with open(filename, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    stripped_line = line.strip()
                    splitted_line = stripped_line.split()
                    if len(splitted_line) == 1 and cell is None:
                        cell = Cell()
                        if last_cell is not None:
                            cell = last_cell
                        yield cell
                    elif len(splitted_line) == 4:
                        cell = Cell(float(splitted_line[1]), float(
                            splitted_line[2]), float(splitted_line[3]))
                        yield cell
                    elif len(splitted_line) == 7:
                        cell = Cell(float(splitted_line[1]), float(splitted_line[2]), float(splitted_line[3]),
                                    float(splitted_line[4]), float(splitted_line[5]), float(splitted_line[6]))
                        yield cell
                    else:
                        raise TrajectoryReaderError(
                            f"Invalid number of arguments for box: {len(splitted_line)} encountered in file {filename} {stripped_line}.")

                    last_cell = cell

                    for _ in range(n_atoms+1):
                        next(f)

    def _read_single_frame(self,
                           frame_string: str,
                           topology: Topology | None = None
                           ) -> Frame:
        """
        Reads a single frame from the given string.

        Parameters
        ----------
        frame_string : str
            The string containing the frame information.
        topology : Topology, optional
            The topology of the frame. Default is None.

        Raises
        ------
        TrajectoryReaderError
            If the first atom in the frame is not X for QMCFC.
        """
        return self.frame_reader.read(
            frame_string, format=self.format, topology=topology)
