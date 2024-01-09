"""
A module containing classes for reading a trajectory from a file.

...

Classes
-------
TrajectoryReader
    A class for reading a trajectory from a file.
"""

import sys

from beartype.typing import List
from tqdm.auto import tqdm

import PQAnalysis.config as config

from . import BaseReader, FrameReader
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

    def __init__(self, filename: str | List[str], format: TrajectoryFormat | str = TrajectoryFormat.XYZ) -> None:
        """
        Initializes the TrajectoryReader with the given filename.

        Parameters
        ----------
        filename : str or list of str
            The name of the file to read from or a list of filenames to read from.
        """
        super().__init__(filename)
        self.frames = []
        self.format = format

    def read(self, md_format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF, constant_topology: bool = True) -> Trajectory:
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
        self.constant_topology = constant_topology
        self.topology = None
        self.frame_reader = FrameReader(md_format=md_format)
        self.md_format = md_format

        if not self.multiple_files:
            return self._read_single_file()

        else:

            traj = Trajectory()
            for filename in self.filenames:
                traj += TrajectoryReader(filename,
                                         self.format).read(self.md_format)

            return traj

    def _read_single_file(self) -> Trajectory:
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

        with open(self.filename, 'r') as f:
            sum_lines = sum(1 for _ in f)

        with open(self.filename, 'r') as f:

            # Concatenate lines of the same frame
            frame_lines = []
            for line in tqdm(f, total=sum_lines, disable=not config.with_progress_bar):
                stripped_line = line.strip()
                if stripped_line == '' or not stripped_line[0].isdigit():
                    frame_lines.append(line)
                else:
                    if frame_lines:
                        self.frames.append(self._read_single_frame(
                            ''.join(frame_lines), self.topology))

                        if self.constant_topology:
                            self.topology = self.frames[-1].topology

                    frame_lines = [line]

            if frame_lines:
                self.frames.append(self._read_single_frame(
                    ''.join(frame_lines), self.topology))

        if self.constant_topology:
            self.topology = self.frames[0].topology

        def get_size(obj):
            size = sys.getsizeof(obj)
            if isinstance(obj, dict):
                size += sum([get_size(v) for v in obj.values()])
                size += sum([get_size(k) for k in obj.keys()])
            elif hasattr(obj, '__dict__'):
                size += get_size(obj.__dict__)
            elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
                size += sum([get_size(i) for i in obj])
            return size

        traj = Trajectory(self.frames)
        print(get_size(self.frames[0]) / 10**6)
        print(get_size(self.frames[1]) / 10**6)
        print(get_size(self.frames[0].system.pos) / 10**6)
        print(
            get_size(self.frames[0].system.topology) / 10**6)
        print(
            get_size(self.frames[1].system.topology) / 10**6)
        print(
            get_size(Topology()) / 10**6)
        self.frames = []

        return traj

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
        frame = self.frame_reader.read(
            frame_string, format=self.format, topology=topology)

        # If the read frame does not have cell information, use the cell information of the previous frame
        if len(self.frames) > 0 and frame.cell.is_vacuum:
            frame.cell = self.frames[-1].cell

        return frame
