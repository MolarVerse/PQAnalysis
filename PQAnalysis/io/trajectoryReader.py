"""
A module containing classes for reading a trajectory from a file.

...

Classes
-------
TrajectoryReader
    A class for reading a trajectory from a file.
"""

from .base import BaseReader
from .frameReader import FrameReader
from ..traj.trajectory import Trajectory
from ..traj.formats import TrajectoryFormat, MDEngineFormat
from ..core.cell import Cell


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

    def __init__(self, filename: str, format: TrajectoryFormat | str = TrajectoryFormat.XYZ) -> None:
        """
        Initializes the TrajectoryReader with the given filename.

        Parameters
        ----------
        filename : str
            The name of the file to read from.
        """
        super().__init__(filename)
        self.frames = []
        self.format = format

    def read(self, md_format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF) -> Trajectory:
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
                        self._read_single_frame(
                            frame_string, frame_reader, md_format)

                    frame_string = line
                else:
                    frame_string += line

            self._read_single_frame(frame_string, frame_reader, md_format)

        return Trajectory(self.frames)

    def _read_single_frame(self, frame_string: str, frame_reader: FrameReader,  md_format: MDEngineFormat | str) -> None:
        """
        Reads a single frame from the given string.

        Parameters
        ----------
        frame_string : str
            The string containing the frame information.
        """
        frame = frame_reader.read(frame_string, format=self.format)

        # to make sure X particle is not included in the trajectory for QMCFC
        if MDEngineFormat(md_format) == MDEngineFormat.PIMD_QMCF:
            self.frames.append(frame)
        elif MDEngineFormat(md_format) == MDEngineFormat.QMCFC:
            if frame.atoms[0].name.upper() != 'X':
                raise ValueError(
                    "The first atom in one of the frames is not X. Please use pimd_qmcf (default) md engine instead")
            else:
                self.frames.append(frame[1:])

        # If the read frame does not have cell information, use the cell information of the previous frame
        if len(self.frames) > 1 and self.frames[-1].cell == Cell():
            self.frames[-1].cell = self.frames[-2].cell
