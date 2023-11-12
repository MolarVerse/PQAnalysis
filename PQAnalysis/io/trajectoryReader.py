"""
A module containing classes for reading a trajectory from a file.

...

Classes
-------
TrajectoryReader
    A class for reading a trajectory from a file.
"""

from .base import BaseReader
from ..traj.trajectory import Trajectory
from ..traj.formats import TrajectoryFormat
from .frameReader import FrameReader


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
                        self.frames.append(frame_reader.read(
                            frame_string, format=self.format))
                    frame_string = line
                else:
                    frame_string += line

            # Read the last frame and append it to the list of frames
            self.frames.append(frame_reader.read(frame_string))

            # If the read frame does not have cell information, use the cell information of the previous frame
            if self.frames[-1].cell is None:
                self.frames[-1].cell = self.frames[-2].cell

        return Trajectory(self.frames)
