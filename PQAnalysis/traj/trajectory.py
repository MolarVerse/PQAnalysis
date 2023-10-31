"""
A module containing the Trajectory class.

...

Classes
-------
Trajectory
    A trajectory is a sequence of frames.
"""

from typing import List

from PQAnalysis.traj.frame import Frame


class Trajectory:
    """
    A trajectory object is a sequence of frames.

    It can be indexed, iterated over, and added to another trajectory.
    The length of a trajectory is the number of frames in the trajectory.
    A frame can be checked for membership in a trajectory.

    ...

    Attributes
    ----------
    frames : list of Frame
        The list of frames in the trajectory.

    """

    def __init__(self, frames: List[Frame] = None):
        """
        Initializes the Trajectory with the given frames.

        Parameters
        ----------
        frames : list of Frame, optional
            The list of frames in the trajectory.
        """
        if frames is None:
            self.frames = []
        else:
            self.frames = frames

    def check_PBC(self):
        """
        Checks if one cell of the trajectory is not None.

        Returns
        -------
        bool
            False if one cell of the trajectory is None, True otherwise.
        """

        if not all(frame.PBC for frame in self.frames):
            return False
        else:
            return True

    def append(self, frame: Frame):
        """
        Appends a frame to the trajectory.

        Parameters
        ----------
        frame : Frame
            The frame to append.
        """

        if not isinstance(frame, Frame):
            raise TypeError('only Frame can be appended to Trajectory.')

        self.frames.append(frame)

    def __len__(self) -> int:
        return len(self.frames)

    def __getitem__(self, key) -> Frame:
        return self.frames[key]

    def __iter__(self) -> iter:
        return iter(self.frames)

    def __contains__(self, item) -> bool:
        return item in self.frames

    # TODO: check if frames are compatible with each other
    def __add__(self, other) -> 'Trajectory':
        """
        This method allows two trajectories to be added together.

        Parameters
        ----------
        other : Trajectory
            The other trajectory to add.
        """

        if not isinstance(other, Trajectory):
            raise TypeError('only Trajectory can be added to Trajectory.')

        if len(self.frames) != 0 and len(other) != 0 and not self.frames[0].is_combinable(other.frames[0]):
            raise ValueError('Frames are not compatible.')

        return Trajectory(self.frames + other.frames)

    def __eq__(self, other) -> bool:
        """
        This method allows two trajectories to be compared for equality.

        Parameters
        ----------
        other : Trajectory
            The other trajectory to compare.
        """

        if not isinstance(other, Trajectory):
            return False

        if len(self.frames) != len(other.frames):
            return False

        return self.frames == other.frames
