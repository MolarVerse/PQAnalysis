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
