"""
A module containing the Trajectory class.

...

Classes
-------
Trajectory
    A trajectory is a sequence of frames.
"""

from __future__ import annotations

import numpy as np

from beartype.typing import List, Iterator, Any

from . import Frame
from ..types import Np2DNumberArray


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

    def __init__(self, frames: List[Frame] = None) -> None:
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

    @property
    def box_lengths(self) -> Np2DNumberArray:
        """
        Returns the box lengths of the trajectory.

        Returns
        -------
        Np2DNumberArray
            The box lengths of the trajectory.
        """
        return np.array([frame.cell.box_lengths for frame in self.frames])

    def check_PBC(self) -> bool:
        """
        Checks if one cell of the trajectory is Cell().

        Returns
        -------
        bool
            False if one cell of the trajectory is Cell(), True otherwise.
        """

        return all(frame.PBC for frame in self.frames)

    def check_vacuum(self) -> bool:
        """
        Checks if all frames of the trajectory are in vacuum i.e. cell = Cell().

        Returns
        -------
        bool
            True if all frames of the trajectory are in vacuum, False otherwise.
        """

        return not any(frame.cell is None for frame in self.frames)

    def append(self, frame: Frame) -> None:
        """
        Appends a frame to the trajectory.

        Parameters
        ----------
        frame : Frame
            The frame to append.
        """

        self.frames.append(frame)

    def __len__(self) -> int:
        """
        This method allows the length of a trajectory to be computed.

        Returns
        -------
        int
            The number of frames in the trajectory.
        """
        return len(self.frames)

    def __getitem__(self, key: int | slice) -> Frame | Trajectory:
        """
        This method allows a frame or a trajectory to be retrieved from the trajectory.

        For example, if traj is a trajectory, then traj[0] is the first frame of the trajectory.
        If traj is a trajectory, then traj[0:2] is a trajectory containing the first two frames of the trajectory.



        Parameters
        ----------
        key : int | slice
            The index or slice to retrieve the frame or trajectory from.

        Returns
        -------
        Frame | Trajectory
            The frame or trajectory retrieved from the trajectory.
        """
        frames = self.frames[key]
        if np.shape(frames) != ():
            return Trajectory(frames)
        else:
            return frames

    def __iter__(self) -> Iterator:
        """
        This method allows a trajectory to be iterated over.

        Returns
        -------
        Iter
            The iterator over the frames in the trajectory.
        """
        return iter(self.frames)

    def __contains__(self, item: Frame) -> bool:
        """
        This method allows a frame to be checked for membership in a trajectory.

        Parameters
        ----------
        item : Frame
            The frame to check for membership in the trajectory.

        Returns
        -------
        bool
            Whether the frame is in the trajectory.
        """
        return item in self.frames

    def __add__(self, other: Trajectory) -> Trajectory:
        """
        This method allows two trajectories to be added together.

        Parameters
        ----------
        other : Trajectory
            The other trajectory to add.
        """

        return Trajectory(self.frames + other.frames)

    def __eq__(self, other: Any) -> bool:
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
