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

from beartype.typing import List, Iterator, Any, Iterable

from . import Frame
from ..topology import Topology
from ..types import Np2DNumberArray, Np1DNumberArray
from ..core import Cell


class Trajectory:
    """
    A trajectory object is a sequence of frames.

    It can be indexed, iterated over, and added to another trajectory.
    The length of a trajectory is the number of frames in the trajectory.
    A frame can be checked for membership in a trajectory.

    This trajectory class can only handle constant topologies i.e. all frames in the trajectory must have the same topology.
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
            self._frames = []
        else:
            self._frames = frames

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

    @property
    def box_volumes(self) -> Np1DNumberArray:
        """
        Returns the box volumes of the trajectory.

        Returns
        -------
        Np1DNumberArray
            The box volumes of the trajectory.
        """
        return np.array([frame.cell.volume for frame in self.frames])

    def check_PBC(self) -> bool:
        """
        Checks if one cell of the trajectory is Cell().

        Returns
        -------
        bool
            False if one cell of the trajectory is Cell(), True otherwise.
        """

        if len(self.frames) == 0:
            return False

        return all(frame.PBC for frame in self.frames)

    def check_vacuum(self) -> bool:
        """
        Checks if all frames of the trajectory are in vacuum i.e. cell = Cell().

        Returns
        -------
        bool
            True if all frames of the trajectory are in vacuum, False otherwise.
        """

        return not any(frame.PBC for frame in self.frames)

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
            traj = Trajectory(frames)
            traj.topology = self.topology
            return traj
        else:
            frames.topology = self.topology
            return frames

    def __iter__(self) -> Iterable[Frame]:
        """
        This method allows a trajectory to be iterated over.

        Returns
        -------
        Iterable[Frame]
            An Iterable over the frames in the trajectory.
        """
        for frame in self.frames:
            frame.topology = self.topology
            yield frame

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

    @property
    def frames(self) -> List[Frame]:
        """
        The frames in the trajectory.

        Returns
        -------
        List[Frame]
            The frames in the trajectory.
        """
        return self._frames

    @frames.setter
    def frames(self, frames: List[Frame]) -> None:
        """
        Sets the frames in the trajectory.

        Parameters
        ----------
        frames : List[Frame]
            The frames in the trajectory.
        """
        self._frames = frames

    @property
    def topology(self) -> Topology:
        """
        The topology of the trajectory.

        Returns
        -------
        Topology
            The topology of the trajectory.

        Raises
        ------
        TrajectoryError
            If the frames in the trajectory do not have the same topology.
        """

        if len(self.frames) == 0:
            return Topology()

        topology = self.frames[0].topology

        return topology

    @topology.setter
    def topology(self, topology: Topology) -> None:
        """
        Sets the topology of the trajectory.

        Parameters
        ----------
        topology : Topology
            The topology of the trajectory.
        """

        if len(self.frames) != 0:
            self.frames[0].topology = topology

    @property
    def cells(self) -> List[Cell]:
        """
        Returns the cells of the trajectory.

        Returns
        -------
        list of Cell
            The list of cells of the trajectory.
        """
        return [frame.cell for frame in self.frames]
