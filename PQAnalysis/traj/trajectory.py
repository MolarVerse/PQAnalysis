"""
A module containing the Trajectory class.
"""

from __future__ import annotations

import numpy as np

from beartype.typing import List, Any, Iterable

from PQAnalysis.topology import Topology
from PQAnalysis.types import Np2DNumberArray, Np1DNumberArray
from PQAnalysis.core import Cell
from PQAnalysis.atomic_system import AtomicSystem


class Trajectory:
    """
    A trajectory object is a sequence of frames.

    It can be indexed, iterated over, and added to another trajectory.
    The length of a trajectory is the number of frames in the trajectory.
    A frame can be checked for membership in a trajectory.

    This trajectory class can only handle constant topologies
    i.e. all frames in the trajectory must have the same topology.
    """

    def __init__(self,
                 frames: List[AtomicSystem] | AtomicSystem | None = None
                 ) -> None:
        """
        Parameters
        ----------
        frames : AtomicSystem | None, optional
            The list of atomic systems in the trajectory. 
            If frames is an AtomicSystem, it is first converted to list of frames.
            If frames is None, an empty list is created, by default None
        """
        if frames is None:
            frames = []

        self._frames = list(np.atleast_1d(frames))

    @property
    def box_lengths(self) -> Np2DNumberArray:
        """Np2DNumberArray: The box lengths of the trajectory."""
        return np.array([frame.cell.box_lengths for frame in self.frames])

    @property
    def box_volumes(self) -> Np1DNumberArray:
        """Np1DNumberArray: The box volumes of the trajectory."""
        return np.array([frame.cell.volume for frame in self.frames])

    def check_pbc(self) -> bool:
        """
        Checks if one cell of the trajectory is Cell().

        Returns
        -------
        bool
            False if one cell of the trajectory is Cell(), True otherwise.
        """

        if len(self.frames) == 0:
            return False

        return all(frame.pbc for frame in self.frames)

    def check_vacuum(self) -> bool:
        """
        Checks if all frames of the trajectory are in vacuum i.e. cell = Cell().

        Returns
        -------
        bool
            True if all frames of the trajectory are in vacuum, False otherwise.
        """

        return not any(frame.pbc for frame in self.frames)

    def append(self, frame: AtomicSystem) -> None:
        """
        Appends a frame to the trajectory.

        Parameters
        ----------
        frame : AtomicSystem
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

    def __getitem__(self, key: int | slice) -> AtomicSystem | Trajectory:
        """
        This method allows a frame or a trajectory to be retrieved from the trajectory.

        For example, if traj is a trajectory, then traj[0] is the first frame 
        of the trajectory. If traj is a trajectory, then traj[0:2] is a trajectory
        containing the first two frames of the trajectory.

        Parameters
        ----------
        key : int | slice
            The index or slice to retrieve the frame or trajectory from.

        Returns
        -------
        AtomicSystem | Trajectory
            The frame or trajectory retrieved from the trajectory.
        """
        frames = self.frames[key]

        if np.shape(frames) != ():
            traj = Trajectory(frames)
            traj.topology = self.topology
            return traj

        frames.topology = self.topology
        return frames

    def __iter__(self) -> Iterable[AtomicSystem]:
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

    def window(self, step: int, gap: int = 1, start: int = 0, stop: int = 0) -> Iterable[AtomicSystem]:
        """
        This method allows a window of the trajectory to be retrieved.
        Window is a sequence of frames from start to stop with a step size and a gap size.

        Parameters
        ----------
        step : int, optional
            The step of the window
        gap : int, optional
            The gap of the window, by default 1.
        start : int, optional
            The start of the window, by default 0.
        stop : int, optional
            The stop of the window, by default len(self.frames).

        Raises
        ------
        IndexError
            If start is greater than the length of the trajectory.
            If stop is greater than the length of the trajectory.
            If step is less than 1.
            If gap is less than 1.
            If start is less than 0.
            If stop is less than 0.
            If step is greater than stop - start.

        Returns
        -------
        Iterable[Frame]
            An Iterable over the frames in the window.
        """

        # If stop is not provided, set it to the length of the trajectory
        if stop == 0:
            stop = len(self)

        # If start is greater than the length of the trajectory, raise an IndexError
        if start > len(self):
            raise IndexError("start index is greater than the length of the trajectory")
        
        # If stop is greater than the length of the trajectory, raise an IndexError
        if stop > len(self):
            raise IndexError("stop index is greater than the length of the trajectory")
        
        # If step is less than 1, raise an IndexError
        if step == 0:
            raise IndexError("step must be greater than 0")
        
        # If gap is less than 1, raise an IndexError
        if gap < 1:
            raise IndexError("gap must be greater than 0")
        
        # If start is less than 0, raise an IndexError
        if start < 0:
            raise IndexError("start index must be greater than or equal to 0")
        
        # If stop is less than 0, raise an IndexError
        if stop < 0:
            raise IndexError("stop index must be greater than or equal to 0")
        
        # If step is greater than stop - start, raise an IndexError
        if step > stop - start:
            raise IndexError("step must be less than or equal to stop - start")
            
        # generate the window
        for i in range(start, stop, gap):
            if i + step > len(self):
                break
            yield self.frames[i:i + step]

    def __contains__(self, item: AtomicSystem) -> bool:
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
    def frames(self) -> List[AtomicSystem]:
        """List[AtomicSystem]: The frames in the trajectory."""
        return self._frames

    @frames.setter
    def frames(self, frames: List[AtomicSystem]) -> None:
        self._frames = frames

    @property
    def topology(self) -> Topology:
        """Topology: The topology of the trajectory."""

        if len(self.frames) == 0:
            return Topology()

        topology = self.frames[0].topology

        return topology

    @topology.setter
    def topology(self, topology: Topology) -> None:
        if len(self.frames) != 0:
            self.frames[0].topology = topology

    @property
    def cells(self) -> List[Cell]:
        """List[Cell]: The cells of the trajectory."""
        return [frame.cell for frame in self.frames]
