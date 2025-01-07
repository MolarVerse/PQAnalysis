"""
A module containing the Trajectory class.
"""

import logging
import numpy as np

from beartype.typing import List, Any, Iterable

from PQAnalysis.topology import Topology
from PQAnalysis.exceptions import PQIndexError, PQTypeError
from PQAnalysis.core import Cell
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.types import (
    Np2DNumberArray,
    Np1DNumberArray,
    PositiveReal,
    Bool,
)
from PQAnalysis.type_checking import (
    runtime_type_checking,
    runtime_type_checking_setter,
)



class Trajectory:

    """
    A trajectory object is a sequence of frames.

    It can be indexed, iterated over, and added to another trajectory.
    The length of a trajectory is the number of frames in the trajectory.
    A frame can be checked for membership in a trajectory.

    This trajectory class can only handle constant topologies
    i.e. all frames in the trajectory must have the same topology.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)

    @runtime_type_checking
    def __init__(
        self,
        frames: List[AtomicSystem] | AtomicSystem | None = None,
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

        if isinstance(frames, AtomicSystem):
            frames = [frames]

        self._frames = frames.copy()

        self.logger = setup_logger(self.logger)

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

    @runtime_type_checking
    def append(self, frame: AtomicSystem) -> None:
        """
        Appends a frame to the trajectory.

        Parameters
        ----------
        frame : AtomicSystem
            The frame to append.
        """

        self.frames.append(frame)

    @runtime_type_checking
    def pop(self, index: int = -1) -> AtomicSystem:
        """
        Removes a frame from the trajectory at the specified index.

        Parameters
        ----------
        index : int
            The index of the frame to remove.

        Returns
        -------
        AtomicSystem
            The frame removed from the trajectory.
        """

        return self.frames.pop(index)

    def copy(self) -> "Trajectory":
        """
        Returns a copy of the trajectory.

        Returns
        -------
        Trajectory
            A copy of the trajectory.
        """

        return Trajectory(self.frames.copy())

    def __len__(self) -> int:
        """
        This method allows the length of a trajectory to be computed.

        Returns
        -------
        int
            The number of frames in the trajectory.
        """
        return len(self.frames)

    def __getitem__(self, key: int | slice) -> "AtomicSystem | Trajectory":
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

        if isinstance(frames, AtomicSystem):
            frames.topology = self.topology
            return frames

        traj = Trajectory(frames)
        traj.topology = self.topology
        return traj

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

    @runtime_type_checking
    def window(
        self,
        window_size: int,
        window_gap: int = 1,
        trajectory_start: int = 0,
        trajectory_stop: int | None = None,
    ) -> Iterable["Trajectory"]:
        """
        This method allows a window of the trajectory to be retrieved.
        Window is a sequence of frames from start to stop with a window size and a gap size.

        Parameters
        ----------
        window_size : int
            The size of the window.
        window_gap : int, optional
            The gap size between two windows, by default 1
        trajectory_start : int, optional
            The start index of the first window, by default 0
        trajectory_stop : int | None, optional
            Stop index of the window generator, by default None, which then
            set to the length of the trajectory.

        Raises
        ------
        PQIndexError
            If window_start is less than 0 or greater than the length of the trajectory.
            If window_stop is less than 0 or greater than the length of the trajectory.
            If window_size is less than 1 or greater than the length of the trajectory.
            If window_gap is less than 1 or greater than the length of the trajectory.
            If window_size is greater than trajectory_stop - trajectory_start.

        Warning
        -------
        If not all frames are included in the windows, a warning is issued.

        Yields
        ------
        Iterable[Trajectory]
            An iterable over the windows of the trajectory with the specified window size and gap.
        """

        # If trajectory_stop is not provided, set it to the length of the trajectory
        if trajectory_stop is None:
            trajectory_stop = len(self)

        # If trajectory_start is less than 0 or greater than the
        # length of the trajectory, raise an IndexError
        if trajectory_start < 0 or trajectory_start > len(self):
            self.logger.error(
                "start index is less than 0 or greater than the length of the trajectory",
                exception=PQIndexError,
            )

        # If trajectory_stop is less than 0 or greater than the
        # length of the trajectory, raise an IndexError
        if trajectory_stop < 0 or trajectory_stop > len(self):
            self.logger.error(
                "stop index is less than 0 or greater than the length of the trajectory",
                exception=PQIndexError,
            )

        # If window_step is less than 1 or greater than
        # the length of the trajectory, raise an IndexError
        if window_size < 1 or window_size > len(self):
            self.logger.error(
                "window size can not be less than 1 or greater than the length of the trajectory",
                exception=PQIndexError,
            )

        # If window_gap is less than 1 or greater than
        # the length of the trajectory, raise an IndexError
        if window_gap < 1 or window_gap > len(self):
            self.logger.error(
                "window gap can not be less than 1 or greater than the length of the trajectory",
                exception=PQIndexError,
            )

        # If trajectory_start is greater than or equal to trajectory_stop, raise an IndexError
        if trajectory_start >= trajectory_stop:
            self.logger.error(
                "start index is greater than or equal to the stop index",
                exception=PQIndexError,
            )

        # If window_size is greater than trajectory_stop - trajectory_start, raise an IndexError
        if window_size > trajectory_stop - trajectory_start:
            self.logger.error(
                "window size is greater than the trajectory_stop - trajectory_start",
                exception=PQIndexError,
            )

        # Check if all frames are included in the windows
        # Length of the trajectory - window_size should be divisible by window_gap
        if (
            ((trajectory_stop - trajectory_start) - window_size) % window_gap
            != 0
        ):
            self.logger.warning(
                "Not all frames are included in the windows. Check the window size and gap."
            )

        # generate the window of the trajectory
        for i in range(
            trajectory_start,
            trajectory_stop - window_size + 1,
            window_gap,
        ):
            yield self[i:i + window_size]

    @runtime_type_checking
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

    def __add__(self, other: "Trajectory") -> "Trajectory":
        """
        This method allows two trajectories to be added together.

        Parameters
        ----------
        other : Trajectory
            The other trajectory to add.
        """

        if not isinstance(other, Trajectory):
            self.logger.error(
                "Only Trajectory objects can be added to a Trajectory.",
                exception=PQTypeError,
            )

        return Trajectory(self.frames + other.frames)

    def __eq__(self, other: Any) -> Bool:
        """
        This method allows two trajectories to be compared for equality.

        Parameters
        ----------
        other : Trajectory
            The other trajectory to compare.

        Returns
        -------
        Bool
            Whether the two trajectories are equal.
        """

        return self.isclose(other)

    def isclose(
        self,
        other: Any,
        rtol: PositiveReal = 1e-9,
        atol: PositiveReal = 0.0,
    ) -> Bool:
        """
        This method allows two trajectories to be compared for closeness.

        Parameters
        ----------
        other : Any
            The other object to compare with.
        rtol : PositiveReal, optional
            The relative tolerance parameter, by default 1e-9
        atol : PositiveReal, optional
            The absolute tolerance parameter, by default 0.0

        Returns
        -------
        Bool
            Whether the two trajectories are close.
        """

        if not isinstance(other, Trajectory):
            return False

        if len(self.frames) != len(other.frames):
            return False

        return np.all(
            [
                self.frames[i].isclose(other.frames[i], rtol=rtol, atol=atol)
                for i in np.arange(len(self.frames))
            ]
        )

    def __str__(self) -> str:
        """
        This method allows the string representation of a trajectory to be retrieved.

        Returns
        -------
        str
            The string representation of the trajectory.
        """
        return f"Trajectory with {len(self)} frames"

    def __repr__(self) -> str:
        """
        This method allows the string representation of a trajectory to be retrieved.

        Returns
        -------
        str
            The string representation of the trajectory.
        """
        return self.__str__()

    @property
    def frames(self) -> List[AtomicSystem]:
        """List[AtomicSystem]: The frames in the trajectory."""
        return self._frames

    @frames.setter
    @runtime_type_checking_setter
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
    @runtime_type_checking_setter
    def topology(self, topology: Topology) -> None:
        if len(self.frames) != 0:
            self.frames[0].topology = topology

    @property
    def cells(self) -> List[Cell]:
        """List[Cell]: The cells of the trajectory."""
        return [frame.cell for frame in self.frames]

    @property
    def com_residue_traj(self) -> "Trajectory":
        """Trajectory: The trajectory with the center of mass of the residues."""

        frames = [frame.center_of_mass_residues for frame in self.frames]

        return Trajectory(frames)
