"""
A module containing the Momentum class. The Momentum class is used
to calculate the norm of the total linear momentum of a selection
of atoms for every frame of a velocity trajectory. It can be used
to check a simulation for center of mass drift.
"""

import itertools
import logging

# 3rd party imports
import numpy as np

from beartype.typing import Generator
from tqdm.auto import tqdm

# local absolute imports
from PQAnalysis import config
from PQAnalysis.types import (
    Np1DNumberArray,
    Np2DNumberArray,
    PositiveReal,
)
from PQAnalysis.traj import Trajectory, TrajectoryFormat
from PQAnalysis.topology import Selection, SelectionCompatible
from PQAnalysis.utils import timeit_in_class
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.io import RawTrajectoryReader, TrajectoryReader
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

# local relative imports
from .exceptions import MomentumError



class Momentum:

    """
    A class for calculating the norm of the total linear momentum
    of a selection of atoms for every frame of a velocity trajectory.

    For every frame the total linear momentum
    ``P = sum_i m_i * v_i`` is accumulated in float64 over all
    selected atoms and the norm ``|P|`` is multiplied by a scaling
    factor. With velocities in Angstrom/s (PQ velocity trajectories)
    the default scaling factor of 1e-15 converts the momentum norm
    from amu*Angstrom/s to amu*Angstrom/fs.

    Note that the velocities are parsed from file in single precision
    by the TrajectoryReader, so reported norms below roughly
    ``1e-7 * sum_i m_i * |v_i| * scale`` are parsing noise, not
    physical center of mass drift. The legacy ``equipartition.jl``
    tool parses the velocities in double precision and therefore
    resolves correspondingly smaller drift for momentum-conserving
    trajectories.

    The Momentum class can be initialized with either a trajectory
    object or via a TrajectoryReader object. If a trajectory object
    is given, it is assumed to have a constant topology over all
    frames! The main difference between the two is that the
    TrajectoryReader object allows for lazy loading of the
    trajectory, meaning that the trajectory is only loaded frame by
    frame when needed. This can be useful for large trajectories
    that do not fit into memory.

    When initialized with a TrajectoryReader of a velocity
    trajectory, the frames are streamed through the raw fast-path
    reader
    (:py:class:`~PQAnalysis.io.traj_file.raw_frame_reader.RawTrajectoryReader`)
    without building an AtomicSystem per frame. The computed
    momentum norms are bit-identical to the AtomicSystem based
    stream.
    """

    _scale_default = 1e-15
    _use_full_atom_default = False

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        traj: Trajectory | TrajectoryReader,
        selection: SelectionCompatible = None,
        use_full_atom_info: bool | None = False,
        scale: PositiveReal | None = None,
    ):
        """
        Parameters
        ----------
        traj : Trajectory | TrajectoryReader
            The velocity trajectory to analyze. If a TrajectoryReader
            is provided, the trajectory is read frame by frame via a
            frame_generator.
        selection : SelectionCompatible, optional
            The selection of atoms to include in the total momentum,
            by default None (all atoms).
        use_full_atom_info : bool | None, optional
            Whether to use the full atom information of the trajectory
            for the selection or not, by default None (False).
        scale : PositiveReal | None, optional
            The scaling factor applied to the momentum norm before
            output, by default None (1e-15, which converts
            amu*Angstrom/s to amu*Angstrom/fs).

        Raises
        ------
        MomentumError
            If the trajectory is empty.
        MomentumError
            If the selection does not select any atoms.
        MomentumError
            If the mass of an atom of the selection is unknown.
        """

        if use_full_atom_info is None:
            self.use_full_atom_info = self._use_full_atom_default
        else:
            self.use_full_atom_info = use_full_atom_info

        if scale is None:
            self.scale = self._scale_default
        else:
            self.scale = scale

        self.selection = Selection(selection)

        self._raw_reader = None
        self.frame_generator = None
        self._n_frames_total = None

        if (
            isinstance(traj, TrajectoryReader)
            and traj.traj_format == TrajectoryFormat.VEL
        ):
            # additive fast path: stream the raw float32 values of
            # the velocity trajectory without building an
            # AtomicSystem per frame (bit-identical values)
            self._raw_reader = RawTrajectoryReader(
                traj.filenames,
                traj_format=traj.traj_format,
                md_format=traj.md_format,
            )
            self._n_frames_total = self._raw_reader.count_frames()
            self.first_frame = self._raw_reader.read_first_frame()
        elif isinstance(traj, TrajectoryReader):
            # lazy loading of trajectory from file(s)
            self._n_frames_total = sum(
                traj.calculate_number_of_frames_per_file()
            )
            self.frame_generator = traj.frame_generator()
            self.first_frame = next(self.frame_generator)
        elif len(traj) > 0:
            # use trajectory object as iterator
            self._n_frames_total = len(traj)
            self.frame_generator = iter(traj)
            self.first_frame = next(self.frame_generator)
        else:
            self.logger.error(
                "Trajectory cannot be of length 0.",
                exception=MomentumError
            )

        if traj.topology is not None:
            self.topology = traj.topology
        else:
            self.topology = self.first_frame.topology

        self.indices = self.selection.select(
            self.topology,
            self.use_full_atom_info
        )

        if len(self.indices) == 0:
            self.logger.error(
                "The selection does not select any atoms.",
                exception=MomentumError
            )

        masses = [self.topology.atoms[index].mass for index in self.indices]

        if any(mass is None for mass in masses):
            self.logger.error(
                (
                "The mass of at least one selected atom is unknown. "
                "The total momentum cannot be calculated."
                ),
                exception=MomentumError
            )

        self.masses = np.asarray(masses, dtype=np.float64)

        self.momentum_norms = np.array([])

    @timeit_in_class
    def run(self) -> Np1DNumberArray:
        """
        Runs the momentum analysis.

        For every frame of the trajectory the total linear momentum
        of the selected atoms is accumulated in float64 and the
        scaled norm of the momentum vector is stored.

        This method will display a progress bar by default.
        This can be disabled by setting with_progress_bar to
        False.

        Returns
        -------
        Np1DNumberArray
            The scaled norms of the total linear momentum, one value
            per frame.

        Raises
        ------
        MomentumError
            If a frame does not contain velocity information for all
            atoms of the topology.
        """
        norms = []
        selected_masses = self.masses[:, None]

        for velocities in tqdm(
            self._velocities(),
            total=self._n_frames_total,
            disable=not config.with_progress_bar):

            if velocities.shape[0] != self.topology.n_atoms:
                self.logger.error(
                    (
                    "The trajectory does not contain velocity "
                    "information for all atoms. Please provide a "
                    "velocity trajectory."
                    ),
                    exception=MomentumError
                )

            momentum = np.sum(
                selected_masses * velocities[self.indices],
                axis=0
            )

            norms.append(float(np.linalg.norm(momentum)) * self.scale)

        self.momentum_norms = np.array(norms, dtype=np.float64)

        return self.momentum_norms

    def _velocities(self) -> Generator[Np2DNumberArray, None, None]:
        """
        Yields the float64 velocities of all frames.

        Dispatches to the raw fast-path stream if the analysis was
        constructed from a TrajectoryReader of a velocity trajectory
        and to the AtomicSystem based stream otherwise. Both streams
        yield bit-identical float64 arrays.

        Yields
        ------
        Np2DNumberArray
            The velocities of one frame with shape ``(n_atoms, 3)``.
        """
        if self._raw_reader is not None:
            for values, _cell in self._raw_reader.raw_frame_generator():
                yield np.asarray(values, dtype=np.float64)
        else:
            frames = itertools.chain(
                [self.first_frame], self.frame_generator
            )

            for frame in frames:
                yield np.asarray(frame.vel, dtype=np.float64)

    @property
    def n_frames(self) -> int:
        """int: The number of analyzed frames after calling run()."""
        return len(self.momentum_norms)
