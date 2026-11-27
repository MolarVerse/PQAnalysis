"""
A module containing a reader for ASE trajectory files.

ASE trajectories are the binary ``.traj`` files of the Atomic Simulation
Environment, see https://wiki.fysik.dtu.dk/ase/ase/io/trajectory.html.
Reading them needs the optional dependency ``ase``
(``pip install PQAnalysis[ase]``).
"""

import logging

import numpy as np

from beartype.typing import Any, Generator, List

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Cell
from PQAnalysis.topology import Topology
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis import __package_name__

from .exceptions import AseTrajReaderError
from .h5md_reader import _atom_from_name_and_species, _cell_from_box_vectors



class AseTrajReader:

    """
    A class for reading the frames of a single ASE trajectory file.

    ASE stores Angstrom, eV and its own time unit. Values are converted
    to the PQ units used in PQAnalysis: Angstrom, Angstrom s-1,
    kcal mol-1 Angstrom-1, e and kcal mol-1. The atoms are named by the
    chemical symbols of their atomic numbers.

    Generally, this class should not be used directly. Instead, use the
    :py:class:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader`
    with ``traj_format="ase"`` or a ``.traj`` file.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            The name of the ASE trajectory file.

        Raises
        ------
        AseTrajReaderError
            If ase is not installed.
        """
        self.filename = filename
        self._ase = _import_ase()

    @property
    def number_of_frames(self) -> int:
        """int: The number of frames in the file."""
        with self._open() as trajectory:
            return len(trajectory)

    @property
    def steps(self) -> List[int | None]:
        """list of int | None: The step of each frame, None if not stored."""
        with self._open() as trajectory:
            return [atoms.info.get("step") for atoms in trajectory]

    def cell_generator(self) -> Generator[Cell]:
        """
        A generator that yields the cell of each frame.

        Yields
        ------
        Cell
            The cell of the frame, a vacuum cell if the frame is not
            periodic.
        """
        with self._open() as trajectory:
            for atoms in trajectory:
                yield self._read_cell(atoms)

    @runtime_type_checking
    def frame_generator(
        self,
        start: int = 0,
        stop: int | None = None,
        topology: Topology | None = None,
    ) -> Generator[AtomicSystem]:
        """
        A generator that yields the frames of the file.

        Parameters
        ----------
        start : int, optional
            The index of the first frame, by default 0.
        stop : int | None, optional
            The index after the last frame, by default None, meaning
            all frames.
        topology : Topology | None, optional
            The topology of the frames, by default None, meaning it is
            built from the atomic numbers of the first read frame.

        Yields
        ------
        AtomicSystem
            The frames with positions in Angstrom, velocities in
            Angstrom s-1, forces in kcal mol-1 Angstrom-1, charges in e
            and the potential energy in kcal mol-1.
        """
        units = self._ase.units
        velocity_factor = units.second
        energy_factor = units.mol / units.kcal

        with self._open() as trajectory:
            if stop is None:
                stop = len(trajectory)

            for index in range(start, stop):
                atoms = trajectory[index]

                if topology is None:
                    topology = _read_topology(atoms)

                values = {"pos": _to_float64(atoms.get_positions())}

                if "momenta" in atoms.arrays:
                    values["vel"] = (
                        _to_float64(atoms.get_velocities()) * velocity_factor
                    )

                results = atoms.calc.results if atoms.calc is not None else {}

                if "forces" in results:
                    values["forces"] = (
                        _to_float64(results["forces"]) * energy_factor
                    )

                if "charges" in results:
                    values["charges"] = _to_float64(results["charges"])
                elif "initial_charges" in atoms.arrays:
                    values["charges"] = _to_float64(
                        atoms.get_initial_charges()
                    )

                if "energy" in results:
                    values["energy"] = float(results["energy"]) * energy_factor

                yield AtomicSystem(
                    topology=topology,
                    cell=self._read_cell(atoms),
                    **values,
                )

    def _open(self) -> Any:
        """
        Opens the ASE trajectory for reading.
        """
        try:
            return self._ase.io.Trajectory(self.filename, "r")
        except Exception:  # pylint: disable=broad-exception-caught
            self.logger.error(
                f"Could not open the ASE trajectory file {self.filename}.",
                exception=AseTrajReaderError,
            )

        return None  # pragma: no cover

    def _read_cell(self, atoms: Any) -> Cell:
        """
        Reads the cell of a frame, a vacuum cell if it is not periodic.
        """
        pbc = np.asarray(atoms.pbc)

        if not pbc.any():
            return Cell()

        if not pbc.all():
            self.logger.error(
                (
                    f"The ASE trajectory {self.filename} is not periodic in "
                    f"all dimensions ({pbc.tolist()}), which is not "
                    "supported."
                ),
                exception=AseTrajReaderError,
            )

        return _cell_from_box_vectors(_to_float64(atoms.cell[:]))



def _import_ase() -> Any:
    """
    Imports ase, which is an optional dependency.
    """
    try:
        # pylint: disable=import-outside-toplevel
        import ase.io
        import ase.units
    except ModuleNotFoundError:
        AseTrajReader.logger.error(
            (
                "Reading ASE trajectory files requires ase. Install it with "
                "'pip install PQAnalysis[ase]'."
            ),
            exception=AseTrajReaderError,
        )

    return ase



def _read_topology(atoms: Any) -> Topology:
    """
    Builds the topology from the atomic numbers.
    """
    return Topology(
        atoms=[
            _atom_from_name_and_species(symbol, int(number)) for symbol, number
            in zip(atoms.get_chemical_symbols(), atoms.numbers)
        ]
    )



def _to_float64(values: Any) -> np.ndarray:
    """
    Converts ASE values to float64.
    """
    return np.asarray(values, dtype=np.float64)
