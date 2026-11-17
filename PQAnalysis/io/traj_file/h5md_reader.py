"""
A module containing a reader for H5MD trajectory files.

H5MD files are HDF5 files, see https://www.nongnu.org/h5md/. Reading
them needs the optional dependency ``h5py``
(``pip install PQAnalysis[h5md]``).
"""

import logging

import numpy as np

from beartype.typing import Any, Generator

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom, Cell
from PQAnalysis.core.atom.element import atomicNumbers, atomicNumbersReverse
from PQAnalysis.topology import Topology
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis import __package_name__

from .exceptions import H5MDReaderError

#: Factors from H5MD length units to Angstrom.
LENGTH_UNITS = {
    "Angstrom": 1.0,
    "nm": 10.0,
}

#: Factors from H5MD velocity units to Angstrom s-1.
VELOCITY_UNITS = {
    "Angstrom s-1": 1.0,
    "Angstrom ps-1": 1.0e12,
    "Angstrom fs-1": 1.0e15,
    "nm ps-1": 1.0e13,
}

#: Factors from H5MD force units to kcal mol-1 Angstrom-1.
FORCE_UNITS = {"kcal mol-1 Angstrom-1": 1.0}

#: Factors from H5MD charge units to e.
CHARGE_UNITS = {"e": 1.0}

#: Factors from H5MD energy units to kcal mol-1.
ENERGY_UNITS = {"kcal mol-1": 1.0}

# H5MD particle elements and their AtomicSystem keywords
_PARTICLE_KEYWORDS = {
    "position": "pos",
    "velocity": "vel",
    "force": "forces",
    "charge": "charges",
}

_PARTICLE_UNITS = {
    "position": LENGTH_UNITS,
    "velocity": VELOCITY_UNITS,
    "force": FORCE_UNITS,
    "charge": CHARGE_UNITS,
}

_ENERGY_ELEMENT = "observables/potential_energy"



class H5MDReader:

    """
    A class for reading the frames of a single H5MD trajectory file.

    All present time-dependent elements must share the same steps. Values
    are converted to the PQ units used in PQAnalysis: Angstrom, Angstrom
    s-1, kcal mol-1 Angstrom-1, e and kcal mol-1. A unit that is missing
    or not known raises an :py:class:`H5MDReaderError`.

    Generally, this class should not be used directly. Instead, use the
    :py:class:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader`
    with ``traj_format="h5md"`` or a ``.h5md``/``.h5`` file.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            The name of the H5MD file.

        Raises
        ------
        H5MDReaderError
            If h5py is not installed.
        """
        self.filename = filename
        self._h5py = _import_h5py()

    @property
    def number_of_frames(self) -> int:
        """int: The number of frames in the file."""
        with self._open() as file:
            return len(self._steps(file))

    def cell_generator(self) -> Generator[Cell]:
        """
        A generator that yields the cell of each frame.

        Yields
        ------
        Cell
            The cell of the frame, a vacuum cell if the file has no box.
        """
        with self._open() as file:
            n_frames = len(self._steps(file))
            for index in range(n_frames):
                yield self._read_cell(file, index)

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
            built from ``/parameters/atom_names`` and the species.

        Yields
        ------
        AtomicSystem
            The frames with positions in Angstrom, velocities in
            Angstrom s-1, forces in kcal mol-1 Angstrom-1, charges in e
            and the potential energy in kcal mol-1.

        Raises
        ------
        H5MDReaderError
            If the file is not a valid H5MD trajectory or uses a unit that
            is not supported.
        """
        with self._open() as file:
            steps = self._steps(file)

            if stop is None:
                stop = len(steps)

            if topology is None:
                topology = self._read_topology(file)

            particles = self._particle_group(file)
            elements = {}
            for name, keyword in _PARTICLE_KEYWORDS.items():
                if name in particles:
                    dataset = particles[name]["value"]
                    elements[keyword] = (
                        dataset,
                        self._unit_factor(dataset, _PARTICLE_UNITS[name]),
                    )

            energy = None
            if _ENERGY_ELEMENT in file:
                energy_group = file[_ENERGY_ELEMENT]
                energy = (
                    energy_group["value"],
                    self._unit_factor(energy_group["value"], ENERGY_UNITS),
                )

            for index in range(start, stop):
                values = {
                    keyword: _to_float64(dataset[index]) * factor
                    for keyword, (dataset, factor) in elements.items()
                }

                if energy is not None:
                    values["energy"] = float(energy[0][index]) * energy[1]

                yield AtomicSystem(
                    topology=topology,
                    cell=self._read_cell(file, index),
                    **values,
                )

    def _open(self) -> Any:
        """
        Opens the H5MD file for reading.
        """
        try:
            return self._h5py.File(self.filename, "r")
        except OSError:
            self.logger.error(
                f"Could not open the H5MD file {self.filename}.",
                exception=H5MDReaderError,
            )

        return None  # pragma: no cover

    def _particle_group(self, file: Any) -> Any:
        """
        Returns the particle group, ``/particles/all`` or the only one.
        """
        groups = list(file.get("particles", {}).keys())

        if "all" in groups:
            return file["particles/all"]

        if len(groups) == 1:
            return file["particles"][groups[0]]

        self.logger.error(
            (
                f"The H5MD file {self.filename} needs a single particle "
                "group or a group named 'all' in /particles."
            ),
            exception=H5MDReaderError,
        )

        return None  # pragma: no cover

    def _steps(self, file: Any) -> np.ndarray:
        """
        Returns the steps shared by all time-dependent elements.
        """
        particles = self._particle_group(file)

        groups = [
            particles[name] for name in _PARTICLE_KEYWORDS if name in particles
        ]

        if not groups:
            self.logger.error(
                (
                    f"The H5MD file {self.filename} has no position, "
                    "velocity, force or charge data."
                ),
                exception=H5MDReaderError,
            )

        edges = particles.get("box/edges")
        if isinstance(edges, self._h5py.Group):
            groups.append(edges)

        if _ENERGY_ELEMENT in file:
            groups.append(file[_ENERGY_ELEMENT])

        steps = np.asarray(groups[0]["step"])
        for group in groups[1:]:
            if not np.array_equal(np.asarray(group["step"]), steps):
                self.logger.error(
                    (
                        f"The H5MD element {group.name} in {self.filename} "
                        "has other steps than the positions. Elements "
                        "with different steps are not supported."
                    ),
                    exception=H5MDReaderError,
                )

        return steps

    def _unit_factor(self, dataset: Any, units: dict[str, float]) -> float:
        """
        Returns the factor converting the dataset values to PQ units.
        """
        unit = dataset.attrs.get("unit")
        if isinstance(unit, np.ndarray) and unit.size == 1:
            unit = unit.item()
        if isinstance(unit, bytes):
            unit = unit.decode()
        if isinstance(unit, str):
            unit = " ".join(unit.split())

        if not isinstance(unit, str) or unit not in units:
            self.logger.error(
                (
                    f"Unit {unit!r} of the H5MD dataset {dataset.name} in "
                    f"{self.filename} is not supported. Supported units "
                    f"are: {', '.join(units)}."
                ),
                exception=H5MDReaderError,
            )

        return units[unit]

    def _read_cell(self, file: Any, index: int) -> Cell:
        """
        Reads the cell of a frame from ``box/edges``.
        """
        particles = self._particle_group(file)

        if "box/edges" not in particles:
            return Cell()

        boundary = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in particles["box"].attrs.get("boundary", [])
        ]
        if any(value != "periodic" for value in boundary):
            self.logger.error(
                (
                    f"The H5MD box in {self.filename} is not periodic in "
                    f"all dimensions ({boundary}), which is not supported."
                ),
                exception=H5MDReaderError,
            )

        edges = particles["box/edges"]
        if isinstance(edges, self._h5py.Group):
            edges = edges["value"]
            values = edges[index]
        else:
            values = edges[()]

        edges = _to_float64(values) * self._unit_factor(edges, LENGTH_UNITS)

        if edges.shape == (3, ):
            return Cell(*edges)

        return _cell_from_box_vectors(edges)

    def _read_topology(self, file: Any) -> Topology:
        """
        Builds the topology from the atom names and species.
        """
        particles = self._particle_group(file)

        names = None
        if "parameters/atom_names" in file:
            names = [
                name.decode() if isinstance(name, bytes) else str(name)
                for name in file["parameters/atom_names"][()]
            ]

        species = None
        if "species" in particles:
            if not isinstance(particles["species"], self._h5py.Dataset):
                self.logger.error(
                    (
                        f"Time-dependent species in {self.filename} are "
                        "not supported."
                    ),
                    exception=H5MDReaderError,
                )
            species = [int(value) for value in particles["species"][()]]

        if names is None and species is None:
            self.logger.error(
                (
                    f"The H5MD file {self.filename} has neither "
                    "/parameters/atom_names nor species."
                ),
                exception=H5MDReaderError,
            )

        if names is None:
            names = [
                atomicNumbersReverse.get(value, "X") if value > 0 else "X"
                for value in species
            ]

        if species is None:
            species = [0] * len(names)

        return Topology(
            atoms=[
                _atom_from_name_and_species(name, value)
                for name, value in zip(names, species)
            ]
        )



def _import_h5py() -> Any:
    """
    Imports h5py, which is an optional dependency.
    """
    try:
        import h5py  # pylint: disable=import-outside-toplevel
    except ModuleNotFoundError:
        H5MDReader.logger.error(
            (
                "Reading H5MD files requires h5py. Install it with "
                "'pip install PQAnalysis[h5md]'."
            ),
            exception=H5MDReaderError,
        )

    return h5py



def _to_float64(values: Any) -> np.ndarray:
    """
    Converts H5MD values, often float32, to float64.
    """
    return np.asarray(values, dtype=np.float64)



def _cell_from_box_vectors(box_vectors: np.ndarray) -> Cell:
    """
    Builds a cell from the three box vectors, one per row.
    """
    lengths = np.linalg.norm(box_vectors, axis=1)

    def angle(i: int, j: int) -> float:
        cosine = box_vectors[i] @ box_vectors[j] / (lengths[i] * lengths[j])
        return np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0)))

    return Cell(*lengths, angle(1, 2), angle(0, 2), angle(0, 1))



def _atom_from_name_and_species(name: str, species: int) -> Atom:
    """
    Creates an atom like the xyz reader does, with the species as fallback.
    """
    guessed = atomicNumbers.get(name.lower())

    if guessed is not None and species in (0, guessed):
        return Atom(name, disable_type_checking=True)

    if species in atomicNumbersReverse and species > 0:
        return Atom(name, species, disable_type_checking=True)

    return Atom(name, use_guess_element=False, disable_type_checking=True)
