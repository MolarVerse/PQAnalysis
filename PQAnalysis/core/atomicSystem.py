import numpy as np
import numpy.typing as npt


from typeguard import typechecked
from beartype import beartype
from beartype.typing import List
from dataclasses import dataclass

from .atom import Atom
from .cell import Cell
from ..utils.mytypes import Numpy2DFloatArray, Numpy1DFloatArray


def check_atoms_pos(func):
    def wrapper(*args, **kwargs):
        if args[0].pos.shape[0] != args[0].n_atoms:
            raise ValueError(
                "AtomicSystem contains a different number of atoms to positions.")
        return func(*args, **kwargs)
    return wrapper


def check_atoms_has_mass(func):
    def wrapper(*args, **kwargs):
        if not all([atom.mass is not None for atom in args[0].atoms]):
            raise ValueError(
                "AtomicSystem contains atoms without mass information.")
        return func(*args, **kwargs)
    return wrapper


@beartype
class AtomicSystem:
    def __init__(self,
                 atoms: List[Atom] = None,
                 pos: Numpy2DFloatArray = np.zeros((0, 3)),
                 vel: Numpy2DFloatArray = np.zeros((0, 3)),
                 forces: Numpy2DFloatArray = np.zeros((0, 3)),
                 charges: Numpy1DFloatArray = np.zeros(0),
                 cell: Cell | None = None
                 ):

        if atoms is None:
            atoms = []

        self._atoms = atoms
        self._pos = pos
        self._vel = vel
        self._forces = forces
        self._charges = charges
        self._cell = cell

    @property
    def atoms(self) -> List[Atom]:
        return self._atoms

    @property
    def pos(self) -> Numpy2DFloatArray:
        return self._pos

    @property
    def vel(self) -> Numpy2DFloatArray:
        return self._vel

    @property
    def forces(self) -> Numpy2DFloatArray:
        return self._forces

    @property
    def charges(self) -> Numpy1DFloatArray:
        return self._charges

    @property
    def cell(self) -> Cell | None:
        return self._cell

    @property
    def PBC(self) -> bool:
        return self._cell is not None

    @property
    def n_atoms(self) -> int:
        return len(self._atoms)

    @property
    @check_atoms_has_mass
    def atomic_masses(self) -> Numpy1DFloatArray:
        return np.array([atom.mass for atom in self._atoms])

    @property
    def mass(self) -> float:
        return np.sum(self.atomic_masses)

    @property
    @check_atoms_pos
    @check_atoms_has_mass
    def center_of_mass(self) -> Numpy1DFloatArray:

        if self.n_atoms == 0:
            return np.zeros(3)

        if self.cell is not None:
            pos = self.cell.image(self.pos - self.pos[0]) + self.pos[0]
        else:
            pos = self.pos

        pos = np.sum(
            pos * self.atomic_masses[:, None], axis=0) / self.mass

        if self.cell is not None:
            pos = self.cell.image(pos)

        return pos

    @property
    def combined_name(self) -> str:
        return ''.join([atom.name for atom in self.atoms])
