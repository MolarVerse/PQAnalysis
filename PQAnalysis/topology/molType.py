from numbers import Real
from beartype.typing import List

from ..types import Np1DIntArray, Np1DNumberArray
from ..core.atom import Atom


class MolType:
    def __init__(self,
                 name: str,
                 total_charge: Real,
                 elements: List[Atom],
                 atom_types: Np1DIntArray,
                 partial_charges: Np1DNumberArray,
                 ) -> None:

        if not (len(elements) == len(atom_types) == len(partial_charges)):
            raise ValueError(
                "The number of elements, atom_types and partial_charges must be the same.")

        self.name = name
        self.total_charge = total_charge

        # set here the internal variables to avoid setters
        # (which would check the length of the arrays)
        self._elements = elements
        self._atom_types = atom_types
        self._partial_charges = partial_charges

    @property
    def n_atoms(self) -> int:
        return len(self._elements)

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        self._name = name

    @property
    def total_charge(self) -> Real:
        return self._total_charge

    @total_charge.setter
    def total_charge(self, total_charge: Real) -> None:
        self._total_charge = total_charge

    @property
    def elements(self) -> List[Atom]:
        return self._elements

    @elements.setter
    def elements(self, elements: List[Atom]) -> None:
        if len(elements) != self.n_atoms:
            raise ValueError(
                "The number of elements must be the same as the number of atoms.")
        self._elements = elements

    @property
    def atom_types(self) -> Np1DIntArray:
        return self._atom_types

    @atom_types.setter
    def atom_types(self, atom_types: Np1DIntArray) -> None:
        if len(atom_types) != self.n_atoms:
            raise ValueError(
                "The number of atom_types must be the same as the number of atoms.")
        self._atom_types = atom_types

    @property
    def partial_charges(self) -> Np1DNumberArray:
        if self._partial_charges is None:
            raise ValueError("The partial charges have not been set yet.")
        return self._partial_charges
