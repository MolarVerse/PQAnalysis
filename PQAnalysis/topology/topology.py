from ..types import Np1DIntArray


class Topology:
    @property
    def mol_types(self):
        return self._mol_types

    @mol_types.setter
    def mol_types(self, value: Np1DIntArray):
        self._mol_types = value
