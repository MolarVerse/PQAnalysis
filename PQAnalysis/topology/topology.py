from ..types import Np1DIntArray


class Topology:

    def __init__(self):
        self.mol_types = None

    @property
    def mol_types(self):
        return self._mol_types

    @mol_types.setter
    def mol_types(self, value: Np1DIntArray | None):
        self._mol_types = value
