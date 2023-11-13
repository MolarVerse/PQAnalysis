import numpy as np

from PQAnalysis.core.atom import Atom
from PQAnalysis.core.atomicSystem import AtomicSystem


class Test_IndexingMixin:
    atoms = [Atom('H', 1), Atom('O', 8), Atom('H2', 1), Atom('O2', 8)]
    system = AtomicSystem(atoms=atoms)

    def test__indices_by_atom_type_names(self):
        indices = self.system._indices_by_atom_type_names(['H', 'O'])

        assert np.allclose(indices, [0, 1])

    def test__indices_by_atom(self):
        indices = self.system._indices_by_atom([self.atoms[0]])

        assert np.allclose(indices, [0])

    def test__indices_by_element_types(self):
        indices = self.system._indices_by_element_types([Atom(1), Atom(8)])

        assert np.allclose(indices, [0, 1, 2, 3])
