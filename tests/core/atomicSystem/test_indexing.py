import numpy as np
import pytest

from PQAnalysis.core import Atom, AtomicSystem
from PQAnalysis.core.exceptions import AtomicSystemEmptySelectionWarning as EmptySelectionWarning


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

    def test_indices_from_atoms(self):
        indices = self.system.indices_from_atoms(['H', 'O'])

        assert np.allclose(indices, [0, 1])

        with pytest.warns(EmptySelectionWarning):
            self.system.indices_from_atoms(['C'])

        indices = self.system.indices_from_atoms(
            [self.atoms[0]], use_full_atom_info=True)

        assert np.allclose(indices, [0])

        indices = self.system.indices_from_atoms(
            [self.atoms[0]])

        assert np.allclose(indices, [0, 2])

        with pytest.warns(EmptySelectionWarning):
            self.system.indices_from_atoms([Atom('C')])

        indices = self.system.indices_from_atoms([Atom(1), Atom(8)])

        assert np.allclose(indices, [0, 1, 2, 3])

        with pytest.warns(EmptySelectionWarning):
            self.system.indices_from_atoms([Atom(6)])

        indices = self.system.indices_from_atoms(None)

        assert np.allclose(indices, [0, 1, 2, 3])
