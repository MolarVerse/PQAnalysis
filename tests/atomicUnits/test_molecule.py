import numpy as np
import pytest

from PQAnalysis.atomicUnits.molecule import Molecule
from PQAnalysis.atomicUnits.element import Atom, Elements
from PQAnalysis.core.cell import Cell
from PQAnalysis.coordinates.coordinates import Coordinates


class TestMolecule:
    elements = Elements(['C', 'H', 'H', 'O'])
    coords = Coordinates(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, -1]]))

    def test__init__(self):
        molecule = Molecule(atoms=self.elements)
        assert molecule.atoms == self.elements
        assert molecule.coordinates == Coordinates()
        assert molecule.name == str("CHHO").lower()

        molecule = Molecule(atoms=self.elements,
                            coordinates=self.coords, name='my_name')

        assert molecule.atoms == self.elements
        assert molecule.coordinates == self.coords
        assert molecule.name == "my_name"

        with pytest.raises(TypeError) as exception:
            Molecule(atoms=Elements(["H"]), name=1)
        assert str(exception.value) == 'name must be a str.'

    def test_atom_masses(self):
        molecule = Molecule(atoms=self.elements)
        assert np.allclose(molecule.atom_masses, np.array(
            [12.0107, 1.00794, 1.00794, 15.9994]))

    def test_mass(self):
        molecule = Molecule(atoms=self.elements)
        assert np.isclose(molecule.mass, 30.02598)

    def test_com(self):
        molecule = Molecule(atoms=self.elements, coordinates=self.coords)
        assert np.allclose(molecule.center_of_mass.xyz, np.array(
            [0.03356893,  0.03356893, -0.53285188]))

        local_coords = self.coords.copy()
        local_coords.cell = Cell(10, 10, 10)
        molecule = Molecule(atoms=self.elements, coordinates=local_coords)
        assert np.allclose(molecule.center_of_mass.xyz, np.array(
            [0.03356893,  0.03356893, -0.53285188]))

        molecule = Molecule(atoms=Elements(['C', 'H', 'H', 'O']))
        with pytest.raises(ValueError) as exception:
            molecule.center_of_mass
        assert str(
            exception.value) == 'coordinates must be provided when computing com.'
