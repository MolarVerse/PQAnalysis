import numpy as np
import pytest

from PQAnalysis.atom.molecule import Molecule
from PQAnalysis.atom.element import Element
from PQAnalysis.pbc.cell import Cell


def test__init__():
    molecule = Molecule(['C', 'H', 'H', 'O'])
    assert molecule.atoms[0] == Element('C')
    assert molecule.atoms[1] == Element('H')
    assert molecule.atoms[2] == Element('H')
    assert molecule.atoms[3] == Element('O')
    assert molecule.xyz is None
    assert molecule.name == str("CHHO").lower()

    molecule = Molecule(['C', 'H', 'H', 'O'], xyz=np.array([[0, 0, 0], [1, 0, 0], [
                        0, 1, 0], [0, 0, 1]]), name='my_name')

    assert molecule.atoms[0] == Element('C')
    assert molecule.atoms[1] == Element('H')
    assert molecule.atoms[2] == Element('H')
    assert molecule.atoms[3] == Element('O')
    assert np.allclose(molecule.xyz, np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]))
    assert molecule.name == "my_name"


def test_atom_masses():
    molecule = Molecule(['C', 'H', 'H', 'O'])
    assert np.allclose(molecule.atom_masses, np.array(
        [12.0107, 1.00794, 1.00794, 15.9994]))


def test_mass():
    molecule = Molecule(['C', 'H', 'H', 'O'])
    assert np.isclose(molecule.mass, 30.02598)


def test_com():
    molecule = Molecule(['C', 'H', 'H', 'O'], xyz=np.array(
        [[-4, 0, 0], [4, 0, 0], [0, 4, 0], [0, 0, -4]]))
    assert np.allclose(molecule.com(), np.array(
        [-1.46576531,  0.13427572, -2.13140753]))

    molecule = Molecule(['C', 'H', 'H', 'O'], xyz=np.array(
        [[-4, 0, 0], [4, 0, 0], [0, 4, 0], [0, 0, -4]]))
    cell = Cell(10, 10, 10)
    assert np.allclose(molecule.com(cell), np.array(
        [-1.80145461,  0.13427572, -2.13140753]))

    molecule = Molecule(['C', 'H', 'H', 'O'])
    with pytest.raises(ValueError) as exception:
        molecule.com()

    assert str(exception.value) == 'xyz must be provided when computing com.'
