import numpy as np

from PQAnalysis.atom.element import Element


class Molecule:
    def __init__(self, atoms, xyz=None, name=None):
        self.atoms = [Element(atom) for atom in atoms]
        self.xyz = xyz

        if name is None:
            self.name = ''.join(atoms)
        else:
            self.name = name

    @property
    def atom_masses(self):
        return [atom.mass for atom in self.atoms]

    @property
    def mass(self):
        return sum(self.atom_masses)

    def com(self, cell=None):
        if self.xyz is None:
            raise ValueError('xyz must be provided when computing com.')

        if cell is not None:
            xyz_imaged = cell.image(self.xyz - self.xyz[0]) + self.xyz[0]
        else:
            xyz_imaged = self.xyz

        com = sum(np.array([self.atom_masses]).T *
                  xyz_imaged) / sum(self.atom_masses)

        if cell is not None:
            com = cell.image(com)

        return com
