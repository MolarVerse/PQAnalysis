from PQAnalysis.atom.element import Element


class Molecule:
    def __init__(self, atoms, xyz=None, name=None):
        self.atoms = [Element(atom) for atom in atoms]
        self.xyz = xyz

        if name is None:
            self.name = ''.join(atoms)
        else:
            self.name = name

    def atom_masses(self):
        return [atom.mass for atom in self.atoms]

    def compute_com(self):
        _atom_masses = self.atom_masses()

        if self.xyz is None:
            raise ValueError('xyz must be provided when computing com.')

        return sum(_atom_masses * self.xyz) / sum(_atom_masses)
