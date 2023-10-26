import numpy as np
import sys

from PQAnalysis.selection.selection import Selection
from PQAnalysis.pbc.cell import Cell
from PQAnalysis.atom.molecule import Molecule


def read_frame(frame_string):
    '''
    Read a frame from a file.
    '''


class Frame:
    '''
    A class to represent a frame of a trajectory.
    '''

    def __init__(self, n_atoms, xyz, atoms, cell=None):
        self.n_atoms = n_atoms
        self.xyz = xyz
        self.atoms = atoms
        self.cell = cell

        if cell is not None:
            self.PBC = True
        else:
            self.PBC = False

    def __getindex__(self, index):
        if isinstance(index, Selection):
            frame = Frame(
                index.n_atoms, self.xyz[index.selection], index.atoms[index.selection], cell=self.cell)
        else:
            frame = Frame(1, self.xyz[index], np.array(
                [self.atoms[index]]), cell=self.cell)

        if frame.n_atoms == 0:
            raise ValueError('Selection is empty.')

        return frame

    def compute_com(self, group=None):

        if group is None:
            group = self.n_atoms
        elif self.n_atoms % group != 0:
            raise ValueError(
                'Number of atoms in selection is not a multiple of group.')

        com = np.zeros((self.n_atoms // group, 3))
        molecule_names = np.zeros(self.n_atoms // group, dtype=object)

        for i in range(0, self.n_atoms, group):
            molecule = Molecule(self.xyz[i:i+group],
                                self.atoms[i:i+group])

            com[i] = molecule.compute_com()

            molecule_names[i] = molecule.name

        return Frame(self.n_atoms // group, com, molecule_names, cell=self.cell)
