import numpy as np
import sys

from PQAnalysis.selection.selection import Selection
from PQAnalysis.pbc.cell import Cell
from PQAnalysis.atom.molecule import Molecule


def read_frame(frame_string):
    '''
    Read a frame from a file.
    '''

    splitted_frame_string = frame_string.split('\n')

    header_line = splitted_frame_string[0].split()
    if len(header_line) == 4:
        n_atoms = int(header_line[0])
        a, b, c = map(float, header_line[1:4])
        cell = Cell(a, b, c)
    elif len(header_line) == 7:
        n_atoms = int(header_line[0])
        a, b, c, alpha, beta, gamma = map(float, header_line[1:7])
        cell = Cell(a, b, c, alpha, beta, gamma)
    elif len(header_line) == 1:
        n_atoms = int(header_line[0])
        cell = None
    else:
        raise ValueError('Invalid file format in header line of Frame.')

    xyz = np.zeros((n_atoms, 3))
    atoms = []
    for i in range(n_atoms):
        line = splitted_frame_string[2+i]
        xyz[i] = np.array([float(x) for x in line.split()[1:4]])
        atoms.append(line.split()[0])

    return Frame(n_atoms, xyz, np.array(atoms), cell)


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

    def write_xyz_header(self, filename):
        print(
            f"{self.n_atoms} {self.cell.x} {self.cell.y} {self.cell.z} {self.cell.alpha} {self.cell.beta} {self.cell.gamma}", file=filename)

    def write_coordinates(self, filename):
        for i in range(self.n_atoms):
            print(
                f"{self.atoms[i]} {self.xyz[i][0]} {self.xyz[i][1]} {self.xyz[i][2]}", file=filename)

    def write(self, filename=sys.stdout, format=None):
        '''
        Write a frame in to a file.
        '''

        self.write_xyz_header(filename)
        print('', file=filename)
        if format == "qmcfc":
            print("X   0.0 0.0 0.0", file=filename)
        self.write_coordinates(filename)

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
