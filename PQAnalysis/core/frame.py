import numpy as np
from .cell import Cell


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

    return Frame(n_atoms, xyz, atoms, cell)


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
