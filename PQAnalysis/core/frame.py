import numpy as np

def read_frame(filename):
    '''
    Read a frame from a file.
    '''
    file = open(filename, 'r')
    splitted_line = file.readline().split()
    if len(splitted_line) == 4:
        n_atoms, a, b, c = int(splitted_line[0]), map(float, splitted_line[1:4])
        cell = cell.Cell(a, b, c)
    elif len(splitted_line) == 7:
        n_atoms, a, b, c, alpha, beta, gamma = int(splitted_line[0]), map(float, splitted_line[1:7])
        cell = cell.Cell(a, b, c, alpha, beta, gamma)
    elif len(splitted_line) == 1:
        n_atoms= int(splitted_line[0])
        cell = None
    else:
        raise ValueError('Invalid file format.')
    
    file.readline() # Skip the second line

    xyz = np.zeros((n_atoms, 3))
    atoms = []
    for i in range(n_atoms):
        line = file.readline()
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