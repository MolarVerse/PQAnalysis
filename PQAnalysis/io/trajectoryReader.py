import numpy as np

from PQAnalysis.traj.frame import Frame
from PQAnalysis.traj.trajectory import Trajectory
from PQAnalysis.pbc.cell import Cell


class TrajectoryReader:
    def __init__(self, filename):
        self.filename = filename
        self.frames = []

    def read(self):
        frame_reader = FrameReader()
        with open(self.filename, 'r') as f:
            frame_string = ''
            for line in f:
                if line.strip() == '':
                    frame_string += line
                elif line.split()[0].isdigit():
                    if frame_string != '':
                        self.frames.append(frame_reader.read(frame_string))
                    frame_string = line
                else:
                    frame_string += line
            self.frames.append(frame_reader.read(frame_string))
        return Trajectory(self.frames)


class FrameReader:
    def read(self, frame_string):

        splitted_frame_string = frame_string.split('\n')
        header_line = splitted_frame_string[0].split()

        n_atoms, cell = self.read_header_line(header_line)

        xyz, atoms = self.read_xyz(splitted_frame_string, n_atoms)

        return Frame(n_atoms, xyz, np.array(atoms), cell)

    def read_header_line(self, header_line):
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

        return n_atoms, cell

    def read_xyz(self, splitted_frame_string, n_atoms):

        xyz = np.zeros((n_atoms, 3))
        atoms = []
        for i in range(n_atoms):
            line = splitted_frame_string[2+i]
            xyz[i] = np.array([float(x) for x in line.split()[1:4]])
            atoms.append(line.split()[0])

        return xyz, atoms
