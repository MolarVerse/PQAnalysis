import sys

from PQAnalysis.io.base import BaseWriter


def write_trajectory(traj, filename=None, format=None):
    '''
    Wrapper for TrajectoryWriter to write a trajectory to a file.
    '''

    writer = TrajectoryWriter(filename, format)
    writer.write(traj)


class TrajectoryWriter(BaseWriter):

    formats = [None, 'qmcfc']

    def __init__(self, filename=None, format=None, mode='w'):
        super().__init__(filename, mode)
        if format not in self.formats:
            raise ValueError(
                'Invalid format. Has to be either \'qmcfc\' or \'None\'.')

        self.format = format

    def write(self, trajectory):
        self.open()
        for frame in trajectory:
            self._write_header(frame.n_atoms, frame.cell)
            self._write_coordinates(frame.xyz, frame.atoms)
        self.close()

    def _write_header(self, n_atoms, cell):
        if cell is not None:
            print(
                f"{n_atoms} {cell.x} {cell.y} {cell.z} {cell.alpha} {cell.beta} {cell.gamma}\n", file=self.file)
        else:
            print(f"{n_atoms}\n", file=self.file)

    def _write_coordinates(self, xyz, atoms):

        if self.format == "qmcfc":
            print("X   0.0 0.0 0.0", file=self.file)

        for i in range(len(atoms)):
            print(
                f"{atoms[i]} {xyz[i][0]} {xyz[i][1]} {xyz[i][2]}", file=self.file)
