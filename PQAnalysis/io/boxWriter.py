from PQAnalysis.io.base import BaseWriter
from PQAnalysis.utils.decorators import count_decorator


def write_box(traj, filename=None, format=None):
    '''
    Wrapper for BoxWriter to write a trajectory to a file.
    '''

    writer = BoxWriter(filename, format)
    writer.write(traj)


class BoxWriter(BaseWriter):
    formats = [None, 'vmd']

    def __init__(self, filename=None, format=None, mode='w'):
        super().__init__(filename, mode)
        if format not in self.formats:
            raise ValueError(
                'Invalid format. Has to be either \'vmd\' or \'None\'.')

        self.format = format

    def write(self, traj):
        self.open()
        if self.format == "vmd":
            self.write_vmd(traj)
        else:
            self.write_box_file(traj)

        self.close()

    def write_vmd(self, traj):
        for frame in traj:
            cell = frame.cell
            print("8", file=self.file)
            print(
                f"Box   {cell.x} {cell.y} {cell.z}    {cell.alpha} {cell.beta} {cell.gamma}", file=self.file)
            edges = cell.bounding_edges
            for edge in edges:
                print(f"X   {edge[0]} {edge[1]} {edge[2]}", file=self.file)

    @count_decorator
    def write_box_file(self, traj):
        for frame in traj:
            cell = frame.cell
            print(
                f"{self.write_box_file.counter} {cell.x} {cell.y} {cell.z} {cell.alpha} {cell.beta} {cell.gamma}")
