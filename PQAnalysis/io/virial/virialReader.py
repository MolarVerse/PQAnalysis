import numpy as np

from PQAnalysis.io import BaseReader


class _BaseReader(BaseReader):
    def __init__(self, filename: str) -> None:
        super().__init__(filename)

    def read(self):

        with open(self.filename, 'r') as file:

            data = []

            for line in file.readlines():
                matrix = np.zeros((3, 3))
                line_elements = line.split()[1:]
                for i in range(3):
                    for j in range(3):
                        matrix[i, j] = float(line_elements[i * 3 + j])

                data.append(matrix)

            return data


class VirialFileReader(_BaseReader):
    """
    A class to read virial files.
    """
    pass


class StressFileReader(_BaseReader):
    """
    A class to read stress files.
    """
    pass
