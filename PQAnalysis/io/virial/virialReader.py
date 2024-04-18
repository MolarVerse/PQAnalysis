import numpy as np

from PQAnalysis.io import BaseReader


class _BaseReader(BaseReader):
    def __init__(self, filename: str) -> None:
        super().__init__(filename)

    def read(self):
        data = []
        for line in self.file:
            matrix = np.zeros((3, 3))
            for i, row in enumerate(line.split()[1:]):
                matrix[i] = float(row)

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
