"""
A module containing classes to read virial and stress files.
"""

import numpy as np

from PQAnalysis.io import BaseReader
from PQAnalysis.types import Np3x3NumberArray
from PQAnalysis.type_checking import runtime_type_checking



class _BaseReader(BaseReader):

    """
    A base class for reading virial and stress files.
    """

    @runtime_type_checking
    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            The filename of the file to read.
        """
        super().__init__(filename)

    def read(self) -> list[Np3x3NumberArray]:
        """
        Read the file.

        The file is read line by line. Each line is split into its elements and
        the elements are used to create a 3x3 matrix. One line has the following format:

        ```
        _ xx xy xz yx yy yz zx zy zz
        ```

        Returns
        -------
        list[Np3x3NumberArray]
            The data read from the file.
        """
        with open(self.filename, 'r', encoding='utf-8') as file:

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



class StressFileReader(_BaseReader):

    """
    A class to read stress files.
    """
