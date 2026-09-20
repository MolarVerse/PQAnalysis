"""
A module containing classes to read virial and stress files.
"""

import logging

import numpy as np

from PQAnalysis.io import BaseReader
from PQAnalysis.io.exceptions import VirialFileReaderError
from PQAnalysis.types import Np3x3NumberArray
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__



class _BaseReader(BaseReader):

    """
    A base class for reading virial and stress files.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

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

            for line_number, line in enumerate(file, start=1):
                line = line.strip()

                if line == "" or line.startswith("#"):
                    continue

                line_elements = line.split()

                if len(line_elements) != 10:
                    self.logger.error(
                        (
                            f"Invalid number of columns in file "
                            f"{self.filename} line {line_number}. "
                            "Expected 10 columns."
                        ),
                        exception=VirialFileReaderError
                    )

                matrix = np.zeros((3, 3))
                line_elements = line_elements[1:]
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
