"""
A module containing the BoxReader class and its associated methods.
"""

import logging

import numpy as np

from PQAnalysis.io.base import BaseReader
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

from .exceptions import BoxReaderError



class BoxReader(BaseReader):

    """
    A class for reading data-style box files.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def read(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Reads a data-style box file.

        Returns
        -------
        tuple[np.ndarray, np.ndarray, np.ndarray]
            The steps, box lengths and box angles from the box file.
        """
        steps = []
        box_lengths = []
        box_angles = []

        with open(self.filename, "r", encoding="utf-8") as file:
            for line_number, line in enumerate(file, start=1):
                line = line.strip()

                if line == "" or line.startswith("#"):
                    continue

                step, lengths, angles = self._parse_line(line, line_number)

                steps.append(step)
                box_lengths.append(lengths)
                box_angles.append(angles)

        return (
            np.array(steps, dtype=int),
            np.array(box_lengths, dtype=float),
            np.array(box_angles, dtype=float),
        )

    @classmethod
    def _parse_line(
        cls,
        line: str,
        line_number: int,
    ) -> tuple[int, list[float], list[float]]:
        """
        Parses a data-style box file line.
        """
        values = line.split()

        if len(values) != 7:
            cls.logger.error(
                (
                    f"Invalid number of columns in box file line {line_number}. "
                    "Expected 7 columns."
                ),
                exception=BoxReaderError
            )

        try:
            step = int(values[0])
            box_lengths = [float(value) for value in values[1:4]]
            box_angles = [float(value) for value in values[4:]]
        except ValueError:
            cls.logger.error(
                (
                    f"Invalid numeric value in box file line {line_number}: "
                    f"{line}"
                ),
                exception=BoxReaderError
            )

        return step, box_lengths, box_angles


@runtime_type_checking
def read_box(filename: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Reads a data-style box file.

    Parameters
    ----------
    filename : str
        The box file to read.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        The steps, box lengths and box angles from the box file.
    """
    return BoxReader(filename).read()
