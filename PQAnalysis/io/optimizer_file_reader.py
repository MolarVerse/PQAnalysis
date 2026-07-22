"""
A module containing the reader for PQ optimizer output files.
"""

import logging

import numpy as np

from PQAnalysis import __package_name__
from PQAnalysis.physical_data import Energy
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.utils.custom_logging import setup_logger

from .base import BaseReader
from .exceptions import OptimizerReaderError

OPTIMIZER_PARAMETER_UNITS = (
    ("SIMULATION-TIME", "step"),
    ("ABS-ENERGY-CHANGE", "kcal/mol"),
    ("REL-ENERGY-CHANGE", "-"),
    ("MAX-FORCE", "kcal/mol/A"),
    ("RMS-FORCE", "kcal/mol/A"),
    ("REL-ENERGY-CONV", "state"),
    ("ABS-ENERGY-CONV", "state"),
    ("MAX-FORCE-CONV", "state"),
    ("RMS-FORCE-CONV", "state"),
    ("REL-ENERGY-LIMIT", "-"),
    ("ABS-ENERGY-LIMIT", "kcal/mol"),
    ("MAX-FORCE-LIMIT", "kcal/mol/A"),
    ("RMS-FORCE-LIMIT", "kcal/mol/A"),
)
OPTIMIZER_COLUMN_COUNT = len(OPTIMIZER_PARAMETER_UNITS)



class OptimizerFileReader(BaseReader):

    """
    A reader for the fixed 13-column ``.opt`` files written by PQ.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(self, filename: str) -> None:
        """
        Parameters
        ----------
        filename : str
            The optimizer output file to read.
        """
        super().__init__(filename)

    @runtime_type_checking
    def read(self) -> Energy:
        """
        Read the optimizer output and return its values with the PQ schema.

        Returns
        -------
        Energy
            Optimizer values arranged as parameters by optimization steps.

        Raises
        ------
        OptimizerReaderError
            If the file is empty or contains an invalid row.
        """
        rows = []

        with open(self.filename, "r", encoding="utf-8") as file:
            for line_number, line in enumerate(file, start=1):
                line = line.strip()

                if line == "" or line.startswith("#"):
                    continue

                rows.append(self._parse_line(line, line_number))

        if not rows:
            self.logger.error(
                (
                    f"Optimizer file {self.filename} does not contain "
                    "optimizer data."
                ),
                exception=OptimizerReaderError
            )

        info = {
            parameter: index
            for index, (parameter, _) in enumerate(OPTIMIZER_PARAMETER_UNITS)
        }
        units = dict(OPTIMIZER_PARAMETER_UNITS)

        return Energy(np.asarray(rows, dtype=float).T, info, units)

    @classmethod
    def _parse_line(cls, line: str, line_number: int) -> list[float]:
        """
        Parse and validate one optimizer output row.
        """
        values = line.split()

        if len(values) != OPTIMIZER_COLUMN_COUNT:
            cls.logger.error(
                (
                    "Invalid number of columns in optimizer file line "
                    f"{line_number}. Expected {OPTIMIZER_COLUMN_COUNT} "
                    "columns."
                ),
                exception=OptimizerReaderError
            )

        parsed_values = []
        try:
            parsed_values = [float(value) for value in values]
        except ValueError:
            cls.logger.error(
                (
                    "Invalid numeric value in optimizer file line "
                    f"{line_number}: {line}"
                ),
                exception=OptimizerReaderError
            )

        return parsed_values



@runtime_type_checking
def read_optimizer_file(filename: str) -> Energy:
    """
    Read a PQ optimizer output file.

    Parameters
    ----------
    filename : str
        The optimizer output file to read.

    Returns
    -------
    Energy
        Optimizer values arranged as parameters by optimization steps.
    """
    return OptimizerFileReader(filename).read()
