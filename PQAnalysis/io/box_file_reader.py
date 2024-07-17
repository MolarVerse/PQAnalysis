"""
A module to read lattice parameter data from a file or from trajectory data.
"""


import logging
import numpy as np

from PQAnalysis.utils import instance_function_count_decorator
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking, runtime_type_checking_setter
from PQAnalysis.traj import Trajectory
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.exceptions import PQFileNotFoundError
from .exceptions import BoxReaderError
from PQAnalysis.physical_data.box import Box
from .base import BaseReader


class BoxFileReader(BaseReader):
    """
    A class to read lattice parameter data from a file or from trajectory data.
    """

    def __init__(
            self,
            filename: str | None = None,
            trajectory: Trajectory | None = None,
            unit: str | None = None,
            engine_format: MDEngineFormat | str = MDEngineFormat.PQ
    ):
        """
        Initialize the BoxReader object. 
        Either a filename or a trajectory must be provided
        depending which engine was used.
        The lattice parameter data are stored in a box object.

        For more information about the box object, see
        :py:class:`PQAnalysis.physical_data.box.Box`.

        Parameters
        ----------
        filename : str, optional
            The file to read the lattice parameter data from.
        trajectory : object, optional
            The trajectory data to read the lattice parameter data from.
        unit : str, optional
            The unit of the lattice parameters.
        engine_format : object or str, optional
        """

        self.engine_format = engine_format

        if filename is MDEngineFormat.PQ and filename is not None:
            super().__init__(filename)

        elif trajectory is not None:
            self.trajectory = trajectory
        else:
            raise PQFileNotFoundError(
                "Either a filename or a trajectory must be provided depending on the engine format.")
        if unit is not None:
            self.unit = unit

    def read(self):
        """
        Read the lattice parameter data from the file or trajectory data.

        Returns
        -------
        object
            The lattice parameter data.
        """

        if self.engine_format is MDEngineFormat.PQ:
            return self._read_from_file()
        else:
            return self._read_from_trajectory()

    def _read_from_file(self):
        """
        Read the lattice parameter data from a file.

        Returns
        -------
        Box
            The lattice parameter data.
        """
        with open(self.filename, 'r', encoding='utf-8') as file:
            a = []
            b = []
            c = []
            alpha = []
            beta = []
            gamma = []

            for line in file:
                if line.startswith("#"):
                    continue
                line = line.split()
                a.append(float(line[0]))
                b.append(float(line[1]))
                c.append(float(line[2]))
                alpha.append(float(line[3]))
                beta.append(float(line[4]))
                gamma.append(float(line[5]))
            return Box(np.array(a), np.array(b), np.array(c), np.array(alpha), np.array(beta), np.array(gamma), self.unit)

    def _read_from_trajectory(self):
        """
        Read the lattice parameter data from trajectory data.

        Returns
        -------
        Box
            The lattice parameter data.
        """
        self.__check_pbc__(self.trajectory)

        a = []
        b = []
        c = []
        alpha = []
        beta = []
        gamma = []

        for i, frame in enumerate(self.trajectory):
            cell = frame.cell
            a.append(cell.a)
            b.append(cell.b)
            c.append(cell.c)
            alpha.append(cell.alpha)
            beta.append(cell.beta)
            gamma.append(cell.gamma)

        return Box(np.array(a), np.array(b), np.array(c), np.array(alpha), np.array(beta), np.array(gamma), self.unit)

    def __check_pbc__(self, traj: Trajectory) -> None:
        """
        Checks if the cell of the trajectory is not None.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to check.

        Raises
        ------
        ValueError
            If the cell of a frame of the trajectory is None.
        """

        if not traj.check_pbc():
            self.logger.error(
                "At least on cell of the trajectory is None. Cannot write box file.",
                exception=BoxReaderError
            )
