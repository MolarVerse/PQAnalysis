"""
A module to read lattice parameter data from a file or from trajectory data.
"""
import numpy as np
from PQAnalysis.core.cell import Cell, Cells

from PQAnalysis.traj import Trajectory
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.exceptions import PQFileNotFoundError

from .exceptions import BoxReaderError
from .base import BaseReader


class BoxFileReader(BaseReader):

    """
    A class to read lattice parameter data from a file or from trajectory data.
    """

    def __init__(
        self,
        filename: str | None = None,
        trajectory: Trajectory | None = None,
        engine_format: MDEngineFormat | str = MDEngineFormat.PQ
    ):
        """
        Initialize the BoxReader object. 
        Either a filename or a trajectory must be provided
        depending which engine was used.
        The lattice parameter data are stored in a Cells object.

        For more information about the cell object, see
        :py:class:`PQAnalysis.core.cell`.

        Parameters
        ----------
        filename : str, optional
            The file to read the lattice parameter data from.
        trajectory : object, optional
            The trajectory data to read the lattice parameter data from.
        engine_format : object or str, optional
        """

        self.engine_format = engine_format

        if engine_format is MDEngineFormat.PQ and filename is not None:
            super().__init__(filename)
            self.trajectory = None

        elif trajectory is not None:
            self.trajectory = trajectory
            self.filename = None

        else:
            raise PQFileNotFoundError(
                (
                    "Either a filename or a trajectory must be "
                    "provided depending on the engine format."
                )
            )

    def read(self):
        """
        Read the lattice parameter data from the file or trajectory data.

        Returns
        -------
        object
            The lattice parameter data.
        """

        if self.filename is not None:
            return self._read_from_file()
        return self._read_from_trajectory()

    def _read_from_file(self):  # -> list:
        """
        Read the lattice parameter data from a file.

        Returns
        -------
        Cells
            The lattice parameter data.
        """

        with open(self.filename, 'r', encoding='utf-8') as file:
            cell = []
            i = 0
            while True:
                line = file.readline()
                if not line:
                    break
                if line.startswith("#"):
                    continue
                line = [float(l) for l in line.split()]
                if len(line) == 3:
                    a, b, c = line
                    cell.append(Cell(a, b, c))
                elif len(line) == 4:
                    N, a, b, c = line
                    cell.append(Cell(a, b, c))
                elif len(line) == 7:
                    N, a, b, c, alpha, beta, gamma = line
                    cell.append(Cell(a, b, c, alpha, beta, gamma))
                elif len(line) == 6:
                    a, b, c, alpha, beta, gamma = line
                    cell.append(Cell(a, b, c, alpha, beta, gamma))
                else:
                    self.logger.error(
                        (
                            f"Line {i} in file {self.filename} has "
                            f"an invalid number of columns: {len(line)}"
                        ),
                        exception=BoxReaderError
                    )
                i += 1
            print(f"Read {i} lines from file {self.filename}")
            return cell

    def _read_from_trajectory(self):
        """
        Read the lattice parameter data from trajectory data.

        Returns
        -------
        Cells
            The lattice parameter data.
        """
        self.__check_pbc__(self.trajectory)

        return [frame.cell for frame in self.trajectory]

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
