"""
A module containing the BoxWriter class and its associated methods.
"""

import logging

from PQAnalysis.traj import Trajectory
from PQAnalysis.utils import instance_function_count_decorator
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking, runtime_type_checking_setter

from .base import BaseWriter
from .formats import BoxFileFormat, FileWritingMode
from .exceptions import BoxWriterError



class BoxWriter(BaseWriter):

    """
    A class for writing a trajectory to a box file.
    Inherits from BaseWriter

    It can write a trajectory to a box file in either a data file
    format or a VMD file format. For more information see
    :py:class:`~PQAnalysis.io.formats.BoxFileFormat`.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(
        self,
        filename: str | None = None,
        output_format: str | BoxFileFormat = 'data',
        mode: str | FileWritingMode = 'w'
    ) -> None:
        """
        Parameters
        ----------
        filename : str, optional
            The name of the file to write to.
            If None, the output is printed to stdout.
        output_format : str | BoxFileFormat, optional
            The format of the file.
            The default is 'data' i.e. BoxFileFormat.DATA.
        mode : str | FileWritingMode, optional
            The mode of the file. Either 'w' for write,
            'a' for append or 'o' for overwrite. The default is 'w'.

        Raises
        ------
        ValueError
            If the given format is not in :py:class:`~PQAnalysis.io.formats.BoxFileFormat`.
        """

        super().__init__(filename, FileWritingMode(mode))
        self.output_format = BoxFileFormat(output_format)

    @runtime_type_checking
    def write(self, traj: Trajectory, reset_counter: bool = True) -> None:
        """
        Wrapper to write the given trajectory to the file.
        Depending on the format, either write_vmd or write_box_file is called.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to write.
        reset_counter : bool, optional
            If True, the function execution counter of write_box_file
            is reset to 0, otherwise it is not reset.
        """

        self.open()
        if self.output_format == BoxFileFormat.VMD:
            self.write_vmd(traj)
        else:
            self.write_box_file(traj, reset_counter=reset_counter)

        self.close()

    @runtime_type_checking
    def write_vmd(self, traj: Trajectory) -> None:
        """
        Writes the given trajectory to the file in VMD format.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to write.
        """
        self.__check_pbc__(traj)

        for frame in traj:
            cell = frame.cell

            print("8", file=self.file)
            print(
                (
                f"Box   "
                f"{cell.x} {cell.y} {cell.z}    "
                f"{cell.alpha} {cell.beta} {cell.gamma}"
                ),
                file=self.file
            )

            edges = cell.bounding_edges

            for edge in edges:
                print(f"X   {edge[0]} {edge[1]} {edge[2]}", file=self.file)

    @instance_function_count_decorator
    @runtime_type_checking
    def write_box_file(
        self,
        traj: Trajectory,
        reset_counter: bool = True  # pylint: disable=unused-argument # is needed for the decorator
    ) -> None:
        """
        Writes the given trajectory to the file in data file format.

        The @count_decorator is used to count the number of frames written to the file. The default
        way is that the counter is reset to 0 after each call of the function. This can be changed
        by setting the reset_counter parameter to False.


        Parameters
        ----------
        traj : Trajectory
            The trajectory to write.
        reset_counter : bool, optional
            If True, the counter is reset to 0, otherwise it is not reset.

        Raises
        ------
        BoxWriterError
            If the cell of a frame of the trajectory is None.
        """
        self.__check_pbc__(traj)

        counter = self.counter[BoxWriter.write_box_file.__name__]  # pylint: disable=no-member # is added via decorator
        counter = len(traj) * (counter - 1)

        for i, frame in enumerate(traj):
            cell = frame.cell
            print(
                (
                f"{counter + i+1} "
                f"{cell.x} {cell.y} {cell.z} "
                f"{cell.alpha} {cell.beta} {cell.gamma}"
                ),
                file=self.file
            )

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
                exception=BoxWriterError
            )

    @property
    def output_format(self) -> BoxFileFormat:
        """BoxFileFormat: The format of the file."""
        return self._output_format

    @output_format.setter
    @runtime_type_checking_setter
    def output_format(self, output_format: str | BoxFileFormat) -> None:
        self._output_format = BoxFileFormat(output_format)
