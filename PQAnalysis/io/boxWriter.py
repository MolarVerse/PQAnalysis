"""
A module containing the BoxWriter class and its associated methods.

...

Classes
-------
BoxWriter
    A class for writing a trajectory to a box file.
"""

from . import BaseWriter, BoxWriterError
from ..utils import instance_function_count_decorator
from ..traj import Trajectory, BoxFileFormat


class BoxWriter(BaseWriter):
    """
    A class for writing a trajectory to a box file.
    Inherits from BaseWriter. See BaseWriter for more information.

    It can write a trajectory to a box file in either a data file format or a VMD file format.
    """

    def __init__(self, filename: str | None = None, output_format: str | BoxFileFormat = 'data', mode='w') -> None:
        """
        Parameters
        ----------
        filename : str, optional
            The name of the file to write to. If None, the output is printed to stdout.
        output_format : str | BoxFileFormat, optional
            The format of the file. The default is 'data' i.e. BoxFileFormat.DATA.
        mode : str, optional
            The mode of the file. Either 'w' for write, 'a' for append or 'o' for overwrite.

        Raises
        ------
        ValueError
            If the given format is not in :py:class:`~PQAnalysis.traj.formats.BoxFileFormat`.
        """

        super().__init__(filename, mode)
        self.output_format = BoxFileFormat(output_format)

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
        if self.output_format == "vmd":
            self.write_vmd(traj)
        else:
            self.write_box_file(traj, reset_counter=reset_counter)

        self.close()

    def write_vmd(self, traj: Trajectory) -> None:
        """
        Writes the given trajectory to the file in VMD format.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to write.
        """
        self.__check_PBC__(traj)

        for frame in traj:
            cell = frame.cell

            print("8", file=self.file)
            print(
                f"Box   {cell.x} {cell.y} {cell.z}    {cell.alpha} {cell.beta} {cell.gamma}", file=self.file)
            edges = cell.bounding_edges
            for edge in edges:
                print(f"X   {edge[0]} {edge[1]} {edge[2]}", file=self.file)

    @instance_function_count_decorator
    def write_box_file(self, traj: Trajectory, reset_counter: bool = True) -> None:
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
        self.__check_PBC__(traj)

        counter = self.counter[BoxWriter.write_box_file.__name__]
        counter = len(traj)*(counter - 1)

        for i, frame in enumerate(traj):
            cell = frame.cell
            print(
                f"{counter + i+1} {cell.x} {cell.y} {cell.z} {cell.alpha} {cell.beta} {cell.gamma}", file=self.file)

    def __check_PBC__(self, traj: Trajectory) -> None:
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

        if not traj.check_PBC():
            raise BoxWriterError(
                "At least on cell of the trajectory is None. Cannot write box file.")
