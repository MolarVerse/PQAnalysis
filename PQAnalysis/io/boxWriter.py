"""
A module containing the BoxWriter class and its associated methods.

...

Classes
-------
BoxWriter
    A class for writing a trajectory to a box file.
"""

from .base import BaseWriter
from ..utils.decorators import instance_function_count_decorator
from ..traj.trajectory import Trajectory


def write_box(traj, filename: str | None = None, format: str | None = None) -> None:
    '''
    Wrapper for BoxWriter to write a trajectory to a file.

    Parameters
    ----------
    traj : Trajectory
        The trajectory to write.
    filename : str, optional
        The name of the file to write to. If None, the output is printed to stdout.
    format : str, optional
        The format of the file. If None, the format is inferred as a data file format.
        (see BoxWriter.formats for available formats)
    '''

    writer = BoxWriter(filename, format)
    writer.write(traj)


class BoxWriter(BaseWriter):
    """
    A class for writing a trajectory to a box file.
    Inherits from BaseWriter. See BaseWriter for more information.

    It can write a trajectory to a box file in either a data file format or a VMD file format.

    ...

    Class Attributes
    ----------------
    formats : list of str
        The available formats for the box file.

    Attributes
    ----------
    format : str
        The format of the file. If None, the format is inferred as a data file format.
        (see BoxWriter.formats for available formats)
    """
    formats = [None, 'data', 'vmd']

    def __init__(self, filename: str | None = None, format: str | None = None, mode='w') -> None:
        """
        Initializes the BoxWriter with the given filename, format and mode.

        Parameters
        ----------
        filename : str, optional
            The name of the file to write to. If None, the output is printed to stdout.
        format : str, optional
            The format of the file. If None, the format is inferred as a data file format.
            (see BoxWriter.formats for available formats)
        mode : str, optional
            The mode of the file. Either 'w' for write or 'a' for append.

        Raises
        ------
        ValueError
            If the given mode is not 'w' or 'a'.
        """

        super().__init__(filename, mode)
        if format not in self.formats:
            raise ValueError(
                'Invalid format. Has to be either \'vmd\', \'data\' or \'None\'.')

        if format is None:
            format = 'data'

        self.format = format

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
        if self.format == "vmd":
            self.write_vmd(traj)
        else:
            self.write_box_file(traj, reset_counter=reset_counter)

        self.close()

    def write_vmd(self, traj: Trajectory) -> None:
        """
        Writes the given trajectory to the file in VMD format.

        The format looks in general like this:

                8
                Box  1.0 1.0 1.0    90.0 90.0 90.0
                X   0.0 0.0 0.0
                X   1.0 0.0 0.0
                X   0.0 1.0 0.0
                X   1.0 1.0 0.0
                X   0.0 0.0 1.0
                X   1.0 0.0 1.0
                X   0.0 1.0 1.0
                X   1.0 1.0 1.0
                8
                Box  1.0 1.0 1.0    90.0 90.0 90.0
                X   0.0 0.0 0.0
                ...

        where all X represent the vertices of the box. The first line contains the number of vertices.
        The second line contains the box dimensions and box angles as the comment line for a xyz file.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to write.

        Raises
        ------
        ValueError
            If the cell of a frame of the trajectory is None.
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

        The format looks in general like this:

                1 1.0 1.0 1.0 90.0 90.0 90.0
                2 1.0 1.0 1.0 90.0 90.0 90.0
                ...
                n 1.1 1.1 1.1 90.0 90.0 90.0

        where the first column represents the step starting from 1, the second to fourth column
        represent the box vectors a, b, c, the fifth to seventh column represent the box angles.

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
        ValueError
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
            raise ValueError(
                "At least on cell of the trajectory is None. Cannot write box file.")
