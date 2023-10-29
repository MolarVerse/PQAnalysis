"""
A module containing the BoxWriter class and its associated methods.

...

Classes
-------
BoxWriter
    A class for writing a trajectory to a box file.

Methods
-------
write_box(traj, filename=None, format=None)
    Wrapper for BoxWriter to write a trajectory to a file.
"""

from typing import Union

from PQAnalysis.io.base import BaseWriter
from PQAnalysis.traj.trajectory import Trajectory


def write_box(traj, filename: Union[str, None] = None, format: Union[str, None] = None):
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

    Attributes
    ----------
    formats : list of str
        The available formats for the box file.
    format : str
        The format of the file. If None, the format is inferred as a data file format.
        (see BoxWriter.formats for available formats)

    Methods
    -------
    write(traj)
        Writes the given trajectory to the file.
    write_vmd(traj)
        Writes the given trajectory to the file in VMD format.
    write_box_file(traj)
        Writes the given trajectory to the file in data file format.
    """
    formats = [None, 'vmd']

    def __init__(self, filename: Union[str, None] = None, format: Union[str, None] = None, mode='w'):
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
                'Invalid format. Has to be either \'vmd\' or \'None\'.')

        self.format = format

    def write(self, traj: Trajectory):
        """
        Wrapper to write the given trajectory to the file.
        Depending on the format, either write_vmd or write_box_file is called.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to write.
        """

        self.open()
        if self.format == "vmd":
            self.write_vmd(traj)
        else:
            self.write_box_file(traj)

        self.close()

    def write_vmd(self, traj: Trajectory):
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
        self.__check_if_PBC__(traj)

        for frame in traj:
            cell = frame.cell

            print("8", file=self.file)
            print(
                f"Box   {cell.x} {cell.y} {cell.z}    {cell.alpha} {cell.beta} {cell.gamma}", file=self.file)
            edges = cell.bounding_edges
            for edge in edges:
                print(f"X   {edge[0]} {edge[1]} {edge[2]}", file=self.file)

    def write_box_file(self, traj: Trajectory):
        """
        Writes the given trajectory to the file in data file format.

        The format looks in general like this:

                1 1.0 1.0 1.0 90.0 90.0 90.0
                2 1.0 1.0 1.0 90.0 90.0 90.0
                ...
                n 1.1 1.1 1.1 90.0 90.0 90.0

        where the first column represents the step starting from 1, the second to fourth column
        represent the box vectors a, b, c, the fifth to seventh column represent the box angles.

        Parameters
        ----------
        traj : Trajectory
            The trajectory to write.

        Raises
        ------
        ValueError
            If the cell of a frame of the trajectory is None.
        """
        self.__check_if_PBC__(traj)

        for i, frame in enumerate(traj):
            cell = frame.cell
            print(
                f"{i+1} {cell.x} {cell.y} {cell.z} {cell.alpha} {cell.beta} {cell.gamma}")

    def __check_if_PBC__(self, traj: Trajectory):
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

        if not all(frame.PBC for frame in traj):
            raise ValueError(
                "Cell of trajectory is None. Cannot write box file.")
