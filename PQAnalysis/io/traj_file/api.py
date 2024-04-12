from beartype.typing import Generator

from PQAnalysis.io import (
    TrajectoryWriter,
    TrajectoryReader,
    FileWritingMode,
)
from PQAnalysis.traj import (
    Trajectory,
    MDEngineFormat,
    TrajectoryFormat,
    Frame,
)
from PQAnalysis.topology import Topology


def write_trajectory(traj,
                     filename: str | None = None,
                     format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF,
                     type: TrajectoryFormat | str = TrajectoryFormat.XYZ,
                     mode: FileWritingMode | str = FileWritingMode.WRITE,
                     ) -> None:
    """Wrapper for TrajectoryWriter to write a trajectory to a file.

    if format is None, the default PIMD-QMCF format is used. (see TrajectoryWriter.formats for available formats)
    if format is 'qmcfc', the QMCFC format is used (see TrajectoryWriter.formats for more information).

    Parameters
    ----------
    traj : Trajectory
        The trajectory to write.
    filename : str, optional
        The name of the file to write to. If None, the output is printed to stdout.
    format : MDEngineFormat | str, optional
        The format of the md engine for the output file. The default is MDEngineFormat.PIMD_QMCF.
    type : TrajectoryFormat | str, optional
        The type of the data to write to the file. Default is TrajectoryFormat.XYZ.
    mode  : FileWritingMode | str, optional
        The mode of the file. Either 'w' for write, 'a' for append or 'o' for overwrite. The default is 'w'.

    """

    writer = TrajectoryWriter(filename, format=format, mode=mode)
    writer.write(traj, type=type)


def read_trajectory(filename: str,
                    md_format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF,
                    traj_format: TrajectoryFormat | str = TrajectoryFormat.AUTO,
                    topology: Topology | None = None,
                    constant_topology: bool = True
                    ) -> Trajectory:
    """
    API function for reading a trajectory from a file.

    Parameters
    ----------
    filename : str
        The name of the file to read from.
    md_format : MDEngineFormat | str, optional
        The format of the trajectory, by default MDEngineFormat.PIMD_QMCF
    traj_format : TrajectoryFormat | str, optional
        The format of the trajectory, by default TrajectoryFormat.AUTO. The format is inferred from the file extension.
    topology : Topology | None, optional
        The topology of the trajectory, by default None
    constant_topology : bool, optional
        Whether the topology is constant over the trajectory or does change, by default True

    Returns
    -------
    Trajectory
        The trajectory read from the file.
    """

    reader = TrajectoryReader(
        filename,
        traj_format=traj_format,
        md_format=md_format,
        topology=topology,
        constant_topology=constant_topology
    )

    return reader.read()


def read_trajectory_generator(filename: str,
                              md_format: MDEngineFormat | str = MDEngineFormat.PIMD_QMCF,
                              traj_format: TrajectoryFormat | str = TrajectoryFormat.AUTO,
                              topology: Topology | None = None,
                              constant_topology: bool = True
                              ) -> Generator[Frame]:
    """
    API function for building a frame generator from a trajectory file.

    Parameters
    ----------
    filename : str
        The name of the file to read from.
    md_format : MDEngineFormat | str, optional
        The format of the trajectory, by default MDEngineFormat.PIMD_QMCF
    traj_format : TrajectoryFormat | str, optional
        The format of the trajectory, by default TrajectoryFormat.AUTO. The format is inferred from the file extension.
    topology : Topology | None, optional
        The topology of the trajectory, by default None
    constant_topology : bool, optional
        Whether the topology is constant over the trajectory or does change, by default True

    Returns
    -------
    Generator[Frame]
        A generator for the frames in the trajectory.
    """

    reader = TrajectoryReader(
        filename,
        traj_format=traj_format,
        md_format=md_format,
        topology=topology,
        constant_topology=constant_topology
    )

    return reader.frame_generator()
