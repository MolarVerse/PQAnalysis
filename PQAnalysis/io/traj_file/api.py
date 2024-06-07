"""
A module containing different API functions for reading and writing trajectory files.
"""

from beartype.typing import Generator, List

from PQAnalysis.topology import Topology
from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.type_checking import runtime_type_checking

from PQAnalysis.io.traj_file import (
    TrajectoryWriter,
    TrajectoryReader,
)

from PQAnalysis.traj import (
    Trajectory,
    MDEngineFormat,
    TrajectoryFormat,
)



@runtime_type_checking
def write_trajectory(
    traj: Trajectory | AtomicSystem,
    filename: str | None = None,
    engine_format: MDEngineFormat | str = MDEngineFormat.PQ,
    traj_type: TrajectoryFormat | str = TrajectoryFormat.XYZ,
    mode: FileWritingMode | str = FileWritingMode.WRITE,
) -> None:
    """Wrapper for TrajectoryWriter to write a trajectory to a file.

    if format is None, the default PQ format is used
    (see TrajectoryWriter.formats for available formats).
    if format is 'qmcfc', the QMCFC format is used
    (see TrajectoryWriter.formats for more information).

    Parameters
    ----------
    traj : Trajectory
        The trajectory to write.
    filename : str, optional
        The name of the file to write to. If None, the output is printed to stdout.
    engine_format : MDEngineFormat | str, optionalsssss
        The format of the md engine for the output file.
        The default is MDEngineFormat.PQ.
    traj_type : TrajectoryFormat | str, optional
        The type of the data to write to the file.
        Default is TrajectoryFormat.XYZ.
    mode : FileWritingMode | str, optional
        The mode of the file. Either 'w' for write, 
        'a' for append or 'o' for overwrite. The default is 'w'.

    """

    writer = TrajectoryWriter(filename, engine_format=engine_format, mode=mode)
    writer.write(traj, traj_type=traj_type)



@runtime_type_checking
def read_trajectory(
    filename: str | List[str],
    md_format: MDEngineFormat | str = MDEngineFormat.PQ,
    traj_format: TrajectoryFormat | str = TrajectoryFormat.AUTO,
    topology: Topology | None = None,
    constant_topology: bool = True
) -> Trajectory:
    """
    API function for reading a trajectory from a file.

    Parameters
    ----------
    filename : str | List[str]
        The name of the file to read from or a list of files to read from.
    md_format : MDEngineFormat | str, optional
        The format of the trajectory, by default MDEngineFormat.PQ
    traj_format : TrajectoryFormat | str, optional
        The format of the trajectory, by default TrajectoryFormat.AUTO.
        The format is inferred from the file extension.
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



@runtime_type_checking
def calculate_frames_of_trajectory_file(filename: str):
    """
    Calculate the number of frames in a trajectory file.

    Parameters
    ----------
    filename : str
        The name of the file to read.

    Returns
    -------
    int
        The number of frames in the trajectory file.
    """
    reader = TrajectoryReader(filename)
    return sum(reader.calculate_number_of_frames_per_file())



@runtime_type_checking
def read_trajectory_generator(
    filename: str,
    md_format: MDEngineFormat | str = MDEngineFormat.PQ,
    traj_format: TrajectoryFormat | str = TrajectoryFormat.AUTO,
    topology: Topology | None = None,
    constant_topology: bool = True
) -> Generator[AtomicSystem]:
    """
    API function for building a frame generator from a trajectory file.

    Parameters
    ----------
    filename : str
        The name of the file to read from.
    md_format : MDEngineFormat | str, optional
        The format of the trajectory, by default MDEngineFormat.PQ
    traj_format : TrajectoryFormat | str, optional
        The format of the trajectory, by default TrajectoryFormat.AUTO.
        The format is inferred from the file extension.
    topology : Topology | None, optional
        The topology of the trajectory, by default None
    constant_topology : bool, optional
        Whether the topology is constant over the trajectory or does change, by default True

    Returns
    -------
    Generator[AtomicSystem]
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
