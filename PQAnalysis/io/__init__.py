"""
A package containing classes and functions to handle input
and output of molecular dynamics simulations.
"""

# import the formats from the formats module
from .formats import BoxFileFormat, FileWritingMode, OutputFileFormat

# import the classes from the base module
from .base import BaseReader, BaseWriter

# import the classes from the moldescriptorReader module
from .moldescriptor_reader import MoldescriptorReader

# import the classes from the restart_file subpackage
from .restart_file.restart_writer import RestartFileWriter
from .restart_file.restart_reader import RestartFileReader
from .restart_file.api import read_restart_file

# import the classes from the traj_file subpackage
from .traj_file.trajectory_reader import TrajectoryReader
from .traj_file.trajectory_writer import TrajectoryWriter
from .traj_file.frame_reader import _FrameReader
from .traj_file.api import (
    read_trajectory,
    write_trajectory,
    read_trajectory_generator,
    calculate_frames_of_trajectory_file,
)

# import the classes from the gen_file subpackage
from .gen_file.gen_file_reader import GenFileReader
from .gen_file.gen_file_writer import GenFileWriter
from .gen_file.api import read_gen_file, write_gen_file

# import the classes from the topology_file subpackage
from .topology_file import TopologyFileReader
from .topology_file import TopologyFileWriter
from .topology_file.api import read_topology_file, write_topology_file

from .info_file_reader import InfoFileReader
from .energy_file_reader import EnergyFileReader
from .box_writer import BoxWriter

from .input_file_reader import InputFileParser
from .input_file_reader import PQInputFileReader
from .input_file_reader import PQAnalysisInputFileReader
from .input_file_reader import InputFileFormat

from .api import continue_input_file
from .conversion_api import (gen2xyz, xyz2gen, rst2xyz, traj2box, traj2qmcfc)
from .write_api import write, write_box
