"""
A package containing classes and functions to handle input
and output of molecular dynamics simulations.
"""

from typing import TYPE_CHECKING, Any

from PQAnalysis._lazy_import import public_dir, resolve_export

if TYPE_CHECKING:
    from .api import continue_input_file
    from .base import BaseReader, BaseWriter
    from .box_reader import BoxReader, read_box
    from .box_writer import BoxWriter
    from .conversion_api import (
        gen2xyz,
        rst2xyz,
        traj2box,
        traj2extxyz,
        traj2qmcfc,
        xyz2gen,
        xyz2rst,
    )
    from .energy_file_reader import EnergyFileReader
    from .formats import (
        BoxFileFormat,
        ExtXYZProfile,
        FileWritingMode,
        OutputFileFormat,
    )
    from .gen_file.api import read_gen_file, write_gen_file
    from .gen_file.gen_file_reader import GenFileReader
    from .gen_file.gen_file_writer import GenFileWriter
    from .info_file_reader import InfoFileReader
    from .input_file_reader.formats import InputFileFormat
    from .input_file_reader.input_file_parser import InputFileParser
    from .input_file_reader.pq.pq_input_file_reader import PQInputFileReader
    from .input_file_reader.pq_analysis.pqanalysis_input_file_reader import (
        PQAnalysisInputFileReader,
    )
    from .moldescriptor_reader import MoldescriptorReader
    from .optimizer_file_reader import OptimizerFileReader, read_optimizer_file
    from .restart_file.api import read_restart_file
    from .restart_file.restart_reader import RestartFileReader
    from .restart_file.restart_writer import RestartFileWriter
    from .topology_file.api import read_topology_file, write_topology_file
    from .topology_file.topology_file_reader import TopologyFileReader
    from .topology_file.topology_file_writer import TopologyFileWriter
    from .traj_file.api import (
        calculate_frames_of_trajectory_file,
        read_trajectory,
        read_trajectory_generator,
        write_trajectory,
    )
    from .traj_file.frame_reader import (
        BaseFrameReader,
        ExtXYZFrameReader,
        XYZFrameReader,
        _FrameReader,
        get_frame_reader,
    )
    from .traj_file.raw_frame_reader import RawTrajectoryReader
    from .traj_file.trajectory_reader import TrajectoryReader
    from .traj_file.trajectory_writer import TrajectoryWriter
    from .write_api import write, write_box

_EXPORTS = {
    "BoxFileFormat": ".formats",
    "ExtXYZProfile": ".formats",
    "FileWritingMode": ".formats",
    "OutputFileFormat": ".formats",
    "BaseReader": ".base",
    "BaseWriter": ".base",
    "MoldescriptorReader": ".moldescriptor_reader",
    "RestartFileWriter": ".restart_file.restart_writer",
    "RestartFileReader": ".restart_file.restart_reader",
    "read_restart_file": ".restart_file.api",
    "TrajectoryReader": ".traj_file.trajectory_reader",
    "TrajectoryWriter": ".traj_file.trajectory_writer",
    "RawTrajectoryReader": ".traj_file.raw_frame_reader",
    "BaseFrameReader": ".traj_file.frame_reader",
    "XYZFrameReader": ".traj_file.frame_reader",
    "ExtXYZFrameReader": ".traj_file.frame_reader",
    "_FrameReader": ".traj_file.frame_reader",
    "get_frame_reader": ".traj_file.frame_reader",
    "read_trajectory": ".traj_file.api",
    "write_trajectory": ".traj_file.api",
    "read_trajectory_generator": ".traj_file.api",
    "calculate_frames_of_trajectory_file": ".traj_file.api",
    "GenFileReader": ".gen_file.gen_file_reader",
    "GenFileWriter": ".gen_file.gen_file_writer",
    "read_gen_file": ".gen_file.api",
    "write_gen_file": ".gen_file.api",
    "TopologyFileReader": ".topology_file.topology_file_reader",
    "TopologyFileWriter": ".topology_file.topology_file_writer",
    "read_topology_file": ".topology_file.api",
    "write_topology_file": ".topology_file.api",
    "InfoFileReader": ".info_file_reader",
    "EnergyFileReader": ".energy_file_reader",
    "BoxReader": ".box_reader",
    "read_box": ".box_reader",
    "OptimizerFileReader": ".optimizer_file_reader",
    "read_optimizer_file": ".optimizer_file_reader",
    "BoxWriter": ".box_writer",
    "InputFileParser": ".input_file_reader.input_file_parser",
    "PQInputFileReader": ".input_file_reader.pq.pq_input_file_reader",
    "PQAnalysisInputFileReader": (
        ".input_file_reader.pq_analysis.pqanalysis_input_file_reader"
    ),
    "InputFileFormat": ".input_file_reader.formats",
    "continue_input_file": ".api",
    "gen2xyz": ".conversion_api",
    "xyz2gen": ".conversion_api",
    "rst2xyz": ".conversion_api",
    "xyz2rst": ".conversion_api",
    "traj2box": ".conversion_api",
    "traj2qmcfc": ".conversion_api",
    "traj2extxyz": ".conversion_api",
    "write": ".write_api",
    "write_box": ".write_api",
}

__all__ = [name for name in _EXPORTS if not name.startswith("_")]



def __getattr__(name: str) -> Any:
    return resolve_export(__name__, globals(), _EXPORTS, name)



def __dir__() -> list[str]:
    return public_dir(globals(), _EXPORTS)
