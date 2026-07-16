"""
A subpackage to handle trajectory files.
"""

from .trajectory_reader import TrajectoryReader
from .trajectory_writer import TrajectoryWriter
from .raw_frame_reader import RawTrajectoryReader
from .frame_reader import (
    BaseFrameReader,
    XYZFrameReader,
    ExtXYZFrameReader,
    _FrameReader,
    get_frame_reader,
)
