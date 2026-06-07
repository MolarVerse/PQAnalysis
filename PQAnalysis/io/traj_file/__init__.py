"""
A subpackage to handle trajectory files.
"""

from .trajectory_reader import TrajectoryReader
from .trajectory_writer import TrajectoryWriter
from .frame_reader import (
    BaseFrameReader,
    XYZFrameReader,
    _FrameReader,
    get_frame_reader,
)
