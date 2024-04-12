"""
A package containing classes and functions to handle molecular dynamics trajectories.
"""

from .frame.exceptions import FrameError
from .frame.frame import Frame

from .exceptions import TrajectoryError, TrajectoryFormatError, MDEngineFormatError

from .formats import TrajectoryFormat, MDEngineFormat
from .trajectory import Trajectory

from .api import *
