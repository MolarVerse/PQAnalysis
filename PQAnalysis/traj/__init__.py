"""
A package containing classes and functions to handle molecular dynamics trajectories.
"""

from .exceptions import FrameError, TrajectoryError, TrajectoryFormatError, MDEngineFormatError


from .formats import TrajectoryFormat, MDEngineFormat
from .frame import Frame
from .trajectory import Trajectory

from .api import *
