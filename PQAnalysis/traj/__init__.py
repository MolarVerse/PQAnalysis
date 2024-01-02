"""
A package containing classes and functions to handle molecular dynamics trajectories.

The traj package contains the following submodules:
    
        - frame
        - trajectory
        - formats
        
The traj package contains the following classes:
            
        - Frame
        - Trajectory
        - TrajectoryFormat
        - MDEngineFormat
            
The traj package contains the following exceptions:
        
        - FrameError
        - TrajectoryError
        - TrajectoryFormatError
        - MDEngineFormatError
"""

from .exceptions import FrameError, TrajectoryError, TrajectoryFormatError, MDEngineFormatError


from .formats import TrajectoryFormat, MDEngineFormat
from .frame import Frame
from .trajectory import Trajectory
