"""
A package containing classes and functions to handle molecular dynamics trajectories.
"""

from typing import TYPE_CHECKING, Any

from PQAnalysis._lazy_import import public_dir, resolve_export

if TYPE_CHECKING:  # pragma: no cover
    from .api import check_trajectory_pbc, check_trajectory_vacuum
    from .exceptions import (
        MDEngineFormatError,
        TrajectoryError,
        TrajectoryFormatError,
    )
    from .formats import MDEngineFormat, TrajectoryFormat
    from .trajectory import Trajectory

_EXPORTS = {
    "TrajectoryError": ".exceptions",
    "TrajectoryFormatError": ".exceptions",
    "MDEngineFormatError": ".exceptions",
    "TrajectoryFormat": ".formats",
    "MDEngineFormat": ".formats",
    "Trajectory": ".trajectory",
    "check_trajectory_pbc": ".api",
    "check_trajectory_vacuum": ".api",
}

__all__ = list(_EXPORTS)



def __getattr__(name: str) -> Any:
    return resolve_export(__name__, globals(), _EXPORTS, name)



def __dir__() -> list[str]:
    return public_dir(globals(), _EXPORTS)
