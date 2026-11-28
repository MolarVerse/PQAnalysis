"""
A subpackage to handle trajectory files.
"""

from typing import TYPE_CHECKING, Any

from PQAnalysis._lazy_import import public_dir, resolve_export

if TYPE_CHECKING:  # pragma: no cover
    from .frame_reader import (
        BaseFrameReader,
        ExtXYZFrameReader,
        XYZFrameReader,
        _FrameReader,
        get_frame_reader,
    )
    from .ase_traj_reader import AseTrajReader
    from .h5md_reader import H5MDReader
    from .raw_frame_reader import RawTrajectoryReader
    from .trajectory_reader import TrajectoryReader
    from .trajectory_writer import TrajectoryWriter
    from .zarr_traj_reader import ZarrTrajReader

_EXPORTS = {
    "TrajectoryReader": ".trajectory_reader",
    "TrajectoryWriter": ".trajectory_writer",
    "RawTrajectoryReader": ".raw_frame_reader",
    "H5MDReader": ".h5md_reader",
    "AseTrajReader": ".ase_traj_reader",
    "ZarrTrajReader": ".zarr_traj_reader",
    "BaseFrameReader": ".frame_reader",
    "XYZFrameReader": ".frame_reader",
    "ExtXYZFrameReader": ".frame_reader",
    "_FrameReader": ".frame_reader",
    "get_frame_reader": ".frame_reader",
}

__all__ = [name for name in _EXPORTS if not name.startswith("_")]



def __getattr__(name: str) -> Any:
    return resolve_export(__name__, globals(), _EXPORTS, name)



def __dir__() -> list[str]:
    return public_dir(globals(), _EXPORTS)
