"""
A subpackage to handle trajectory files.
"""

from typing import TYPE_CHECKING, Any

from PQAnalysis._lazy_import import public_dir, resolve_export

if TYPE_CHECKING:
    from .frame_reader import (
        BaseFrameReader,
        ExtXYZFrameReader,
        XYZFrameReader,
        _FrameReader,
        get_frame_reader,
    )
    from .raw_frame_reader import RawTrajectoryReader
    from .trajectory_reader import TrajectoryReader
    from .trajectory_writer import TrajectoryWriter

_EXPORTS = {
    "TrajectoryReader": ".trajectory_reader",
    "TrajectoryWriter": ".trajectory_writer",
    "RawTrajectoryReader": ".raw_frame_reader",
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
