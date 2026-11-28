"""
A module containing a reader for Zarr trajectory stores.

A Zarr trajectory is a directory ending in ``.zarr`` with the H5MD layout
of PQ, see https://zarr.dev/. Reading it needs the optional dependency
``zarr`` (``pip install PQAnalysis[zarr]``).
"""

import contextlib
import logging

from beartype.typing import Any

from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__

from .exceptions import ZarrTrajReaderError
from .h5md_reader import H5MDReader



class ZarrTrajReader(H5MDReader):

    """
    A class for reading the frames of a single Zarr trajectory store.

    The store has the groups, arrays, attributes and units of the H5MD
    files written by PQ. Values are converted to the same PQ units as
    the :py:class:`~PQAnalysis.io.traj_file.h5md_reader.H5MDReader` does.
    Frames are read lazily, one shard or chunk of frames at a time.

    Generally, this class should not be used directly. Instead, use the
    :py:class:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader`
    with ``traj_format="zarr"`` or a ``.zarr`` store.
    """

    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    _format_name = "Zarr"

    _error = ZarrTrajReaderError

    @staticmethod
    def _import_backend() -> Any:
        """
        Imports zarr, which is an optional dependency.
        """
        try:
            import zarr  # pylint: disable=import-outside-toplevel
        except ModuleNotFoundError:
            ZarrTrajReader.logger.error(
                (
                    "Reading Zarr trajectories requires zarr. Install it "
                    "with 'pip install PQAnalysis[zarr]'."
                ),
                exception=ZarrTrajReaderError,
            )

        return zarr

    def _open(self) -> Any:
        """
        Opens the Zarr store for reading.
        """
        try:
            return contextlib.nullcontext(
                self._backend.open_group(self.filename, mode="r")
            )
        except (OSError, ValueError):
            self.logger.error(
                f"Could not open the Zarr store {self.filename}.",
                exception=ZarrTrajReaderError,
            )

        return None  # pragma: no cover

    def _is_dataset(self, node: Any) -> bool:
        """
        Whether a node of the store is an array.
        """
        return isinstance(node, self._backend.Array)

    @staticmethod
    def _read_atom_names(file: Any) -> list[str] | None:
        """
        Reads the atom names, which PQ stores as an attribute of
        ``/parameters``.
        """
        if "parameters" in file:
            names = file["parameters"].attrs.get("atom_names")
            if names is not None:
                return [str(name) for name in names]

        return H5MDReader._read_atom_names(file)

    @staticmethod
    def _frames(dataset: Any) -> "_ChunkedFrames":
        """
        Returns the array indexed per frame, read one chunk at a time.
        """
        return _ChunkedFrames(dataset)



class _ChunkedFrames:

    """
    Reads the frames of a Zarr array one shard or chunk at a time.
    """

    def __init__(self, array: Any) -> None:
        self.array = array
        self.chunk_size = int((array.shards or array.chunks)[0])
        self.start = None
        self.values = None

    def __getitem__(self, index: int) -> Any:
        start = index - index % self.chunk_size
        if start != self.start:
            self.values = self.array[start:start + self.chunk_size]
            self.start = start

        return self.values[index - start]
