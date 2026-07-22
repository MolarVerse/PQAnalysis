"""
Tests of the header-only cell scan of the ADF raw-frame fast path.

The fast path collects the per-frame cells with a cheap header-only
scan (:py:meth:`~PQAnalysis.analysis.adf.adf.ADF._scan_cells`) instead
of building an AtomicSystem per frame. These tests assert that the scan
reproduces the cells of the full
:py:class:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader`
scan for constant, NPT, vacuum-inheriting and pure-vacuum boxes (also
across file boundaries) and that an invalid box header raises the same
error as the reader.
"""

import numpy as np
from beartype.roar import BeartypeCallHintReturnViolation

from PQAnalysis.analysis import ADF
from PQAnalysis.io import TrajectoryReader
from PQAnalysis.io.traj_file.exceptions import TrajectoryReaderError

from .. import pytestmark  # pylint: disable=unused-import
from ...conftest import assert_logging_with_exception

# pylint: disable=protected-access



def _write_trajectory(path, positions, headers):
    """
    Writes an xyz trajectory file with the given per-frame header box
    strings.

    Parameters
    ----------
    path : os.PathLike
        The path of the trajectory file to write.
    positions : np.ndarray
        The per-frame atom positions of shape (n_frames, n_atoms, 3).
    headers : List[str]
        The per-frame box header suffixes appended after the atom count.

    Returns
    -------
    str
        The string path of the written trajectory file.
    """
    with open(path, "w", encoding="utf-8") as file:
        for frame_positions, header in zip(positions, headers):
            file.write(f"{len(frame_positions)}{header}\n\n")

            for i, (x, y, z) in enumerate(frame_positions):
                name = "O" if i % 2 == 0 else "H"
                file.write(f"{name} {x:.6f} {y:.6f} {z:.6f}\n")

    return str(path)


def _write_npt_trajectory(path, n_frames=24, n_atoms=16, seed=99):
    """
    Writes an NPT-like xyz trajectory file with per-frame changing
    boxes, mixing orthorhombic and triclinic headers.

    Parameters
    ----------
    path : os.PathLike
        The path of the trajectory file to write.
    n_frames : int, optional
        The number of frames to write, by default 24.
    n_atoms : int, optional
        The number of atoms per frame, by default 16.
    seed : int, optional
        The random seed for the positions, by default 99.

    Returns
    -------
    str
        The string path of the written trajectory file.
    """
    rng = np.random.default_rng(seed)

    positions = rng.uniform(0.0, 12.0, size=(n_frames, n_atoms, 3))

    headers = []
    for i in range(n_frames):
        factor = 1.0 + 0.01 * (i % 4)
        if (i // 4) % 2 == 0:
            headers.append(f" {12.0 * factor} {13.0 * factor} {14.0}")
        else:
            headers.append(
                f" {12.0 * factor} {13.0 * factor} {14.0} 80.0 95.0 103.0"
            )

    return _write_trajectory(path, positions, headers)


def _write_random_trajectory(path, n_frames=25, n_atoms=16, seed=4242):
    """
    Writes a constant orthorhombic box xyz trajectory file.

    Parameters
    ----------
    path : os.PathLike
        The path of the trajectory file to write.
    n_frames : int, optional
        The number of frames to write, by default 25.
    n_atoms : int, optional
        The number of atoms per frame, by default 16.
    seed : int, optional
        The random seed for the positions, by default 4242.

    Returns
    -------
    str
        The string path of the written trajectory file.
    """
    rng = np.random.default_rng(seed)

    positions = rng.uniform(0.0, 12.0, size=(n_frames, n_atoms, 3))
    headers = [" 12.0 12.0 12.0"] * n_frames

    return _write_trajectory(path, positions, headers)



class TestADFScanCells:

    """
    Tests of the header-only cell scan of the fast path against the
    cells full scan of the TrajectoryReader.
    """

    @staticmethod
    def _assert_cells_match(scanned, reference):
        """
        Asserts that two cell lists have identical box lengths and angles.

        Parameters
        ----------
        scanned : Cells
            The cells produced by the header-only scan.
        reference : Cells
            The reference cells of the TrajectoryReader.
        """
        assert len(scanned) == len(reference)

        for scanned_cell, reference_cell in zip(scanned, reference):
            assert np.array_equal(
                scanned_cell.box_lengths, reference_cell.box_lengths
            )
            assert np.array_equal(
                scanned_cell.box_angles, reference_cell.box_angles
            )

    def test_matches_trajectory_reader_cells_npt(self, tmp_path):
        """The scan reproduces the NPT reader cells and deduplicates."""
        filename = _write_npt_trajectory(tmp_path / "traj.xyz")

        cells, unique_cells = ADF._scan_cells([filename])

        self._assert_cells_match(cells, TrajectoryReader(filename).cells)

        # every distinct box appears exactly once in the unique list
        assert len({id(cell) for cell in cells}) == len(unique_cells)

    def test_matches_trajectory_reader_cells_dedup(self, tmp_path):
        """A constant box collapses to a single shared unique cell."""
        filename = _write_random_trajectory(tmp_path / "traj.xyz")

        cells, unique_cells = ADF._scan_cells([filename])

        self._assert_cells_match(cells, TrajectoryReader(filename).cells)

        assert len(unique_cells) == 1
        assert all(cell is unique_cells[0] for cell in cells)

    def test_matches_trajectory_reader_cells_vacuum_inheritance(
        self, tmp_path
    ):
        """Frames without a box header inherit the previous box."""
        rng = np.random.default_rng(11)
        positions = rng.uniform(0.0, 10.0, size=(6, 4, 3))

        # frames without box information inherit the last box
        headers = [
            " 10.0 10.0 10.0",
            "",
            " 11.0 11.0 11.0 80.0 95.0 103.0",
            "",
            "",
            " 10.0 10.0 10.0",
        ]
        filename = _write_trajectory(tmp_path / "traj.xyz", positions, headers)

        cells, unique_cells = ADF._scan_cells([filename])

        self._assert_cells_match(cells, TrajectoryReader(filename).cells)

        assert len(unique_cells) == 2
        assert cells[1] is cells[0]
        assert cells[3] is cells[2]
        assert cells[4] is cells[2]
        assert cells[5] is cells[0]

    def test_matches_trajectory_reader_cells_pure_vacuum(self, tmp_path):
        """A trajectory without any box header maps to one vacuum cell."""
        rng = np.random.default_rng(12)
        positions = rng.uniform(0.0, 10.0, size=(3, 4, 3))

        filename = _write_trajectory(
            tmp_path / "traj.xyz", positions, ["", "", ""]
        )

        cells, unique_cells = ADF._scan_cells([filename])

        self._assert_cells_match(cells, TrajectoryReader(filename).cells)

        assert len(unique_cells) == 1
        assert unique_cells[0].is_vacuum

    def test_multiple_files_inherit_across_boundaries(self, tmp_path):
        """A missing box header inherits the last box across files."""
        rng = np.random.default_rng(13)
        positions = rng.uniform(0.0, 10.0, size=(2, 4, 3))

        filename1 = _write_trajectory(
            tmp_path / "traj1.xyz", positions, [" 10.0 10.0 10.0", ""]
        )
        filename2 = _write_trajectory(
            tmp_path / "traj2.xyz", positions, ["", " 11.0 11.0 11.0"]
        )

        cells, unique_cells = ADF._scan_cells([filename1, filename2])

        reader = TrajectoryReader([filename1, filename2])
        self._assert_cells_match(cells, reader.cells)

        assert len(cells) == 4
        assert cells[2] is cells[0]
        assert len(unique_cells) == 2

    def test_invalid_box_header(self, tmp_path, caplog):
        """An invalid box header raises the reader error via the scan."""
        rng = np.random.default_rng(14)
        positions = rng.uniform(0.0, 10.0, size=(2, 4, 3))

        filename = _write_trajectory(
            tmp_path / "traj.xyz",
            positions,
            [" 10.0 10.0 10.0", " 10.0 10.0"],
        )

        assert_logging_with_exception(
            caplog=caplog,
            logging_name=ADF.__qualname__,
            logging_level="ERROR",
            message_to_test=(
                "Invalid number of arguments for box: 3 encountered "
                f"in file {filename}:2 = 4 10.0 10.0"
            ),
            exception=TrajectoryReaderError,
            function=ADF._scan_cells,
            filenames=[filename],
        )

    def test_resolve_header_cell_returns_none_when_error_suppressed(
        self, monkeypatch
    ):
        """
        The invalid-header branch of :py:meth:`ADF._resolve_header_cell`
        returns ``None`` after logging.

        ``logger.error`` normally raises for an invalid box header, so
        the trailing ``return None`` (a defensive fall-through) is only
        reached when the logger is silenced. Patching the logger to a
        no-op exercises that final statement. In RELEASE mode the call
        returns ``None``; under the DEBUG beartype claw hook the same
        executed ``return None`` is rejected by the wrapper's ``-> Cell``
        return check (which still means the statement ran).
        """
        monkeypatch.setattr(ADF.logger, "error", lambda *args, **kwargs: None)

        state = {
            "cells": [],
            "unique_cells": [],
            "cell_cache": {},
            "last_cell": None,
            "vacuum_cell": None,
        }

        # a two-token header is neither a bare atom count (1 token) nor a
        # valid 4-/7-token box, so it takes the invalid-header branch
        try:
            result = ADF._resolve_header_cell(
                "4 10.0 10.0", "traj.xyz", 2, state
            )
        except BeartypeCallHintReturnViolation:
            # DEBUG mode: the executed ``return None`` fails the
            # runtime ``-> Cell`` return-type check of the claw hook
            pass
        else:
            assert result is None
