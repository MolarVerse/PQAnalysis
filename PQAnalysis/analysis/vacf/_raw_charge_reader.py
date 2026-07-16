"""
A module containing a raw fast-path reader for charge (.chrg)
trajectory files.

The :py:class:`RawChargeTrajectoryReader` streams the charge values of
a charge trajectory as plain float64 numpy arrays, without building
:py:class:`~PQAnalysis.atomic_system.atomic_system.AtomicSystem`
objects for every frame. It is the scalar counterpart of
:py:class:`~PQAnalysis.io.traj_file.raw_frame_reader.RawTrajectoryReader`
and is used by the VACF analysis to read charge trajectories in
lockstep with the velocity trajectory. The charge values are bitwise
identical to the ones produced by
:py:meth:`~PQAnalysis.io.traj_file.trajectory_reader.TrajectoryReader.frame_generator`
(both parse the charges as correctly rounded float64 values).
"""

import logging

from beartype.typing import List

from PQAnalysis.io.base import BaseReader
from PQAnalysis.io.traj_file.exceptions import FrameReaderError
from PQAnalysis.io.traj_file.raw_frame_reader import RawTrajectoryReader
from PQAnalysis.traj import MDEngineFormat, TrajectoryFormat
from PQAnalysis.types import Np1DNumberArray
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

try:
    from ._vacf_kernel import parse_charge_lines  # pylint: disable=import-error
except ModuleNotFoundError:
    from ._vacf_kernel_py import parse_charge_lines



class RawChargeTrajectoryReader(RawTrajectoryReader):

    """
    A fast-path reader that streams the raw per-frame charge values of
    charge trajectory files.

    The reader shares the frame layout handling of
    :py:class:`~PQAnalysis.io.traj_file.raw_frame_reader.RawTrajectoryReader`
    (multiple files, QMCFC dummy atom stripping, cell caching and
    vacuum cell inheritance), but parses the scalar 'name charge' body
    lines of the
    :py:class:`~PQAnalysis.traj.formats.TrajectoryFormat.CHARGE`
    format instead of the xyz-family vector lines. The
    :py:meth:`~PQAnalysis.io.traj_file.raw_frame_reader.RawTrajectoryReader.raw_frame_generator`
    method therefore yields ``(values, cell)`` tuples with ``values``
    being the ``(n_atoms,)`` float64 array of the charge values of the
    frame.
    """

    # Set up the logger
    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

    @runtime_type_checking
    def __init__(  # pylint: disable=super-init-not-called
        self,
        filename: str | List[str],
        md_format: MDEngineFormat | str = MDEngineFormat.PQ,
    ) -> None:
        """
        Parameters
        ----------
        filename : str or list of str
            The name of the file to read from or a list of filenames
            to read from.
        md_format : MDEngineFormat | str, optional
            The format of the MD engine. Default is MDEngineFormat.PQ.
        """
        # Deliberately not calling RawTrajectoryReader.__init__ here:
        # the parent restricts itself to the xyz-family vector formats
        # while this reader always reads the CHARGE format.
        # pylint: disable-next=non-parent-init-called
        BaseReader.__init__(self, filename)  # pylint: disable=super-init-not-called

        if not self.multiple_files:
            self.filenames = [self.filename]

        self.traj_format = TrajectoryFormat.CHARGE
        self.md_format = MDEngineFormat(md_format)

        self._cell_cache = {}

    def _parse_body_lines(
        self,
        body_lines: List[str],
        n_atoms: int,
    ) -> Np1DNumberArray:
        """
        Parses the charge values of the frame body.

        Parameters
        ----------
        body_lines : List[str]
            The lines of the frame body (comment line + atom lines).
        n_atoms : int
            The number of atoms in the frame.

        Returns
        -------
        Np1DNumberArray
            The (n_atoms,) float64 array of charge values parsed from
            the atom lines.

        Raises
        ------
        FrameReaderError
            If an atom line cannot be parsed.
        """

        try:
            return parse_charge_lines(body_lines[1:], n_atoms)
        except ValueError:
            self.logger.error(
                'Invalid file format in scalar values of Frame.',
                exception=FrameReaderError,
            )

        return None  # pragma: no cover - logger.error raises

    def _strip_dummy_atom(
        self,
        body_lines: List[str],
        values: Np1DNumberArray,
    ) -> Np1DNumberArray:
        """
        Strips the leading QMCFC dummy atom row from the values.

        Same semantics as the base class method, but typed for the
        one-dimensional charge values of this reader.

        Parameters
        ----------
        body_lines : List[str]
            The lines of the frame body (comment line + atom lines).
        values : Np1DNumberArray
            The parsed charge values of the frame body.

        Returns
        -------
        Np1DNumberArray
            The charge values without the leading dummy atom entry.

        Raises
        ------
        FrameReaderError
            If the first atom of the frame is not X.
        """

        first_atom = body_lines[1].split(None, 1)[0]

        if first_atom.upper() != 'X':
            self.logger.error(
                (
                    'The first atom in one of the frames is not X. '
                    'Please use PQ (default) md engine instead'
                ),
                exception=FrameReaderError,
            )

        return values[1:]
