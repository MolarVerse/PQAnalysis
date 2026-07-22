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
from PQAnalysis.io.traj_file._slab_parser_py import MODE_CHARGE
from PQAnalysis.io.traj_file.raw_frame_reader import RawTrajectoryReader
from PQAnalysis.traj import MDEngineFormat, TrajectoryFormat
from PQAnalysis.utils.custom_logging import setup_logger
from PQAnalysis import __package_name__
from PQAnalysis.type_checking import runtime_type_checking

# The charge body lines are parsed by the shared slab parser of
# RawTrajectoryReader (MODE_CHARGE). parse_charge_lines remains the
# reference scalar line parser of the vacf kernels and is re-exported
# here for the kernel wiring tests.
try:
    from ._vacf_kernel import parse_charge_lines  # pylint: disable=import-error,unused-import
except ModuleNotFoundError:
    from ._vacf_kernel_py import parse_charge_lines  # pylint: disable=unused-import



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

    #: The slab parser body mode of this reader (a name token plus
    #: exactly one float64 value per atom line).
    _SLAB_MODE = MODE_CHARGE

    #: The error message used when a frame body line cannot be parsed.
    _BODY_ERROR_MESSAGE = 'Invalid file format in scalar values of Frame.'

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
