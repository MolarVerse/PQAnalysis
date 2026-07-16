"""
API functions for the momentum analysis.
"""

from beartype.typing import List

from PQAnalysis.io import TrajectoryReader
from PQAnalysis.io.formats import FileWritingMode
from PQAnalysis.topology import SelectionCompatible
from PQAnalysis.traj import MDEngineFormat
from PQAnalysis.type_checking import runtime_type_checking
from PQAnalysis.types import Np1DNumberArray, PositiveReal

from .momentum import Momentum
from .momentum_output_file_writer import MomentumDataWriter



@runtime_type_checking
def check_momentum(
    trajectory_files: str | List[str],
    output: str | None = None,
    selection: SelectionCompatible = None,
    use_full_atom_info: bool = False,
    scale: PositiveReal | None = None,
    md_format: MDEngineFormat | str = MDEngineFormat.PQ,
    mode: str | FileWritingMode = "w",
) -> Np1DNumberArray:
    """
    Calculate the total linear momentum norm per frame and write it.

    Reads the given velocity trajectory file(s) frame by frame,
    calculates the norm of the total linear momentum
    ``P = sum_i m_i * v_i`` of the selected atoms for every frame and
    writes one row per frame containing the one-based frame index and
    the scaled momentum norm. This is a port of the legacy
    ``equipartition.jl`` tool.

    Parameters
    ----------
    trajectory_files : str | List[str]
        The velocity trajectory file(s) to read.
    output : str | None, optional
        The output file. If None, the output is printed to stdout,
        by default None.
    selection : SelectionCompatible, optional
        The selection of atoms to include in the total momentum,
        by default None (all atoms).
    use_full_atom_info : bool, optional
        Whether to use the full atom information of the trajectory
        for the selection or not, by default False.
    scale : PositiveReal | None, optional
        The scaling factor applied to the momentum norm before
        output, by default None (1e-15, which converts
        amu*Angstrom/s to amu*Angstrom/fs).
    md_format : MDEngineFormat | str, optional
        The format of the trajectory, by default MDEngineFormat.PQ.
        Use ``qmcfc`` for legacy trajectories with a leading dummy
        'X' atom.
    mode : str | FileWritingMode, optional
        The writing mode of the output file, by default "w".

    Returns
    -------
    Np1DNumberArray
        The scaled norms of the total linear momentum, one value
        per frame.
    """
    reader = TrajectoryReader(trajectory_files, md_format=md_format)

    momentum = Momentum(
        reader,
        selection=selection,
        use_full_atom_info=use_full_atom_info,
        scale=scale,
    )

    data_writer = MomentumDataWriter(output, mode=mode)

    momentum_norms = momentum.run()

    data_writer.write(momentum_norms)

    return momentum_norms
