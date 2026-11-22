"""
Public functions of the ParmEd-based force-field import.
"""

from numbers import Real
from pathlib import Path

from beartype.typing import Any, Sequence

from PQAnalysis.type_checking import runtime_type_checking

from .exceptions import Parm2PQError
from .parmed_converter import convert_structure
from .pq_setup import PQForceFieldSetup



def import_parmed() -> Any:
    """
    Imports ParmEd, which is an optional dependency.

    Returns
    -------
    module
        The parmed module.

    Raises
    ------
    ImportError
        If ParmEd is not installed.
    """
    try:
        import parmed  # pylint: disable=import-outside-toplevel
    except ImportError as error:
        raise ImportError(
            "parm2pq needs ParmEd. Install it with "
            "'pip install PQAnalysis[ff]' or 'pip install parmed'."
        ) from error
    return parmed



@runtime_type_checking
def load_structure(
    topology_file: str | Path,
    coordinate_file: str | Path | None = None,
    parameter_files: Sequence[str | Path] = (),
) -> Any:
    """
    Loads a parametrized structure with ParmEd.

    Parameters
    ----------
    topology_file : str | Path
        AMBER prmtop, GROMACS top or CHARMM psf file.
    coordinate_file : str | Path | None, optional
        Coordinates (inpcrd/rst7, gro, crd, pdb, ...), by default None.
    parameter_files : Sequence[str | Path], optional
        CHARMM parameter and stream files (prm, rtf, str), needed for a
        psf topology, by default ().

    Returns
    -------
    parmed.Structure
        The structure with parameters and, if given, coordinates.

    Raises
    ------
    Parm2PQError
        If a psf topology is given without parameter files.
    """
    parmed = import_parmed()
    topology_file = str(topology_file)

    structure = parmed.load_file(topology_file)
    is_psf = isinstance(structure, parmed.charmm.CharmmPsfFile)
    if is_psf:
        if not parameter_files:
            raise Parm2PQError(
                "A CHARMM psf topology needs its parameter files "
                "(--parameters)."
            )
        parameters = parmed.charmm.CharmmParameterSet(
            *[str(file) for file in parameter_files]
        )
        structure.load_parameters(parameters)

    if coordinate_file is not None:
        coordinates = parmed.load_file(str(coordinate_file))
        structure.coordinates = coordinates.coordinates
        if getattr(coordinates, "box", None) is not None:
            structure.box = coordinates.box
        velocities = getattr(coordinates, "velocities", None)
        if velocities is not None:
            structure.velocities = velocities

    return structure



@runtime_type_checking
def parm2pq(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    topology_file: str | Path,
    coordinate_file: str | Path | None = None,
    parameter_files: Sequence[str | Path] = (),
    output_dir: str | Path = ".",
    shake_h: bool = False,
    rcoulomb: Real = 12.0,
    overwrite: bool = False,
) -> PQForceFieldSetup:
    """
    Converts an AMBER, CHARMM or GROMACS force field into PQ setup files.

    The function writes ``moldescriptor.dat``, ``topology.top``,
    ``parameter.prm``, ``intra_nonbonded.dat``, ``start.rst`` and the
    input snippet ``pq_input.in`` into ``output_dir``, and ``mshake.dat``
    if ``shake_h`` keeps rigid water rigid with M-SHAKE.

    Parameters
    ----------
    topology_file : str | Path
        AMBER prmtop, GROMACS top or CHARMM psf file.
    coordinate_file : str | Path | None, optional
        The coordinates; required unless the topology contains them.
    parameter_files : Sequence[str | Path], optional
        CHARMM parameter and stream files for a psf topology.
    output_dir : str | Path, optional
        The output directory, by default the working directory.
    shake_h : bool, optional
        Constrain bonds to hydrogen with SHAKE, by default False. Rigid
        three-site water (with an H-H bond) gets M-SHAKE instead.
    rcoulomb : float, optional
        The Coulomb cutoff in Angstrom for the input snippet and the
        vacuum box of a structure without box, by default 12.0.
    overwrite : bool, optional
        Replace existing output files, by default False.

    Returns
    -------
    PQForceFieldSetup
        The converted setup.

    Raises
    ------
    UnsupportedForceFieldTermError
        If the force field contains a term PQ cannot evaluate.
    Parm2PQError
        If the input cannot be converted.
    FileExistsError
        If an output file exists and overwrite is False.
    """
    structure = load_structure(topology_file, coordinate_file, parameter_files)
    source = Path(topology_file).name
    if coordinate_file is not None:
        source += f" + {Path(coordinate_file).name}"

    setup = convert_structure(
        structure,
        shake_h=shake_h,
        rcoulomb=rcoulomb,
        source=source,
    )
    setup.write(output_dir, rcoulomb=rcoulomb, overwrite=overwrite)
    return setup
