"""
API functions for vibrational analysis.
"""

from pathlib import Path

import numpy as np

from PQAnalysis.io import MoldescriptorReader, RestartFileReader, read_trajectory
from PQAnalysis.type_checking import runtime_type_checking

from .exceptions import VibrationalAnalysisError
from .vibrational_analysis import (
    calculate_from_system,
    read_hessian_file,
    write_calculate_output,
    write_extxyz_modes,
    write_normal_modes,
    write_xyz_modes,
)
from .vibrational_input_file_reader import VibrationalAnalysisInputFileReader



@runtime_type_checking
def vibrations(input_file: str) -> None:
    """
    Run vibrational analysis from an input file.
    """
    input_reader = VibrationalAnalysisInputFileReader(input_file)
    input_reader.read()

    system = _read_structure_file(
        input_reader.structure_file,
        input_reader.moldescriptor_file,
    )
    hessian = read_hessian_file(input_reader.hessian_file)
    atom_charges = _read_atom_charges(system, input_reader.moldescriptor_file)

    result = calculate_from_system(
        system,
        hessian,
        unit=input_reader.unit,
        hessian_sign=input_reader.hessian_sign,
        atom_charges=atom_charges,
    )

    write_calculate_output(result, input_reader.out_file)

    if input_reader.normal_modes_file is not None:
        write_normal_modes(result.normal_modes, input_reader.normal_modes_file)

    if input_reader.modes_prefix is not None:
        atom_names = [atom.name for atom in system.atoms]
        write_xyz_modes(
            result.normal_modes,
            system.pos,
            atom_names,
            filename=input_reader.modes_prefix,
            wavenumbers=result.wavenumbers,
            modes=input_reader.modes,
            n_frames=input_reader.modes_frames,
            amplitude=input_reader.modes_amplitude,
            temperature=input_reader.modes_temperature,
            threshold=input_reader.modes_threshold,
        )

    if input_reader.modes_file is not None:
        atom_names = [atom.name for atom in system.atoms]
        write_extxyz_modes(
            result.normal_modes,
            system.pos,
            atom_names,
            filename=input_reader.modes_file,
            wavenumbers=result.wavenumbers,
            intensities=result.intensities,
            modes=input_reader.modes,
            threshold=input_reader.modes_threshold,
        )



def _read_structure_file(
    structure_file: str,
    moldescriptor_file: str | None = None,
):
    """
    Read a restart file or single-frame XYZ structure file.
    """
    if Path(structure_file).suffix.lower() == ".xyz":
        trajectory = read_trajectory(structure_file)

        if len(trajectory) != 1:
            raise VibrationalAnalysisError(
                "XYZ structure input must contain exactly one frame."
            )

        return trajectory[0]

    return RestartFileReader(
        structure_file,
        moldescriptor_filename=moldescriptor_file,
    ).read()



def _read_atom_charges(
    system, moldescriptor_file: str | None
) -> np.ndarray | None:
    """
    Read atom charges from a moldescriptor file if one was provided.
    """
    if moldescriptor_file is None:
        return None

    residues = MoldescriptorReader(moldescriptor_file).read()

    if Path(moldescriptor_file).is_file() and not system.topology.residues:
        if len(residues) != 1:
            raise VibrationalAnalysisError(
                "XYZ input requires a moldescriptor file with exactly one molecule type."
            )

        residue = residues[0]
        if residue.n_atoms != system.n_atoms:
            raise VibrationalAnalysisError(
                "The moldescriptor molecule size must match the XYZ structure."
            )

        return np.asarray(residue.partial_charges, dtype=float)

    charges = np.zeros(system.n_atoms, dtype=float)
    residue_by_id = {residue.id: residue for residue in residues}

    for residue_id in np.unique(system.topology.residue_ids):
        if residue_id not in residue_by_id:
            raise VibrationalAnalysisError(
                f"Residue id {residue_id} has no moldescriptor entry."
            )

    for residue, atom_indices in zip(
        system.topology.residues,
        system.topology.residue_atom_indices,
    ):
        if residue.n_atoms != len(atom_indices):
            raise VibrationalAnalysisError(
                "The moldescriptor residue size does not match the structure."
            )

        charges[atom_indices] = residue.partial_charges

    return charges
